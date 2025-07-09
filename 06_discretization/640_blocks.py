#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon
# venv: aa

from compas.datastructures import Mesh
from compas.geometry import Frame
from compas.geometry import Polygon
from compas.geometry import Plane
from compas.geometry import bestfit_frame_numpy
from compas.geometry import Line
from compas.itertools import pairwise
from compas.scene import Scene
from compas_libigl.mapping import map_mesh
from compas_cgal.meshing import trimesh_dual
from compas_libigl.parametrisation import trimesh_lsc_mapping
from compas_rhino.objects import select_mesh
from compas_rhino.objects import select_points
from compas_rhino.objects import select_curves
from compas_rhino.conversions import meshobject_to_compas
from compas_rhino.conversions import pointobject_to_compas
from compas_rhino.conversions import curveobject_to_compas

# mesh from rhino
guid = select_mesh()
mesh = meshobject_to_compas(guid)

# points form rhino
guids = select_points()
points = []
for guid in guids:
    point = pointobject_to_compas(guid)
    points.append(point)

# polygons from rhino
guids = select_curves()
polygons = []

for guid in guids:
    curve_points = curveobject_to_compas(guid).points
    polygon = Polygon(curve_points[:-1])
    polygons.append(polygon)
pattern = Mesh.from_polygons(polygons)

# =============================================================================
# Remesh and Flatten
# =============================================================================

fixed_vertices = []
for point in points:
    for v in mesh.vertices():
        mesh_point = mesh.vertex_point(v)
        if mesh_point.distance_to_point(point) < 0.01:
            fixed_vertices.append(v)
            break

V, F, DV, DF = trimesh_dual(
    mesh.to_vertices_and_faces(True),
    length_factor=1.0,
    number_of_iterations=10,
    angle_radians=0.9,
    scale_factor=1.0,
    fixed_vertices=fixed_vertices,
)

trimesh = Mesh.from_vertices_and_faces(V, F)

# ==============================================================================
# Least-squares conformal map - to see where the pattern is mapped.
# ==============================================================================

uv = trimesh_lsc_mapping((V, F))
flat_mesh = trimesh.copy()
for i in range(trimesh.number_of_vertices()):
    flat_mesh.vertex_attributes(i, "xyz", [uv[i][0], uv[i][1], 0])

# ==============================================================================
# Select Pattern and Map
# ==============================================================================

mv, mf, mn, mb, mg = map_mesh(
    (V, F),
    pattern.to_vertices_and_faces(),
    clip_boundaries=False,
    simplify_borders=False,
)
pattern = Mesh.from_vertices_and_faces(mv, mf)
polygons = pattern.to_polygons()

# ==============================================================================
# Create Blocks
# ==============================================================================


def pattern_inverse_height_thickness(pattern: Mesh, tmin=None, tmax=None):
    x: list[float] = pattern.vertices_attribute(name="x")  # type: ignore
    xmin = min(x)
    xmax = max(x)
    dx = xmax - xmin

    y: list[float] = pattern.vertices_attribute(name="y")  # type: ignore
    ymin = min(y)
    ymax = max(y)
    dy = ymax - ymin

    d = (dx**2 + dy**2) ** 0.5

    tmin = tmin or 3 * d / 1000
    tmax = tmax or 50 * d / 1000

    pattern.update_default_vertex_attributes(thickness=0)
    zvalues: list[float] = pattern.vertices_attribute(name="z")  # type: ignore
    zmin = min(zvalues)
    zmax = max(zvalues)

    for vertex in pattern.vertices():
        point = pattern.vertex_point(vertex)
        z = (point.z - zmin) / (zmax - zmin)
        thickness = (1 - z) * (tmax - tmin) + tmin
        pattern.vertex_attribute(vertex, name="thickness", value=thickness)


def pattern_idos(pattern: Mesh):
    idos: Mesh = pattern.copy()
    for vertex in idos.vertices():
        point = pattern.vertex_point(vertex)
        normal = pattern.vertex_normal(vertex)
        thickness = pattern.vertex_attribute(vertex, name="thickness")
        idos.vertex_attributes(
            vertex, names="xyz", values=point - normal * (0.5 * thickness)
        )  # type: ignore
    return idos


def pattern_face_block(pattern: Mesh, idos: Mesh, face: int) -> Mesh:
    vertices = pattern.face_vertices(face)
    normals = [pattern.vertex_normal(vertex) for vertex in vertices]
    thickness = pattern.vertices_attribute("thickness", keys=vertices)
    bottom = idos.vertices_points(vertices)
    top = [point + vector * t for point, vector, t in zip(bottom, normals, thickness)]  # type: ignore
    frame = Frame(*bestfit_frame_numpy(top))
    plane = Plane.from_frame(frame)
    flattop = []
    for a, b in zip(bottom, top):
        b = plane.intersection_with_line(Line(a, b))
        flattop.append(b)
    sides = []
    for (a, b), (aa, bb) in zip(
        pairwise(bottom + bottom[:1]), pairwise(flattop + flattop[:1])
    ):
        sides.append([a, b, bb, aa])
    polygons = [bottom[::-1], flattop] + sides
    block: Mesh = Mesh.from_polygons(polygons)
    return block


def pattern_blocks(pattern: Mesh, idos: Mesh) -> dict[int, Mesh]:
    face_block: dict[int, Mesh] = {}
    face: int
    for face in pattern.faces():  # type: ignore
        face_block[face] = pattern_face_block(pattern, idos, face)
    return face_block


pattern_inverse_height_thickness(pattern, 0.1, 0.4)
idos = pattern_idos(pattern)
blocks = pattern_blocks(pattern, idos)

# =============================================================================
# Vizualize
# =============================================================================

scene = Scene()
for key, block in blocks.items():
    print(block)
    scene.add(block)
scene.draw()

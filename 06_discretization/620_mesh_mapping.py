#! python 3
# r: compas_model, Tessagon, compas_cgal==0.9.1, compas_libigl==0.7.4
# venv: brg-csd

from compas.datastructures import Mesh
from compas.geometry import Polygon
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
    clip_boundaries=True,
    simplify_borders=False,
)
pattern = Mesh.from_vertices_and_faces(mv, mf)
polygons = pattern.to_polygons()

# =============================================================================
# Vizualize
# =============================================================================

scene = Scene()
for polygon in polygons:
    scene.add(Polygon(polygon).translated([0, 0, 1]))
scene.draw()

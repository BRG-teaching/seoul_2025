#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon
# venv: aa

from compas.datastructures import Mesh
from compas.scene import Scene
from compas_cgal.meshing import trimesh_dual
from compas_libigl.parametrisation import trimesh_lsc_mapping
from compas_rhino.objects import select_mesh
from compas_rhino.objects import select_points
from compas_rhino.conversions import meshobject_to_compas
from compas_rhino.conversions import pointobject_to_compas

guid = select_mesh()
mesh = meshobject_to_compas(guid)

guids = select_points()
points = []
for guid in guids:
    point = pointobject_to_compas(guid)
    points.append(point)

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

# =============================================================================
# Vizualize
# =============================================================================

scene = Scene()
scene.add(trimesh.translated([0, 0, 0.5]))
scene.add(flat_mesh)
scene.draw()

#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon
# venv: aa

import pathlib
import compas
from compas.datastructures import Mesh
from compas.scene import Scene

# =============================================================================
# RhinoVault mesh
# =============================================================================

mesh = compas.json_load(pathlib.Path(__file__).parent /  "ThrustDiagram.json")

points = []
for v in mesh.vertices():
    if mesh.vertex_attribute(v, "is_support"):
        points.append(mesh.vertex_point(v))

V, F = mesh.to_vertices_and_faces()
mesh = Mesh.from_vertices_and_faces(V, F)

# =============================================================================
# Vizualize
# =============================================================================

scene = Scene()
scene.add(mesh)
for p in points:
    scene.add(p)
scene.draw()

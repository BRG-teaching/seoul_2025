#! python 3
# r: compas_model, Tessagon, compas_cgal==0.9.1, compas_libigl==0.7.4
# venv: brg-csd

import pathlib
import compas
from compas.datastructures import Mesh
from compas.scene import Scene

# =============================================================================
# FormFinder Mesh
# =============================================================================

session = compas.json_load(pathlib.Path(__file__).parent / "FormFinder.json")
scene = session["scene"]
mesh = scene.find_by_name("CableMesh").mesh

points = []
for v in mesh.vertices():
    if mesh.vertex_attribute(v, "is_support"):
        points.append(mesh.vertex_point(v))

V, F = mesh.to_vertices_and_faces()
mesh = Mesh.from_vertices_and_faces(V, F)

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

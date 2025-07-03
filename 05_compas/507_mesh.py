#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon
# venv: aa

import compas
from compas.scene import Scene
from compas.datastructures import Mesh
import compas

# Create a mesh from an obj file
mesh = Mesh.from_obj(compas.get("tubemesh.obj"))
print(mesh)

# Create a scene
scene = Scene()
scene.clear_context()

# Add the mesh to the scene
scene.add(mesh)

# Draw the scene
scene.draw()

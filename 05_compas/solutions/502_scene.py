#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon
# venv: aa

from compas.geometry import Box
from compas.scene import Scene

# Create a box
box = Box(1, 2, 3)

# Create a scene and add the box to it
scene = Scene()
scene.clear_context()
scene.add(box)

# Draw the scene
scene.draw()

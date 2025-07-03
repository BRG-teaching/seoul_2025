#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon
# venv: aa

from compas.scene import Scene
from compas.geometry import Box
from compas.geometry import Frame

scene = Scene()
scene.clear_context()

# Create a frame
frame = Frame([0, 0, 0], [1, 1, 1], [0, 1, 0])
scene.add(frame)

# Create a box
box = Box(1, 2, 3, frame=frame)
scene.add(box)

# Draw the scene
scene.draw()

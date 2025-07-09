#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon
# venv: aa

from compas.colors import Color
from compas.geometry import Box
from compas.scene import Scene

# Create a box
box = Box(1)

# Convert the box to different representations
poly = box.to_polyhedron()
mesh = box.to_mesh()
brep = box.to_brep()

# Print the box
print(box)
print("=" * 100)

# Print the polyhedron
print(poly)
print("=" * 100)

# Print the mesh
print(mesh)
print("=" * 100)

# Print the brep
print(brep)
print("=" * 100)

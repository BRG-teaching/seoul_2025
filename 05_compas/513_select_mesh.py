#! python 3
# r: compas_model, Tessagon, compas_cgal==0.9.1, compas_libigl==0.7.4
# venv: brg-csd

import compas
from compas_rhino.objects import select_mesh
from compas_rhino.conversions import meshobject_to_compas
from compas.datastructures import Mesh
from compas.geometry import Polygon
from compas.scene import Scene

# mesh from rhino
guid = select_mesh()
mesh = meshobject_to_compas(guid)
print(mesh)

# =============================================================================
# Vizualize
# =============================================================================

scene = Scene()
scene.add(mesh.translated([0,0,1]))
scene.draw()

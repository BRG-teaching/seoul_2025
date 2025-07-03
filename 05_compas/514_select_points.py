#! python 3
# r: compas_model, Tessagon, compas_cgal==0.9.1, compas_libigl==0.7.4
# venv: brg-csd

from compas.datastructures import Mesh
from compas.geometry import Polygon
from compas.scene import Scene
from compas_rhino.objects import select_points
from compas_rhino.conversions import pointobject_to_compas

# points form rhino
guids = select_points()
points = []
for guid in guids:
    point = pointobject_to_compas(guid)
    points.append(point)

# =============================================================================
# Vizualize
# =============================================================================

scene = Scene()
for p in points:
    scene.add(p.translated([0,0,1]))
scene.draw()

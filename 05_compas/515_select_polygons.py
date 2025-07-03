#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon
# venv: aa

from compas.geometry import Polygon
from compas.scene import Scene
from compas_rhino.objects import select_curves
from compas_rhino.conversions import curveobject_to_compas

# polygons from rhino
guids = select_curves()
polygons = []

for guid in guids:
    curve_points = curveobject_to_compas(guid).points
    polygon = Polygon(curve_points[:-1])
    polygons.append(polygon)

# =============================================================================
# Vizualize
# =============================================================================

scene = Scene()
for polygon in polygons:
    scene.add(polygon.translated([0,0,1]))
scene.draw()

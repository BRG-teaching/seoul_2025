#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon
# venv: aa

from compas.geometry import Box
from compas import json_dump
import pathlib

# Create a box
box = Box(1)

# Define the file path to save the box to JSON
filepath = pathlib.Path(__file__).parent / "box.json"

# Save the box to JSON file
json_dump(box, filepath)

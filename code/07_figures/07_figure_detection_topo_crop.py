# Crops the topo_{I,H,F,E}.png files written by code/07_figures/07_figure_detection_topo.R
# and saves them as detection_topo_scenario{1,2,3,4}.png (I=1, H=2, F=3, E=4).
# The crop removes 245 px on the left and 246 px on the right (1653 -> 1162 px wide).
import os
import sys
from PIL import Image

# Run from the project root (PROJECT_ROOT overrides it).
project = os.path.abspath(os.environ.get("PROJECT_ROOT", os.getcwd()))
src = os.path.join(project, "results", "heatmaps")
out = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("THESIS_DIR", os.path.join(os.pardir, "thesis")) + "/diagnostics"
scenario = {"I": 1, "H": 2, "F": 3, "E": 4}
box = (245, 0, 245 + 1162, 2125)

for letter, k in scenario.items():
    im = Image.open(os.path.join(src, f"topo_{letter}.png"))
    assert im.size == (1653, 2125), im.size
    im.crop(box).save(os.path.join(out, f"detection_topo_scenario{k}.png"))
    print("saved scenario", k)

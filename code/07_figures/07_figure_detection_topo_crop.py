# Thesis Figures C.1 to C.4 (fig:detection_topo_1..4); run after 07_figure_detection_topo.R.
# Crops the topo_{1,2,3,4}.png files written by code/07_figures/07_figure_detection_topo.R
# and saves them as detection_topo_scenario{1,2,3,4}.png (scenario numbers 1 to 4).
# The crop removes 245 px on the left and 246 px on the right (1653 -> 1162 px wide).
import os
import sys
from PIL import Image

# Run from the project root (PROJECT_ROOT overrides it).
project = os.path.abspath(os.environ.get("PROJECT_ROOT", os.getcwd()))
# FIG_OUT, if set, holds the topo_*.png files and receives the cropped images.
fig_out = os.environ.get("FIG_OUT", "")
src = fig_out or os.path.join(project, "results", "heatmaps")
out = sys.argv[1] if len(sys.argv) > 1 else (fig_out or os.path.join(os.environ.get("THESIS_DIR", os.path.join(project, "results")), "diagnostics"))
os.makedirs(out, exist_ok=True)
scenario = {"1": 1, "2": 2, "3": 3, "4": 4}
box = (245, 0, 245 + 1162, 2125)

for key, k in scenario.items():
    im = Image.open(os.path.join(src, f"topo_{key}.png"))
    assert im.size == (1653, 2125), im.size
    im.crop(box).save(os.path.join(out, f"detection_topo_scenario{k}.png"))
    print("saved scenario", k)

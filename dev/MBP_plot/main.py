import os
import json
from pymatgen.core import Element
from pymatgen.core import periodic_table as pt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.metalig.mol import Ligand
import numpy as np
from matplotlib import rcParams
rcParams['svg.fonttype'] = 'none'

# Path to cache file
cache_file = "mbp_cache.json"

# Try to load precomputed mbp_dict_total
if os.path.exists(cache_file):
    print("Loading cached data from", cache_file)
    with open(cache_file, "r") as f:
        mbp_dict_total = json.load(f)
else:
    print("Cache not found. Parsing LigandDB...")
    ligand_db = LigandDB.from_json(n_max=560000)

    # Collect total metal binding counts
    mbp_dict_total = {}
    for lig_name, lig_data in ligand_db.db.items():
        lig_data: Ligand
        mbp_dict = lig_data.metal_counts
        for metal, count in mbp_dict.items():
            mbp_dict_total[metal] = mbp_dict_total.get(metal, 0) + 1

    # Save to cache
    with open(cache_file, "w") as f:
        json.dump(mbp_dict_total, f)
    print("Cached data written to", cache_file)

# Load periodic table data
with open(pt.ptable_json.name) as f:
    data = json.load(f)

# Prepare data for plotting
x, y, text, color = [], [], [], []

for el in data:
    try:
        element = Element(el)
        group = element.group
        period = element.row
        val = mbp_dict_total.get(element.symbol, 0)

        # Lanthanides (Period 8) and Actinides (Period 9)
        if group == 3 and period == 6:
            group = element.number - Element("La").number + Element("La").group
            period = 8
        elif group == 3 and period == 7:
            group = element.number - Element("Ac").number + Element("Ac").group
            period = 9

        x.append(group)
        y.append(period)
        text.append(element.symbol)
        color.append(val)
    except:
        continue

# Log scale for color mapping
color_log = [np.log10(c) if c > 0 else np.nan for c in color]

# Create figure
fig, ax = plt.subplots(figsize=(12, 8))
ax.set_facecolor("white")

# Define custom colormap
colors = ["#AFCEFF", "#A6ECFE", "#9FFDDD" ,"#9FFDA4", "#D7FBA3", "#FCEBA3", "#FBC6A3", "#FFA8A6"]  # From smallest to largest
bounds = [10**1.375, 10**1.75, 10**2.125, 10**2.5, 10**2.875, 10**3.25, 10**3.625, 10**4]
log_bounds = np.log10(bounds)
norm_values = (log_bounds - log_bounds.min()) / (log_bounds.max() - log_bounds.min())
cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", list(zip(norm_values, colors)))

# Normalization
norm = mcolors.LogNorm(vmin=1, vmax=10 ** 4)


# Helper: get brightness for text color
def get_brightness(rgb):
    r, g, b = rgb[:3]
    return 0.299 * r + 0.587 * g + 0.114 * b


# Plot circles
for xi, yi, label, raw_val in zip(x, y, text, color):
    facecolor = cmap(norm(raw_val)) if raw_val > 0 else (0.7, 0.7, 0.7, 1)
    circ = plt.Circle((xi, yi), 0.4, facecolor=facecolor, edgecolor='none')
    ax.add_patch(circ)

    brightness = get_brightness(facecolor)
    text_color = 'white' if brightness < 0.5 else 'black'
    ax.text(xi, yi, label, ha='center', va='center', fontsize=10, color=text_color)

ax.set_xlim(0.5, 18.5)
ax.set_ylim(0.5, 10.5)
ax.set_xticks(range(1, 19))
ax.set_yticks(range(1, 10))
ax.set_xlabel("Group")
ax.set_ylabel("Period")
ax.set_title("Periodic Table Heatmap", fontsize=16)
ax.set_aspect('equal', adjustable='datalim')  # 🔥 Fix squished circles
ax.invert_yaxis()

# Colorbar
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.03, pad=0.04)
cbar.set_label("Value (log scale)")

# Explicit major ticks: 10¹, 10², 10³, 10⁴
cbar_ticks = [10 ** 1, 10 ** 2, 10 ** 3, 10 ** 4]
cbar.ax.set_yticks(cbar_ticks)
cbar.ax.set_yticklabels([r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])
cbar.ax.tick_params(which='both', length=5)

# Save figure
plt.savefig("periodic_table_matplotlib.svg")
plt.show()

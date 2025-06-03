from pymatgen.core import Element
from pymatgen.core import periodic_table as pt
import math
import json
import plotly.graph_objects as go
import plotly.io as pio
from DARTassembler.src.ligand_extraction.DataBase import RCA_Ligand, LigandDB
import numpy as np

ligand_db = LigandDB.load_from_json(n_max=500000)

mbp_dict_total = {}
for lig_name, lig_data in ligand_db.db.items():
    lig_data: RCA_Ligand
    mbp_dict = lig_data.count_metals

    for metal, count in mbp_dict.items():
        if metal not in mbp_dict_total:
            mbp_dict_total[metal] = 0

        # This is the number of UNIQUE metals a ligand can bind to. If ligand A coordinates to Ru 30 times this is still only counted as 1
        mbp_dict_total[metal] += 1
        # This is the number of ligands that can bind to the metal
        # mbp_dict_total[metal] += count

# Load periodic table data from pymatgen's internal JSON
with open(pt.ptable_json.name) as f:
    data = json.load(f)

# Build plotly-compatible data
x = []  # group (column)
y = []  # period (row)
text = []
color = []

for el in data:
    # print(Element(el).symbol, Element(el).number, Element(el).group, Element(el).row)

    try:
        element = Element(el)
        group = element.group
        period = element.row
        val = mbp_dict_total.get(element.symbol, 0)  # Get the value from the dictionary

        # here we have lanthanides
        if group == 3 and period == 6:
            group = element.number - Element("La").number + Element("La").group
            period = 8

        # here we have actinides
        elif group == 3 and period == 7:
            group = element.number - Element("Ac").number + Element("Ac").group
            period = 9

        x.append(group)
        y.append(period)
        text.append(element.symbol)
        color.append(val)
    except:
        continue  # skip if group/row not defined



div = 10  # number of divisions per decade
decades = [10 ** i for i in range(0, 4)]  # e.g. [1, 10, 100, 1000]

tick_val_base_10 = [0]

for i in range(len(decades) - 1):
    start = decades[i]
    end = decades[i + 1]
    vals = [i for i in np.arange(start, end + 1, end / div)]
    tick_val_base_10.extend(vals)

tick_val_base_10.append(decades[-1])  # add final point (e.g. 1000)

# Now compute log values for plot
tick_vals = [np.log10(v) for v in tick_val_base_10]

tick_text = []
debug = False
for v in tick_val_base_10:
    if debug:
        tick_text.append(v)
        continue
    if v == 0:
        tick_text.append("0")
    else:
        order = int(math.log10(v))
        if (v % (2 * 10**order) == 0) or (int(v / (10**order) == 1)) and (v != 1):
            tick_text.append(str(v))
        else:
            tick_text.append(" ")
print(tick_text)

custom_colorscale = [
    [0.0, "rgb(131, 186, 215)"],
    [0.25, "rgb(168, 229, 146)"],
    [0.5, "rgb(238, 225, 54)"],
    [0.75, "rgb(244, 175, 89)"],
    [1.0, "rgb(188, 71, 75)"]
]

# Separate zero and non-zero data
x_zero = [x[i] for i in range(len(color)) if color[i] == 0]
y_zero = [y[i] for i in range(len(color)) if color[i] == 0]
text_zero = [text[i] for i in range(len(color)) if color[i] == 0]

x_nonzero = [x[i] for i in range(len(color)) if color[i] > 0]
y_nonzero = [y[i] for i in range(len(color)) if color[i] > 0]
text_nonzero = [text[i] for i in range(len(color)) if color[i] > 0]
color_log_nonzero = [np.log10(color[i]) for i in range(len(color)) if color[i] > 0]

fig = go.Figure()

# Zero-value elements (fixed grey)
fig.add_trace(go.Scatter(
    x=x_zero,
    y=y_zero,
    mode='markers+text',
    marker=dict(size=40, color='rgb(150,150,150)'),
    text=text_zero,
    textfont=dict(color='white'),
    hovertemplate="<b>%{text}</b><br>Group: %{x}<br>Period: %{y}<extra></extra>"
))

# Non-zero elements (log scale)
fig.add_trace(go.Scatter(
    x=x_nonzero,
    y=y_nonzero,
    mode='markers+text',
    marker=dict(
        size=40,
        color=color_log_nonzero,
        colorscale=custom_colorscale,
        cmin=0,
        cmax=3.5,
        colorbar=dict(title="log10(Value)",
                      tickvals=tick_vals,
                      ticktext=tick_text,
                      ticks="outside",  # style of tick marks
                      ticklen=5,  # length of ticks (optional)
                      tickwidth=1
                      ),
        showscale=True
    ),
    text=text_nonzero,
    textfont=dict(color='white'),
    hovertemplate="<b>%{text}</b><br>Group: %{x}<br>Period: %{y}<extra></extra>"
))

fig.update_layout(
    title="Periodic Table Heatmap",
    xaxis=dict(title="Group", dtick=1),
    yaxis=dict(title="Period", dtick=1, autorange="reversed"),
    plot_bgcolor='white'
)

pio.write_image(fig, "periodic_table.svg", format="svg")
fig.show()

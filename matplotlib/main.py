import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Ensures reproducibility of random numbers

# Here we specify how many unique donor groups per denticity

mono = 3
bi = 4
tri = 4
tetra = 4
penta = 4

dic_2 = {"1": {"C": 5313,
               "N": 4906,
               "P": 2826},

         "2": {"N-N": 8011,
               "N-O": 3540,
               "N-C": 2332,
               "P-P": 1717},
         "3": {"N-N-N": 3526,
               "N-N-O": 1579,
               "N-O-O": 1243,
               "C-N-N": 773},
         "4": {"N-N-N-N": 3655,
               "N-N-N-O": 2559,
               "N-N-O-O": 714,
               "C-C-C-C": 437},
         "5": {"N-N-N-N-N": 482,
               "N-N-N-N-O": 290,
               "N-N-N-O-O": 278,
               "C-C-C-C-C": 170}
         }

dic_3 = {1: {'C': 5313, 'N': 4906, 'P': 2826, 'O': 1434, 'S': 1120, 'others': 858, 'Si': 242, 'Se': 145},
         2: {'N-N': 8011, 'others': 4997, 'N-O': 3540, 'C-N': 2332, 'P-P': 1717, 'C-C': 1672, 'O-O': 1620, 'S-S': 1200},
         3: {'others': 3792, 'N-N-N': 3526, 'N-N-O': 1579, 'N-O-O': 1243, 'C-N-N': 773, 'C-C-C': 623, 'N-N-S': 517, 'N-O-S': 450},
         4: {'N-N-N-N': 3655, 'others': 1722, 'N-N-O-O': 1625, 'C-C-C-C': 714, 'N-N-N-O': 437, 'N-N-S-S': 419, 'O-O-O-O': 215, 'C-C-N-N': 203},
         5: {'others': 526, 'N-N-N-N-N': 482, 'N-N-N-O-O': 290, 'C-C-C-C-C': 278, 'N-N-N-N-O': 170, 'N-N-O-O-O': 77, 'N-N-N-S-S': 60, 'N-N-N-N-S': 36},
         6: {'C-C-C-C-C-C': 1027, 'N-N-N-N-N-N': 651, 'others': 334, 'N-N-N-N-O-O': 222, 'N-N-N-O-O-O': 84, 'O-O-O-O-O-O': 80, 'N-N-O-O-O-O': 69, 'N-N-O-O-S-S': 37},
         7: {'others': 75, 'C-C-C-C-C-C-P': 66, 'N-N-N-N-O-O-O': 52, 'C-C-C-C-C-C-C': 50, 'N-N-N-N-N-N-N': 30, 'C-C-C-C-C-C-N': 28, 'C-C-C-C-C-C-S': 15, 'C-C-C-C-C-P-P': 10},
         8: {'others': 55, 'C-C-C-C-C-C-C-C': 45, 'N-N-N-N-O-O-O-O': 35, 'C-C-C-C-C-C-N-N': 28, 'N-N-N-N-N-N-N-N': 19, 'O-O-O-O-O-O-O-O': 11, 'C-C-C-C-C-C-P-P': 11, 'N-N-O-O-O-O-O-O': 6},
         9: {'others': 12, 'N-N-N-N-N-N-O-O-O': 8, 'N-N-N-O-O-O-O-O-O': 6, 'C-C-C-C-C-C-O-O-O': 3, 'O-O-O-O-O-O-O-O-O': 3, 'C-C-C-C-C-C-N-N-N': 3, 'C-C-C-C-C-C-Ga-Ga-Ga': 2, 'C-C-C-C-C-C-C-N-N': 1},
         10: {'C-C-C-C-C-C-C-C-C-C': 7, 'others': 5, 'C-C-C-C-C-C-C-O-O-O': 3, 'C-C-C-C-C-C-C-C-C-Ge': 3, 'N-N-N-N-O-O-O-O-O-O': 2, 'C-C-C-C-C-C-O-O-O-O': 2, 'Ge-Ge-Ge-Ge-Ge-Ge-Ge-Ge-Ge-Ge': 1, 'Sn-Sn-Sn-Sn-Sn-Sn-Sn-Sn-Sn-Sn': 1},
         11: {'C-C-C-C-C-C-C-C-C-C-Si': 2, 'C-C-C-C-C-C-C-C-O-O-O': 2, 'C-C-C-C-C-C-C-C-C-C-C': 1, 'Pb-Pb-Pb-Pb-Pb-Pb-Pb-Pb-Pb-Pb-Pb': 1, 'Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi': 1},
         12: {'C-C-C-C-C-C-C-C-C-C-C-C': 19, 'B-B-C-C-C-C-C-C-C-C-C-C': 4, 'others': 3, 'Pb-Pb-Pb-Pb-Pb-Pb-Pb-Pb-Pb-Pb-Pb-Pb': 2, 'Ge-Ge-Ge-Ge-Ge-Ge-Ge-Ge-Ge-Ge-Ge-Ge': 2, 'C-C-C-C-C-C-C-C-C-C-N-N': 1, 'Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi': 1, 'Sn-Sn-Sn-Sn-Sn-Sn-Sn-Sn-Sn-Sn-Sn-Sn': 1},
         13: {'B-B-C-C-C-C-C-C-C-C-C-C-O': 1, 'C-C-C-C-C-C-C-C-C-C-C-N-N': 1, 'Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Sn-Sn-Sn-Sn': 1, 'Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Pb-Pb-Pb': 1, 'Bi-Bi-Bi-Bi-Bi-Bi-Pb-Pb-Pb-Pb-Pb-Pb-Pb': 1},
         14: {'Sb-Sb-Sb-Sb-Sb-Sb-Sb-Sb-Sn-Sn-Sn-Sn-Sn-Sn': 2, 'Sb-Sb-Sb-Sb-Sb-Sb-Sn-Sn-Sn-Sn-Sn-Sn-Sn-Sn': 2, 'C-C-C-C-C-C-C-C-C-C-C-C-N-N': 1, 'Bi-Bi-Bi-Bi-Bi-Bi-Bi-Sn-Sn-Sn-Sn-Sn-Sn-Sn': 1, 'Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Sn-Sn-Sn-Sn-Sn-Sn': 1, 'Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Bi-Pb-Pb-Pb-Pb-Pb': 1, 'C-C-C-C-C-C-C-C-C-C-C-C-C-C': 1},
         16: {'C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C': 1}}

dic_4 = {1: {'C': 5313, 'N': 4906, 'P': 2826, 'O': 1434, 'S': 1120, 'others': 858},
         2: {'N-N': 8011, 'others': 4997, 'N-O': 3540, 'C-N': 2332, 'P-P': 1717, 'C-C': 1672, 'O-O': 1620, 'S-S': 1200},
         3: {'others': 3792, 'N-N-N': 3526, 'N-N-O': 1579, 'N-O-O': 1243, 'C-N-N': 773, 'C-C-C': 623, 'N-N-S': 517, 'N-O-S': 450},
         4: {'N-N-N-N': 3655, 'others': 1722, 'N-N-O-O': 1625, 'C-C-C-C': 714, 'N-N-N-O': 437, 'N-N-S-S': 419, 'O-O-O-O': 215, 'C-C-N-N': 203},
         5: {'others': 526, 'N-N-N-N-N': 482, 'N-N-N-O-O': 290, 'C-C-C-C-C': 278, 'N-N-N-N-O': 170},
         }
DICT = dic_4
# The number of empty bars between the groups
PAD = 1


colour_1 = "#EF6461"
colour_2 = "#61A0AF"
colour_3 = "#074F57"
colour_4 = "#F6D8AE"
colour_5 = "#FF9F1C"
colour_else ="#6C9A8B"
rows = []
COLORS = []
for dent in DICT.keys():
    for key, value in DICT[dent].items():
        tmp_dic = {"name":key,
                   "value":value,
                   "group":dent}
        rows.append(tmp_dic)
        if dent == "1" or dent == 1:
            COLORS.append(colour_1)
        elif dent == "2" or dent == 2:
            COLORS.append(colour_2)
        elif dent == "3" or dent == 3:
            COLORS.append(colour_3)
        elif dent == "4" or dent == 4:
            COLORS.append(colour_4)
        elif dent == "5" or dent == 5:
            COLORS.append(colour_5)
        else:
            COLORS.append(colour_else)
df = pd.DataFrame(rows)

# Show 3 first rows
print(df.head(3))


# helper function that given the angle at which the bar is positioned and the offset used in the barchart, determines the rotation and alignment of the labels.
def get_label_rotation(angle, offset):
    # Rotation must be specified in degrees :(
    rotation = np.rad2deg(angle + offset)
    if angle <= np.pi:
        alignment = "right"
        rotation = rotation + 180
    else:
        alignment = "left"
    return rotation, alignment


def add_labels(angles, values, labels, offset, ax):
    # This is the space between the end of the bar and the label
    padding = 500

    # Iterate over angles, values, and labels, to add all of them.
    for angle, value, label, in zip(angles, values, labels):
        angle = angle

        # Obtain text rotation and alignment
        rotation, alignment = get_label_rotation(angle, offset)

        # And finally add the text
        ax.text(
            x=angle,
            y=value + padding,
            s=label,
            ha=alignment,
            va="center",
            rotation=rotation,
            rotation_mode="anchor",
            fontsize=8,
        )


OFFSET = np.pi / 2

# Reorder the dataframe
df_sorted = (
    df
    .groupby(["group"])
    .apply(lambda x: x.sort_values(["value"], ascending=False))
    .reset_index(drop=True)
)

#################################################################################################################################################################


# All this part is like the code above
VALUES = df_sorted["value"].values
LABELS = df_sorted["name"].values
GROUP = df_sorted["group"].values

ANGLES_N = len(VALUES) + PAD * len(np.unique(GROUP))
ANGLES = np.linspace(0, 2 * np.pi, num=ANGLES_N, endpoint=False)
WIDTH = (2 * np.pi) / len(ANGLES)

offset = 0
IDXS = []
GROUPS_SIZE = []
for key, value in DICT.items():
    GROUPS_SIZE.append(len(value))
for size in GROUPS_SIZE:
    IDXS += list(range(offset + PAD, offset + size + PAD))
    offset += size + PAD

fig, ax = plt.subplots(figsize=(20, 10), subplot_kw={"projection": "polar"})
ax.set_axisbelow(True)
ax.grid(color='gray', linestyle='dashed')
ax.set_theta_offset(OFFSET)
ax.set_ylim(-5000, 8500)
ax.set_frame_on(False)
ax.xaxis.grid(False)
ax.yaxis.grid(True)
ax.set_xticks([])
ax.set_yticks([100, 2000, 4000, 6000, 8000, 10000])


ax.bar(
    ANGLES[IDXS], VALUES, width=WIDTH, color=COLORS,
    edgecolor="white", linewidth=2
)

add_labels(ANGLES[IDXS], VALUES, LABELS, OFFSET, ax)

# Extra customization below here --------------------

# This iterates over the sizes of the groups adding reference
# lines and annotations.

offset = 0
for group, size in zip(DICT.keys(), GROUPS_SIZE):
    # Add line below bars
    x1 = np.linspace(ANGLES[offset + PAD], ANGLES[offset + size + PAD - 1], num=sum(GROUPS_SIZE))
    # ax.plot(x1, [-5] * 50, color="#333333")

    # Add text to indicate group
    """
    ax.text(
        np.mean(x1), -500, group, color="#333333", fontsize=14,
        fontweight="bold", ha="center", va="center"
    )
    """

    # Add reference lines at 20, 40, 60, and 80
    x2 = np.linspace(ANGLES[offset], ANGLES[offset + PAD - 1], num=sum(GROUPS_SIZE))
    x2 = np.linspace(ANGLES[offset], ANGLES[offset + PAD - 1], num=sum(GROUPS_SIZE))
    ax.plot(x2, [0] * (sum(GROUPS_SIZE)), color="#bebebe", lw=0.8)
    ax.plot(x2, [2500] * (sum(GROUPS_SIZE)), color="#bebebe", lw=0.8)
    ax.plot(x2, [5000] * (sum(GROUPS_SIZE)), color="#bebebe", lw=0.8)
    ax.plot(x2, [7500] * (sum(GROUPS_SIZE)), color="#bebebe", lw=0.8)

    offset += size + PAD
#plt.show()
plt.savefig("radial_graph.svg")

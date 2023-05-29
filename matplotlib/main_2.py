import matplotlib.pyplot as plt
import numpy as np

dic_3 = {1: {'C': 5313, 'N': 4906, 'P': 2826, 'O': 1434, 'S': 1120, 'others': 858, 'Si': 242, 'Se': 145},
         2: {'N-N': 8011, 'others': 4997, 'N-O': 3540, 'C-N': 2332, 'P-P': 1717, 'C-C': 1672, 'O-O': 1620, 'S-S': 1200},
         3: {'others': 3792, 'N-N-N': 3526, 'N-N-O': 1579, 'N-O-O': 1243, 'C-N-N': 773, 'C-C-C': 623, 'N-N-S': 517, 'N-O-S': 450},
         4: {'N-N-N-N': 3655, 'others': 1722, 'N-N-O-O': 1625, 'C-C-C-C': 714, 'N-N-N-O': 437, 'N-N-S-S': 419, 'O-O-O-O': 215, 'C-C-N-N': 203},
         5: {'others': 526, 'N-N-N-N-N': 482, 'N-N-N-O-O': 290, 'C-C-C-C-C': 278, 'N-N-N-N-O': 170, 'N-N-O-O-O': 77, 'N-N-N-S-S': 60, 'N-N-N-N-S': 36},
         }

dic_color_hex = {1: "#EF6461",
                 2: "#61A0AE",
                 3: "#075057",
                 4: "#F6D8AE",
                 5: "#FFA01C"}
dic_color_rgba = {1: (239.0, 100.0, 97.0, 1.0),
                  2: (97, 160, 174, 1.0),
                  3: (7, 80, 87, 1.0),
                  4: (246, 216, 174, 1.0),
                  5: (255, 160, 28, 1.0)}

width = 1
multiplier = 0

fig, ax = plt.subplots()
tic_list = []
number_of_bars = 0
for denticity, distribution in dic_3.items():
    num_entries = len(distribution)
    try:
        colour = dic_color_hex[denticity]
    except:
        colour = "#E9EFEF"
    for donor, amount in distribution.items():
        offset = width * multiplier
        rects = ax.bar(1 + offset, amount, width, color=colour, edgecolor="white")
        ax.annotate('{}'.format(('\n'.join(donor).replace("-", "|"))),
                    xy=(1 + offset, amount),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    fontsize=7,
                    color="dimgrey",
                    ha='center', va='bottom')
        multiplier += 1
        number_of_bars += 1

    tic_list.append(((number_of_bars * width) - (num_entries * width)) + ((num_entries * width) / 2) + (0.5 * width))

ax.set_xticks(tic_list, dic_3.keys())
plt.ylim(0, 10000)
#plt.savefig("test.svg")
plt.show()

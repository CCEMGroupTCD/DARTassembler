# Import libraries
import pandas as pd
import matplotlib.pyplot as plt

# Create dictionary
plot_dict = {'A': {'Apples': 3, 'Bananas': 5, 'Oranges': 6, 'Kiwis': 9}, 'B': {'Apples': 1, 'Bananas': 9, 'Oranges': 3, 'Kiwis': 1}, 'C': {'Apples': 6, 'Bananas': 9, 'Oranges': 3, 'Kiwis': 3}}

dic_4 = {1: {'C': 5313, 'N': 4906, 'P': 2826, 'O': 1434, 'S': 1120, 'others_1': 858},
         2: {'N-N': 8011, 'others_2': 4997, 'N-O': 3540, 'C-N': 2332, 'P-P': 1717, 'C-C': 1672, 'O-O': 1620, 'S-S': 1200},
         3: {'others_3': 3792, 'N-N-N': 3526, 'N-N-O': 1579, 'N-O-O': 1243, 'C-N-N': 773, 'C-C-C': 623, 'N-N-S': 517, 'N-O-S': 450},
         4: {'N-N-N-N': 3655, 'others_4': 1722, 'N-N-O-O': 1625, 'C-C-C-C': 714, 'N-N-N-O': 437, 'N-N-S-S': 419, 'O-O-O-O': 215, 'C-C-N-N': 203},
         5: {'others_5': 526, 'N-N-N-N-N': 482, 'N-N-N-O-O': 290, 'C-C-C-C-C': 278, 'N-N-N-N-O': 170}, }
df = pd.DataFrame(dic_4)
ax = df.transpose().plot(kind="bar", stacked=True, figsize=(20,10))
ax.legend(
    loc='center right', ncol=2, title="Year of Eating", bbox_to_anchor=(1.0, 0.5))
plt.savefig("test.svg")
plt.show()






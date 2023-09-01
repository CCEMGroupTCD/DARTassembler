import json
import networkx as nx
import matplotlib.pyplot as plt

# cusgee = '/Users/timosommer/PhD/projects/RCA/projects/CreateTMC/database/tmQMg_original/graphs/CUSGEE.gml'
# G = nx.read_gml(cusgee)
# nx.draw_networkx(G, node_size=100, with_labels=True,
#                  labels={node: G.nodes[node]["node_label"] for node in G.nodes})
# plt.show()

# with open('/Users/timosommer/PhD/projects/RCA/projects/CreateTMC/data/tmQMG_Jsons_original/tmQMG.json', 'r') as file:
#     d = json.load(file)
# c = d['CUSGEE']
#
with open('/Users/timosommer/PhD/projects/RCA/projects/CreateTMC/data/tmQMG_Jsons_original/test.json', 'r') as file:
    d1= json.load(file)
c1 = d1['CUSGEE']
# for lig in c1['ligands.py']:
#         pc = [val[1][0] for val in lig['atomic_props']['partial_charge'].values()]
#         spc = sum(pc)
#         print(spc)
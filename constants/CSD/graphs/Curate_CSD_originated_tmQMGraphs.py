"""
This script was meant to curate the graphs (extracted from the .mol2 files from the CSD) for all tmQM entries.
However, the graph nodes (even after removing all unconnected atoms and reshifting the graph nodes) dont match
the 3D coordinates from the tmQM .xyz files and thus we can - at current time - not employ the graphs for our
extraction method.
"""
import os
import json
from constants.Serverpath.serverpath import serverpath
from tqdm import tqdm
from src01.utilities_graph import graph_from_graph_dict, graph_to_dict_with_node_labels, view_graph
import networkx as nx
from pymatgen.core.periodic_table import Element as Pymatgen_Element
from src01.constants import metals_in_pse


TM_element_set = set([Pymatgen_Element.from_Z(m).symbol for m in metals_in_pse])


def split_them_into_seperate_files(d = None):
    if d is None:
        with open(f"{serverpath}/CSD/graphs/curated_tmqm_graphs_from_csd.json", "r") as f:
            d = json.load(f)

    for id, Gd in d.items():
        with open(f"{serverpath}/CSD/graphs/GraphJsonsTMQM_curated/{id}_g.json", "w+") as f:
            json.dump(Gd, f)


def get_metal_containing_subgraph(gr: nx.Graph):
    comps = list(nx.connected_components(gr))

    for c in comps:
        atom_types_in_comp = set([gr.nodes[i]["node_label"] for i in c])

        if len(atom_types_in_comp.intersection(TM_element_set)) > 0:
            return gr.subgraph(c).copy()


def main():

    with open(f"{serverpath}/CSD/graphs/CSD_Graphs.json", "r") as file:
        tmqm_graph_dict = json.load(file)

    with open(f"../identifierLists/tmqm_identifier_list.json", "r") as file:
        id_list = json.load(file)

    curated_graphs = {}

    for id, d in tqdm(tmqm_graph_dict.items()):
        if id in id_list:
            G = graph_from_graph_dict(d)
            S = get_metal_containing_subgraph(G)
            nx.relabel_nodes(S,
                             mapping={node: i for i, node in enumerate(sorted(S.nodes))},
                             copy=False)
            curated_graphs[id] = graph_to_dict_with_node_labels(S)

    with open(f"{serverpath}/CSD/graphs/curated_tmqm_graphs_from_csd.json", "w+") as f:
        json.dump(curated_graphs, f)

    split_them_into_seperate_files(curated_graphs)


if __name__ == "__main__":

    main()

    """
    import os

    a = os.listdir(f"{serverpath}/CSD/graphs/GraphJsonsTMQM_curated/")

    for p in tqdm(a):
        with open(f"{serverpath}/CSD/graphs/GraphJsonsTMQM_curated/{p}", "r") as f:
            d = json.load(f)

        G = graph_from_graph_dict(d)
        if 0 not in G.nodes:
            print(f"0 not contained in the graph {p}")
            breakpoint()
    """

    print("done")


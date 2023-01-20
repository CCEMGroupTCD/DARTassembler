"""
We want to see if the graphs of the tmQMG based RCA Molecules look as they should look like
"""
import json
import random
import networkx as nx
from tqdm import tqdm
import gc

from src01.utilities_graph import make_graph_labels_integers, graphs_are_equal
from src01.DataBase import MoleculeDB

from networkx import weisfeiler_lehman_graph_hash as graph_hash

graph_path = "../database/tmQMg/graphs"
data_path = "../data/tmQMG_Jsons"


if __name__ == "__main__":

    # Testsize = 10000
    Testsize = None

    # identifier_list = r_ind zu complexen uebersetzen
    with open(f"{data_path}/tmQMG.json", "r") as file:
        dict_ = json.load(file)

    full_csd_codes = list(dict_.keys())

    del dict_
    gc.collect()

    if Testsize is None:
        r_id_list = full_csd_codes
    else:
        r_id_list = random.sample(full_csd_codes, Testsize)

    tmQMG = MoleculeDB.from_json(json_=f'{data_path}/tmQMG.json', type_="Molecule", identifier_list=r_id_list)

    for csd_code, mol in tqdm(tmQMG.db.items(), desc="Testing Graphs"):
        G_created = mol.graph

        G_real = nx.read_gml(f"{graph_path}/{csd_code}.gml")
        G_real = nx.Graph(G_real)
        G_real = make_graph_labels_integers(G_real)

        # gr_equal = graphs_are_equal(G_created, G_real)        # takes too much cpu
        gr_equal = graph_hash(G_real, node_attr='node_label', iterations=3, digest_size=16) == graph_hash(G_created, node_attr='node_label', iterations=3, digest_size=16)

        if gr_equal is True:
            pass
        else:
            print(f"Graphs are not equal for: {csd_code}")


from src01.utilities_graph import view_graph
import networkx as nx
import pandas as pd
import json

if __name__ == "__main__":

    G = nx.read_gml("../database/tmQMg/graphs/GIDLEM.gml")

    df = pd.read_csv("../database/tmQM/data/global_mol_properties.csv")

    with open("../database/tmQMg/atomic_properties/atomic_properties.json", "r") as handle:
        data = json.load(handle)

    print("done")

"""
todo: to b deleted
"""
from src01.DataBase import LigandDB
from src02_Pre_Ass_Filtering.FilteringStage import FilterStage

# only for testing
from src02_Pre_Ass_Filtering.test import id_list

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

if __name__ == "__main__":
    """
    Is the graph reading in deterministic? -> Apparantely not
    """

    import json
    with open('../data/tmQMG_Jsons/tmQMG_Ligands_unique.json', "r") as handle:
        dict_ = json.load(handle)

    ligand_dict = dict_["CSD-DUCVIG-03-a"]

    from src01.Molecule import RCA_Ligand

    lig = RCA_Ligand.read_from_mol_dict(ligand_dict)

    G = lig.graph

    import networkx as nx
    print(nx.adjacency_matrix(G).todense())

    #print(G.nodes)

    tmQM_unique_Ligands = LigandDB.from_json(json_='../data/tmQMG_Jsons/tmQMG_Ligands_unique.json',
                                             type_="Ligand",
                                             # only for testing
                                             identifier_list=id_list
                                             )

    l = list(tmQM_unique_Ligands.db.values())[0]

    G = l.graph

    import networkx as nx
    A = nx.adjacency_matrix(G).todense()

    print(A)
    print("Done")

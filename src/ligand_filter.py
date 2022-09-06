# creating the ligand filter methods
from ASE_Molecule import *
import collections
import networkx as nx


def same_sum_formula(lig1: ASE_Ligand, lig2: ASE_Ligand):
    sum_formula_1 = [a[0] for a in lig1.xyz.coordinates.values()]
    sum_formula_2 = [a[0] for a in lig2.xyz.coordinates.values()]

    return collections.Counter(sum_formula_1) == collections.Counter(sum_formula_2)


def node_check(dict1, dict2):
    return dict1["label"] == dict2["label"]


def same_structure(lig1: ASE_Ligand, lig2: ASE_Ligand):
    g1, g2 = lig1.get_graph(), lig2.get_graph()
    return nx.is_isomorphic(g1, g2, node_match=node_check)


def remove_duplicants(ligand_list: list):
    new_ligand_list = ligand_list.copy()

    while len(ligand_list) > 1:
        lig = ligand_list.pop()

        for lig2 in ligand_list:
            if same_sum_formula(lig, lig2) is True:
                if same_structure(lig, lig2) is True and lig2 in new_ligand_list:
                    lig.view_3d()
                    lig2.view_3d()
                    input("Press Enter to continue")
                    new_ligand_list.remove(lig2)

    return new_ligand_list

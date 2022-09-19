# creating the ligand filter methods
from ASE_Molecule import *
import collections
import networkx as nx
import pickle
from tqdm import tqdm


def same_sum_formula(lig1: ASE_Ligand, lig2: ASE_Ligand):
    sum_formula_1 = [a[0] for a in lig1.xyz.coordinates.values()]
    sum_formula_2 = [a[0] for a in lig2.xyz.coordinates.values()]

    return collections.Counter(sum_formula_1) == collections.Counter(sum_formula_2)


def node_check(dict1, dict2):
    return dict1["label"] == dict2["label"]


def same_structure(lig1: ASE_Ligand, lig2: ASE_Ligand):
    g1, g2 = lig1.get_graph(), lig2.get_graph()
    return nx.is_isomorphic(g1, g2, node_match=node_check)


def remove_duplicants(ligand_list: list, denticity: int):
    new_ligand_list = ligand_list.copy()

    while len(ligand_list) > 1:
        print(f"remaining (denticity {denticity}) : {len(ligand_list)}")
        lig = ligand_list.pop()

        for lig2 in ligand_list:
            if same_sum_formula(lig, lig2) is True:
                if same_structure(lig, lig2) is True and lig2 in new_ligand_list:
                    new_ligand_list.remove(lig2)

    return new_ligand_list


def duplicant_filter(ligand_dict: dict):
    for denticity, ligand_list in ligand_dict.items():
        if denticity == 2:
            num_batches = int(len(ligand_list) / 1000)

            batch_list = [ligand_list[n*1000:(n+1)*1000] for n in range(num_batches)]
            batch_list.append(ligand_list[num_batches*100:])

            reduced_list = list()

            for _list in batch_list:
                reduced_list += remove_duplicants(ligand_list=_list, denticity=denticity)

            ligand_dict[denticity] = remove_duplicants(ligand_list=reduced_list, denticity=denticity)

    return ligand_dict


if __name__ == "__main__":
    with open("../data/ligand_dict.pickle", "rb") as handle:
        ligand_dict = pickle.load(handle)

    ligand_dict_new = duplicant_filter(ligand_dict)

    with open("../data/ligand_dict_filtered.pickle", "wb") as handle:
        pickle.dump(ligand_dict_new, handle)
    print("done")
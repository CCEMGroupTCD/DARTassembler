# creating the ligand filter methods
from RCA_Molecule import *
import collections
import networkx as nx
import pickle
from tqdm import tqdm
import gc               # we need the garbage collector


def node_check(dict1, dict2):
    return dict1["label"] == dict2["label"]


def same_structure(lig1: RCA_Ligand, lig2: RCA_Ligand):
    g1, g2 = lig1.get_graph(), lig2.get_graph()
    return nx.is_isomorphic(g1, g2, node_match=node_check)


def remove_duplicants(ligand_list: list, denticity: int):
    """
    batch_filter, very standard way
    """
    new_ligand_list = list()

    while len(ligand_list) > 1:
        lig = ligand_list.pop()
        new_ligand_list.append(lig)

        for lig2 in ligand_list:
            if lig.same_sum_formula(lig2) is True:
                if same_structure(lig, lig2) is True:
                    ligand_list.remove(lig2)
                #if same_structure(lig, lig2) is True and lig2 in new_ligand_list:
                    #new_ligand_list.remove(lig2)

    return new_ligand_list


def duplicant_filter(ligand_dict: dict):

    new_dict = {}

    for denticity, ligand_list in ligand_dict.items():
        if denticity == 2 or denticity == 5:
            size = 1000
            num_batches = int(len(ligand_list) / size)

            batch_list = [ligand_list[n*size:(n+1)*size] for n in range(num_batches)]
            batch_list.append(ligand_list[num_batches*size:])

            reduced_list = list()
            del ligand_list

            while len(batch_list) > 0:
                _list = batch_list.pop()
                reduced_list += remove_duplicants(ligand_list=_list, denticity=denticity)
                del _list
                gc.collect()

            new_dict[denticity] = remove_duplicants(ligand_list=reduced_list, denticity=denticity)
        else:
            new_dict[denticity] = ligand_list

    return new_dict


def new_duplicant_filter(ligand_dict: dict):

    new_dict = {}

    for denticity, ligand_list in ligand_dict.items():
        if denticity == 2 or denticity == 5:
            new_dict[denticity] = list(set(ligand_list))
        else:
            new_dict[denticity] = ligand_list

    return new_dict


if __name__ == "__main__":
    with open("../data/ligand_dict.pickle", "rb") as handle:
        ligand_dict = pickle.load(handle)

    ligand_dict_one = duplicant_filter(ligand_dict)

    ligand_dict_two = new_duplicant_filter(ligand_dict)

    for key in ligand_dict:
        print(f"{key} - filtered (good): {len(ligand_dict_one[key])}, filtered (approx): {len(ligand_dict_two[key])}, unfiltered: {len(ligand_dict[key])}")
        #print(f"{key} - filtered (good): {len(ligand_dict_one[key])}, unfiltered: {len(ligand_dict[key])}")
        #print(f"{key}  unfiltered: {len(ligand_dict[key])}")
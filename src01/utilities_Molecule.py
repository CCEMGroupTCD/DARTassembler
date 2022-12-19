import collections

from pymatgen.core.composition import Composition

# Removed this due to circular import error, and this was only used for type hinting.
# from src01.Molecule import RCA_Ligand


def get_all_ligands_by_graph_hashes(all_ligands: list) -> dict:
    """
    Get dictionary of graph hashes with list of ligands with this graph hash
    :param all_ligands: list of all ligands with graph hashes
    :return: dictionary of graph hash: list_of_ligands
    """
    all_hashes = list(set([lig.graph_hash for lig in all_ligands]))
    all_ligands_by_hashes = {h: [] for h in all_hashes}
    for ligand in all_ligands:
        all_ligands_by_hashes[ligand.graph_hash].append(ligand)

    return all_ligands_by_hashes


def group_list_without_hashing(ligand_list: list) -> list:
    """
    Returns a list of list with unique elements grouped together. Works without hashing, just using equity.
    :param ligand_list: list of elements
    :return: list of lists of grouped elements
    """
    groupings = {}
    counter = 0
    for lig1 in ligand_list:

        tmp_groupings = {}
        equal = False

        for i, lig_list in groupings.items():
            lig_representative = lig_list[0]
            equal = lig_representative == lig1

            if equal:

                if i in tmp_groupings:
                    tmp_groupings[i].append(lig1)
                else:
                    tmp_groupings[i] = lig1

                break

        if not equal:
            tmp_groupings[counter] = [lig1]
            counter += 1

        for i in tmp_groupings.keys():

            if i in groupings:
                groupings[i].append(tmp_groupings[i])
            else:
                groupings[i] = tmp_groupings[i]

    groupings = [group for group in groupings.values()]
    return groupings


def original_metal_ligand(ligand):
    """
    We try to find the original metal of a ligand and return None if we couldnt find any
    """

    if hasattr(ligand, "original_metal"):
        return ligand.original_metal_symbol
    elif ligand.global_props is not None:
        try:
            return ligand.global_props['metal_name']
        except KeyError:
            pass
    else:
        return None

def get_standardized_stoichiometry_from_atoms_list(atoms: list) -> str:
    c = collections.Counter(atoms)
    elements = sorted(el for el in c.keys())
    if "C" in elements:
        if 'H' in elements:
            elements = ["C", 'H'] + [el for el in elements if el not in ["C", 'H']]
        else:
            elements = ["C"] + [el for el in elements if el != "C"]
    else:
        if 'H' in elements:
            elements = ["H"] + [el for el in elements if el != "H"]

    formula = [f"{el}{(c[el]) if c[el] != 1 else ''}" for el in elements]
    return "".join(formula)

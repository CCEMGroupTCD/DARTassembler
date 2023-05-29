import collections

import networkx as nx
from rdkit import Chem

from pymatgen.core.composition import Composition

# Removed this due to circular import error, and this was only used for type hinting.
# from src01.Molecule import RCA_Ligand

unknown_rdkit_bond_orders = [0, 20, 21]

def graph_to_rdkit_mol(graph: nx.Graph, element_label: str='node_label', bond_label: str='bond_type') -> Chem.Mol:
    """
    Create an rdkit mol object from a graph. Note that the bond type must be specified in the graph under the attribute called `edge_label`.
    @param graph: input graph of the molecule
    @param element_label: element label for the node dictionary
    @param bond_label: bond type label for the edge dictionary
    @return: rdkit mol object
    """

    # create empty editable mol object
    mol = Chem.RWMol()

    # add atoms to mol and keep track of index
    node_to_idx = {}
    for idx, atom in graph.nodes(data=True):
        a = Chem.Atom(atom[element_label])
        molIdx = mol.AddAtom(a)
        node_to_idx[idx] = molIdx

    # add bonds between adjacent atoms
    for idx1, idx2, bond in graph.edges(data=True):
        bond = bond[bond_label]

        try:
            bond_type = Chem.rdchem.BondType.values[bond]
        except AttributeError:
            ValueError(f'Unknown bond type {bond} in molecule.')

        mol.AddBond(node_to_idx[idx1], node_to_idx[idx2], bond_type)

    # Convert RWMol to Mol object
    mol = mol.GetMol()

    assert mol.GetNumAtoms() == len(graph), f'Number of atoms in rdkit molecule ({mol.GetNumAtoms()}) does not match number of atoms in molecule ({len(graph)})'
    return mol

def rdkit_mol_to_graph(mol: Chem.Mol, element_label: str='node_label', bond_label: str= 'bond_type') -> nx.Graph:
    """
    Create a graph from an rdkit mol object
    @param mol: RDKit mol object
    @param element_label: element label for node dictionary
    @param bond_label: bond type label for edge dictionary
    @return: networkx graph of the molecule
    """
    G = nx.Graph()

    for atom in mol.GetAtoms():
        node = atom.GetIdx()
        label = atom.GetSymbol()
        G.add_node(node, **{element_label: label})

    for bond in mol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()
        label = bond.GetBondTypeAsDouble()
        G.add_edge(u, v, **{bond_label: label})

    return G

def get_all_ligands_by_graph_hashes(all_ligands: list) -> dict:
    """
    Get dictionary of graph hashes with list of ligands.py with this graph hash
    :param all_ligands: list of all ligands.py with graph hashes
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

    # formula = [f"{el}{(c[el]) if c[el] != 1 else ''}" for el in elements] # drop the 1 if an element occurs only once
    formula = [f"{el}{(c[el])}" for el in elements]
    return "".join(formula)

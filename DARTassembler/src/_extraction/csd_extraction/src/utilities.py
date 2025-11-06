"""
All kind of utilities
"""
from typing import Union
from pathlib import Path

# from molSimplify.Classes.mol3D import mol3D
from rdkit.Chem.rdmolfiles import MolFromMol2File
import networkx as nx
from rdkit import Chem
from pymatgen.core.periodic_table import Element as Pymatgen_Element
from pysmiles import read_smiles

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def view_graph(G: nx.Graph, node_label='element', node_size=150):
    """
    Our standard graph visualizer
    """
    nx.draw_networkx(
        G,
        node_size=node_size,  # 500,
        with_labels=True,
        labels={node: G.nodes[node][node_label] for node in G.nodes}
    )
    plt.show()


def rdkit_mol_to_nx(mol: Chem.rdchem.Mol):
    """
    convert rdkit.chem Mol object to nx.Graph, as there is nothing built in
    But at least so we have full control over how the graphs should actually look lilke
    :param mol: The mol as an rdchem mol object we want to turn into a graph
    """
    G = nx.Graph()

    for atom in mol.GetAtoms():
        # iterating through atoms and adding a node to the graph for each atom
        # You can add each property you want, as a node attribute
        G.add_node(atom.GetIdx(),
                   node_label=Pymatgen_Element.from_Z(
                       atom.GetAtomicNum()).symbol,
                   atomic_num=atom.GetAtomicNum(),
                   # formal_charge=atom.GetFormalCharge(),  #is always set to 0
                   # chiral_tag=atom.GetChiralTag(),
                   # hybridization=atom.GetHybridization(),
                   # num_explicit_hs=atom.GetNumExplicitHs(),
                   # is_aromatic=atom.GetIsAromatic()
                   )

    for bond in mol.GetBonds():
        # similarly we proceed with bond and bond properties
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())

    return G

#
# def mol2_to_mol3D(mol2: Union[str, Path]) -> mol3D:
#     """
#     Convert mol2 file into molsimplify molecule
#     :param mol2: path to mol2 file
#     :return:
#     """
#
#     mol = mol3D()
#     return mol.readfrommol2(filename=mol2)


def mol2_to_rdkitmol(mol2: Union[str, Path]) -> Chem.rdchem.Mol:
    """
    Convert mol2 file into rdkit.Chem molecule molecule
    :param mol2: path to mol2 file
    :return:
    """

    return MolFromMol2File(mol2)


def smiles_to_graph(smiles: str) -> nx.Graph:
    """
    Convert smiles string into graph by employing readsmiles
    """

    G = read_smiles(smiles, explicit_hydrogen=True)
    # turns smiles into graph wit "element" as node_label

    return G

import warnings
import numpy as np
import json
import jsonlines
from pathlib import Path
from typing import Union
import networkx as nx
import matplotlib.pyplot as plt
import pysmiles
from copy import deepcopy



def load_json(path: Union[str, Path], n_max: int=None) -> dict:
    """
    Load a JSON or JSON Lines file. If the file is a JSON Lines file, it is converted to a dictionary.
    :param path: Path to the JSON or JSON Lines file
    :return: Dictionary with the contents of the file
    """
    # Accept False, None or np.inf to disable n_max
    if n_max is None or n_max is False:
        n_max = np.inf

    db = {}
    for i, (key, value) in enumerate(iterate_over_json(path)):
        if i >= n_max:
            break
        db[key] = value

    return db

def iterate_over_json(path: Union[str, Path]) -> tuple[str, dict]:
    """
    Iterate over a JSON or JSON Lines file and yield the key and value of each entry.
    :param path: Path to the JSON or JSON Lines file
    :return: Tuple with the key and value of each entry
    """
    try:
        # Try to load as normal JSON file first
        with open(path, 'r') as file:
            db = json.load(file)
            for key, value in db.items():
                yield key, value
    except json.JSONDecodeError:
        # If normal JSON fails, try to load as JSON Lines
        with jsonlines.open(path, 'r') as reader:
            for line in reader:
                # Since 'line' is a dictionary, no need to use json.loads
                yield line['key'], line['value']

    return

def view_graph(G, node_label='node_label', node_size=150):
    nx.draw_networkx(
        G,
        node_size=node_size,  # 500,
        with_labels=True,
        labels={node: G.nodes[node][node_label] for node in G.nodes}
    )
    plt.show()

def make_graph_labels_integers(G: [nx.Graph, nx.MultiGraph]):
    str_to_int_node_mapping = {node: int(node) for node in G.nodes}

    assert list(str_to_int_node_mapping.values()) == [int(key) for key in str_to_int_node_mapping.keys()]

    # now relabel them inplace (-> copy=False)
    nx.relabel_nodes(G, mapping=str_to_int_node_mapping, copy=False)

    return G

def graph_from_graph_dict(d):

    G_new = nx.from_dict_of_dicts(d["graph"])
    nx.set_node_attributes(G_new, d["node_attributes"])

    G_new = make_graph_labels_integers(G_new)

    #
    # Bring nodes in canonical order
    H = nx.Graph()
    H.add_nodes_from(sorted(G_new.nodes(data=True)))
    H.add_edges_from(G_new.edges(data=True))

    return H

bond_order_rdkit_to_pysmiles = {
    0: 0,
    1: 1,
    2: 2,
    3: 3,
    4: 4,
    12: 1.5
    }
def graph_to_smiles(graph, element_label='node_label', bond_label='bond_type'):
    """
    Convert networkx graph of a molecule to smiles string.
    @param graph:
    @param element_label:
    @param bond_label:
    @return: smiles string
    """
    # make a deepcopy of the graph to not change the original graph
    graph = deepcopy(graph)

    # convert bond orders to pysmiles format
    for _, _, edge in graph.edges(data=True):
        bo = edge[bond_label]
        if bo == 0 or not bo in bond_order_rdkit_to_pysmiles:
            warnings.warn(f'Bond order {bo} found in graph which is not specified in the code. Cannot convert to smiles.')
            return
        edge['order'] = bond_order_rdkit_to_pysmiles[bo]

    # workaround bug in pysmiles:
    # add zero edges between fragments, otherwise pysmiles will not work
    fragments_connectors = [list(frag)[0] for frag in nx.connected_components(graph)]
    central_node = fragments_connectors.pop(0)
    for idx in fragments_connectors:
        graph.add_edge(central_node, idx, order=0)

    # convert element labels to pysmiles format
    for _, atom in graph.nodes(data=True):
        atom['element'] = atom[element_label]

    # convert graph to smiles
    smiles = pysmiles.write_smiles(graph)

    # Bugfix of pysmiles: add brackets if smiles is not organic and a single atom
    n_atoms = len(graph)
    organic = 'H B C N O P S F Cl Br I'.split()
    if n_atoms == 1 and not smiles in organic:
        if not smiles.startswith('['):
            smiles = '[' + smiles + ']'  # add brackets if smiles is not organic


    # Doublecheck assumptions
    smiles_graph = pysmiles.read_smiles(smiles, explicit_hydrogen=True)
    # Check that number of fragments in graph matches number of fragments in smiles
    n_fragments = len(fragments_connectors) +1
    n_smiles_fragments = len(smiles.split('.'))
    assert n_fragments == n_smiles_fragments, f'Number of molecule fragments ({n_fragments}) does not match number of smiles fragments ({n_smiles_fragments}).'
    # Check that number of heavy atoms in graph matches number of heavy atoms in smiles
    heavy_atoms = sorted(atom[element_label] for _, atom in graph.nodes(data=True) if atom[element_label] != 'H')
    heavy_atoms_smiles = sorted(atom['element'] for _, atom in smiles_graph.nodes(data=True) if atom['element'] != 'H')
    assert heavy_atoms == heavy_atoms_smiles, f'Heavy atoms in graph ({heavy_atoms}) do not match heavy atoms in smiles ({heavy_atoms_smiles}).'
    # Check that the number of hydrogens in graph matches the number of hydrogens in smiles
    n_H = len(graph) - len(heavy_atoms)
    n_H_smiles = len(smiles_graph) - len(heavy_atoms_smiles)
    assert n_H == n_H_smiles, f'Number of hydrogens in graph ({n_H}) does not match number of hydrogens in smiles ({n_H_smiles}).'

    return smiles
import networkx as nx
import pandas as pd
from copy import deepcopy
from networkx import weisfeiler_lehman_graph_hash as graph_hash
from pymatgen.core.periodic_table import Element as Pymatgen_Element
from rdkit import Chem

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def view_graph(G, node_label='node_label', node_size=150):
    nx.draw_networkx(
        G,
        node_size=node_size,  # 500,
        with_labels=True,
        labels={node: G.nodes[node][node_label] for node in G.nodes}
    )
    plt.show()


def get_sorted_atoms_and_indices_from_graph(graph):
    nodes = pd.Series({i: el for i, el in graph.nodes.data('node_label')})
    nodes = nodes.sort_index()
    atoms = nodes.tolist()
    idc = nodes.index.tolist()

    return atoms, idc


def node_check(dict1, dict2):
    return dict1["node_label"] == dict2["node_label"]


def graphs_are_equal(G1, G2):
    """
    Exact comparison of graphs
    """
    return nx.is_isomorphic(G1, G2, node_match=node_check)


def graphs_are_equal_hash_version(G1, G2):
    """
    In theory lower accuracy than "graphs_are_equal", but way lower computational costs as well
    """
    return graph_hash(G1, node_attr='node_label', iterations=3, digest_size=16) == graph_hash(G2, node_attr='node_label', iterations=3, digest_size=16)


def find_node_in_graph_by_label(G: nx.Graph, label_to_find, expected_hits=None):

    dict_ = dict(G.nodes(data="node_label"))
    nodes = [key for key, value in dict_.items() if value == label_to_find]

    if expected_hits is None:
        return nodes
    else:
        assert len(nodes) == expected_hits, "Too many hits in graph search"

        if expected_hits == 1:
            return nodes.pop()

        return nodes


def graph_to_dict_with_node_labels(G, sort_dicts=True):
    """
    Problem: nx.to_dict_of_dicts doesnt preserve node labels
    """

    from src01.utilities import sorted_dict_of_dicts

    graph_dict = nx.to_dict_of_dicts(G)
    node_attributes = {node: G.nodes[node] for node in G.nodes}
    if sort_dicts:
        graph_dict = sorted_dict_of_dicts(graph_dict)
        node_attributes = sorted_dict_of_dicts(node_attributes)

    final_graph_dict = {"graph": graph_dict,
                        "node_attributes": node_attributes
                        }

    return final_graph_dict


def remove_node_features_from_graph(graph, keep: list=[], inplace=True):
    """
    Removes all node features from the given graph except for the ones specified in `keep`.
    :param graph: networkx multigraph with node features
    :param keep: list of node features which will not be removed
    :return:
    """
    if not inplace:
        graph = deepcopy(graph)

    node_attributes = graph.nodes(data=True)
    for _, attrs in node_attributes:
        props = list(attrs.keys())
        for prop in props:
            if not prop in keep:
                del attrs[prop]

    return graph


def remove_edge_features_from_graph(graph, keep=None, inplace=True):
    """
    Removes all edge features from the given graph except for the ones specified in `keep`.
    :param graph: networkx multigraph with node features
    :param keep: list of node features which will not be removed
    :return:
    """
    if keep is None:
        keep = []
    if not inplace:
        graph = deepcopy(graph)

    edge_attributes = graph.edges(data=True)
    for _, _, attrs in edge_attributes:
        if keep == []:
            attrs.clear()
        else:
            props = list(attrs.keys())
            for prop in props:
                if not prop in keep:
                    del attrs[prop]

    return graph


def make_multigraph_to_graph(graph) -> nx.Graph:
    if not isinstance(graph, nx.Graph):
        graph = nx.Graph(graph)
    return graph


def make_graph_labels_integers(G: [nx.Graph, nx.MultiGraph]):
    """
    This method makes the graph labels integers and also checks if they are labelled in the right way,
    i.e. from 0 to len(atom)-1
    or "0" to "len(atom)-1" if they are strings
    """
    # todo This function doesnt look like it really checks whether the graph labels are labeled from 0 to n-1. Makes sense for the relative indexing, but then remove it from the documentation?
    #
    str_to_int_node_mapping = {node: int(node) for node in G.nodes}
    # is required because sometimes (esp. from the readin process) the labels are denoted as str(int) and we need
    # to transform that
    #
    #
    # Das sichert tatsaechlich irgendwie auch die Grunddannahme, was ganz gut ist\
    assert list(str_to_int_node_mapping.values()) == [int(key) for key in str_to_int_node_mapping.keys()]

    # now relabel them inplace (-> copy=False)
    nx.relabel_nodes(G, mapping=str_to_int_node_mapping, copy=False)

    return G


def graph_from_graph_dict(d):
    # assert d hat keys graph, node_attrib sind eigentlich optional
    # assete dann aber dass list(G_new.nodes) == list(final_graph_dict["node_attributes"])

    G_new = nx.from_dict_of_dicts(d["graph"])
    nx.set_node_attributes(G_new, d["node_attributes"])

    G_new = make_graph_labels_integers(G_new)

    #
    # Bring nodes in canonical order
    H = nx.Graph()
    H.add_nodes_from(sorted(G_new.nodes(data=True)))
    H.add_edges_from(G_new.edges(data=True))

    return H


def get_reindexed_graph(graph):
    return nx.relabel.convert_node_labels_to_integers(graph, first_label=0, ordering='sorted')


def unify_graph(G):
    """
    THis method aims to bring graph from the tmQMG format in their .gml into the format we require for our
    process,

    i.e. it assures that the nodes are labelled by integers
    and that the resulting class is a nx.Graph object rather than a nx.MultiGraph object
    """

    # As in the "graph_from_graph_dict" method we need to make the graph nodes to integers
    G = make_graph_labels_integers(G)

    # and convert it to graph rather than Multigraph (as we dont care abount bond orders)
    G = nx.Graph(G)

    return G


def rdchem_mol_to_nx(mol: Chem.rdchem.Mol) -> nx.Graph:
    """
    convert rdkit.chem Mol object to nx.Graph, as there is nothing built in
    But at least so we have full control over how the graphs should actually look lilke
    :param mol: The mol as an rdchem mol object we want to turn into a graph
    """
    G = nx.Graph()

    for atom in mol.GetAtoms():
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
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())

    return G


def mol2_to_graph(filename):
    """
    Convert a mol2 file to a graph
    """
    from rdkit.Chem.rdmolfiles import MolFromMol2File

    mol2 = MolFromMol2File(filename, removeHs=False)
    return rdchem_mol_to_nx(mol2)


def mol2_str_to_graph(str_: str):
    """
    :param str_: The string of the mol2 file
    """
    from rdkit.Chem.rdmolfiles import MolFromMol2Block

    mol2 = MolFromMol2Block(str_, removeHs=False)
    return rdchem_mol_to_nx(mol2)

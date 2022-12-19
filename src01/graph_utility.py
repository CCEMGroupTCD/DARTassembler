import networkx as nx
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
from copy import deepcopy

def view_graph(G, node_label='node_label', node_size=150):
    nx.draw_networkx(
                        G,
                        node_size=node_size,#500,
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
    return nx.is_isomorphic(G1, G2, node_match=node_check)

def sorted_dict_of_dicts(d: dict) -> dict:
    """
    Sorts dictionaries and recursively sorts dictionaries of dictionaries to infinite order.
    :param d: dictionary to sort
    :return: sorted dictionary
    """
    sorted_d = {}
    keys = sorted(d.keys())

    for key in keys:
        value = d[key]

        if (isinstance(value, dict) and (len(value) > 1)):
           value = sorted_dict_of_dicts(value)

        sorted_d[key] = value

    assert (len(d) == len(sorted_d) and all([val == d[key] for key, val in sorted_d.items()])), 'Sorted dictionary is different than original one, there must be a bug.'
    return sorted_d


def graph_to_dict_with_node_labels(G, sort_dicts=True):
    """
    Problem: nx.to_dict_of_dicts doesnt preserve node labels
    """
    graph_dict = nx.to_dict_of_dicts(G)
    node_attributes = {node: G.nodes[node] for node in G.nodes}
    if sort_dicts:
        graph_dict = sorted_dict_of_dicts(graph_dict)
        node_attributes = sorted_dict_of_dicts(node_attributes)

    final_graph_dict = {"graph": graph_dict,
                        "node_attributes": node_attributes
                        }

    return final_graph_dict


def remove_node_features_from_graph(graph, keep: list=['node_label'], inplace=True):
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

def remove_edge_features_from_graph(graph, keep: list=[], inplace=True):
    """
    Removes all edge features from the given graph except for the ones specified in `keep`.
    :param graph: networkx multigraph with node features
    :param keep: list of node features which will not be removed
    :return:
    """
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


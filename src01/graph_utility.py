import networkx as nx
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def view_graph(G):
    nx.draw_networkx(G, node_size=500, with_labels=True,
                     labels={node: G.nodes[node]["node_label"] for node in G.nodes})
    plt.show()


def node_check(dict1, dict2):
    return dict1["node_label"] == dict2["node_label"]


def graphs_are_equal(G1, G2):
    return nx.is_isomorphic(G1, G2, node_match=node_check)


def graph_to_dict_with_node_labels(G):
    """
    Problem: nx.to_dict_of_dicts doesnt preserve node labels
    """
    final_graph_dict = {"graph": nx.to_dict_of_dicts(G),
                        "node_attributes": {node: G.nodes[node] for node in G.nodes}
                        }

    return final_graph_dict


def make_graph_labels_integers(G: [nx.Graph, nx.MultiGraph]):
    """
    This method makes the graph labels integers and also checks if they are labelled in the right way,
    i.e. from 0 to len(atom)-1
    or "0" to "len(atom)-1" if they are strings
    """
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

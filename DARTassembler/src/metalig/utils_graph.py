import networkx as nx
import pandas as pd
from copy import deepcopy
import pysmiles
import warnings
from typing import Union
from networkx import weisfeiler_lehman_graph_hash as graph_hash
from DARTassembler.src.constants.chem import Element

bond_order_rdkit_to_pysmiles = {
    0: 0,
    1: 1,
    2: 2,
    3: 3,
    4: 4,
    12: 1.5
    }

def get_only_complex_graph_connected_to_metal(graph, metal: str, atom_label='node_label') -> nx.Graph:
    """
    Returns the graph of only the metal complex without unconnected ligands.
    """
    if isinstance(metal, int) or isinstance(metal, float):
        raise ValueError(f"metal should be a string but is {metal}!")

    fragment_idc, fragment_els = get_graph_fragments(graph=graph, node_label=atom_label)
    complex_graph = None
    for frag_elements, frac_indices in zip(fragment_els, fragment_idc):
        if metal in frag_elements:
            assert complex_graph is None, "There should only be one fragment with the metal!"
            complex_graph = deepcopy(graph.subgraph(frac_indices))

    if complex_graph is None:
        raise ValueError(f"There is no fragment with the metal {metal}!")

    return complex_graph

def get_heavy_atoms_graph(graph: nx.Graph, element_label='node_label') -> nx.Graph:
    """
    Returns a graph with all heavy atoms (non-hydrogen) of the given graph.
    """
    graph = nx.Graph(graph)     # copy and unfreeze graph
    graph.remove_nodes_from(list(idx for idx, atom in graph.nodes(data=True) if atom[element_label] == 'H'))

    return graph

def assert_graph_and_coordinates_are_consistent(graph: nx.Graph, atoms: list[str], graph_hash: str, donor_idc: list[int]=None, node_label='node_label'):
    """
    Check if the graph and the coordinates are consistent, i.e. if the graph has the same atoms in the same order as the atomic_props and if the graph has the correct coordinating atom indices
    :param graph: networkx graph
    :param atoms: list of atoms in the same order as the atomic_props
    :param graph_hash: hash of the graph, to check if the graph has not changed
    :param donor_idc: list of indices of the coordinating atoms in the ligand, if the graph is a ligand graph. If None, no check is performed.
    :param node_label: the label of the node attribute that contains the atom labels, e.g. 'node_label'.
    :return: None if the graph and coordinates are consistent, raises AssertionError otherwise.
    """
    # Check if graph hashes are the same
    new_graph_hash = get_graph_hash(graph)
    assert new_graph_hash == graph_hash, f'Graph hashes don\'t match: {new_graph_hash} vs {graph_hash}'

    # Check if the relabeled graph has the same atoms in the same order as the atomic_props
    relabeled_atoms = [a for _, a in sorted(graph.nodes(data=node_label))]
    assert relabeled_atoms == atoms, f'Atoms in reindexed graph and atomic_props don\'t match: {relabeled_atoms} vs {atoms}'

    # Check if the relabeled graph has the correct coordinating atom indices, for ligands only
    if not donor_idc is None:
        relabeled_coord_atoms = [i for i, a in graph.nodes(data='metal_neighbor') if a == True]
        if len(relabeled_coord_atoms) > 0:
            assert relabeled_coord_atoms == donor_idc, f'Coordinating atoms in relabeled graph and ligand don\'t match: {relabeled_coord_atoms} vs {donor_idc}'

    return

def get_graph_hash(graph, node_attr='node_label', iterations=3, digest_size=16, edge_attr=None) -> str:
    """
    Returns a Weisfeiler-Lehman graph hash of the given graph, e.g 'd3f4e5b6c7d8e9f0'.
    :param graph: graph
    :param node_attr: node attribute to be used for the hash
    :param iterations: iterations of the hash
    :param digest_size: digest size of the hash
    :return: Weisfeiler-Lehman graph hash of the given graph.
    """
    return graph_hash(graph, node_attr=node_attr, edge_attr=edge_attr, iterations=iterations, digest_size=digest_size)

def smiles2nx(smiles_str: str, explicit_H: bool = True) -> nx.Graph:
    """
    Convert a SMILES string to a networkx graph.
    """
    G = nx.Graph(pysmiles.read_smiles(smiles_str, explicit_hydrogen=explicit_H))
    for node in G.nodes:
        G.nodes[node]["node_label"] = G.nodes[node]["element"]

    return G

def view_graph(G: nx.Graph, node_label='node_label', node_size: int=150, save_path: str=None) -> None:
    """
    Visualizes a networkx graph using matplotlib. The nodes are labeled with the node_label attribute.
    :param G: networkx graph to visualize
    :param node_label: the label of the node attribute that contains the atom labels, e.g. 'node_label'.
    :param node_size: size of the nodes in the visualization.
    :param save_path: path to save the visualization. If None, the visualization is not saved but displayed.
    :return: None, displays the graph using matplotlib.
    """
    import matplotlib
    try:  # Avoid error when running on server
        matplotlib.use('TkAgg')
    except ImportError:
        pass
    import matplotlib.pyplot as plt
    nx.draw_networkx(
        G,
        node_size=node_size,  # 500,
        with_labels=True,
        labels={node: G.nodes[node][node_label] for node in G.nodes}
    )
    if not save_path is None:
        plt.savefig(save_path)
    plt.show()

    return

def get_graph_fragments(graph, node_label='node_label') -> tuple[list[list[int]], list[list[str]]]:
    """
    Returns a list of the fragments (unconnected components) of the molecular graph. Returns both the indices and the elements, if atom_label is provided.
    :param graph: networkx graph
    :param node_label: the label of the node attribute that contains the atom labels, e.g. 'node_label' or 'element'.
    """
    fragment_idc = [sorted(comp) for comp in nx.connected_components(graph)]
    fragment_els = [[graph.nodes[i][node_label] for i in fragment] for fragment in fragment_idc]
    return fragment_idc, fragment_els

def get_sorted_atoms_and_indices_from_graph(graph, atom_label='node_label'):
    nodes = pd.Series({i: el for i, el in graph.nodes.data(atom_label)})
    nodes = nodes.sort_index()
    atoms = nodes.tolist()
    idc = nodes.index.tolist()

    return atoms, idc

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

def graph_to_dict_with_node_labels(G: nx.Graph, sort_dicts: bool=True) -> dict:
    """
    Converts a networkx graph to a dictionary representation, preserving node labels and attributes.
    :param G: networkx graph to convert.
    :param sort_dicts: whether to sort the dictionaries of the graph and node attributes.
    :return: dictionary representation of the graph with node labels and attributes.
    """
    from DARTassembler.src.metalig.utils import sorted_dict_of_dicts
    graph_dict = nx.to_dict_of_dicts(G)
    node_attributes = {node: G.nodes[node] for node in G.nodes}
    if sort_dicts:
        graph_dict = sorted_dict_of_dicts(graph_dict)
        node_attributes = sorted_dict_of_dicts(node_attributes)

    final_graph_dict = {"graph": graph_dict,
                        "node_attributes": node_attributes
                        }

    return final_graph_dict


def remove_node_features_from_graph(graph, keep: list=None, inplace: bool=True):
    """
    Removes all node features from the given graph except for the ones specified in `keep`.
    :param graph: networkx multigraph with node features.
    :param keep: list of node features which will not be removed.
    :return: networkx multigraph with node features.
    """
    if keep is None:
        keep = []
    if not inplace:
        graph = deepcopy(graph)

    node_attributes = graph.nodes(data=True)
    for _, attrs in node_attributes:
        props = list(attrs.keys())
        for prop in props:
            if not prop in keep:
                del attrs[prop]

    return graph


def remove_edge_features_from_graph(graph, keep: list=None, inplace: bool=True):
    """
    Removes all edge features from the given graph except for the ones specified in `keep`.
    :param graph: networkx multigraph with node features.
    :param keep: list of node features which will not be removed.
    :return: networkx multigraph with node features.
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
    str_to_int_node_mapping = {node: int(node) for node in G.nodes}
    assert list(str_to_int_node_mapping.values()) == [int(key) for key in str_to_int_node_mapping.keys()]
    nx.relabel_nodes(G, mapping=str_to_int_node_mapping, copy=False)

    return G

def get_adjacency_matrix(graph):
    assert list(graph.nodes) == sorted(
        graph.nodes), 'Nodes are not sorted, this might lead to problems with the indexing of the adjacency matrix.'

    with warnings.catch_warnings():  # ignore FutureWarning from networkx
        warnings.simplefilter(action='ignore', category=FutureWarning)
        A = nx.to_numpy_array(graph)

    return A

def graph_from_graph_dict(d):

    # The input dictionary has the nodes as strings, convert them to integers because everything else is unintuitive
    d['graph'] = {int(str_node): {int(str_neighbor): d['graph'][str_node][str_neighbor] for str_neighbor in d['graph'][str_node]} for str_node in d['graph']}
    d['node_attributes'] = {int(str_node): d['node_attributes'][str_node] for str_node in d['node_attributes']}

    # Create graph from dictionary
    G = nx.from_dict_of_dicts(d["graph"])
    nx.set_node_attributes(G, d["node_attributes"])

    # Validate graph
    assert sorted(list(G.nodes())) == list(G.nodes()), "Nodes are not sorted"
    assert all(isinstance(node, int) for node in G.nodes), "Nodes are not integers"

    return G

def get_reindexed_graph(graph):
    """
    Reindex the given graph so that the order of the indices stays the same but the indices now go from 0 to n-1. This is a very important function because the graphs of the ligands currently keep the indices from when they were in the complex, which means their indices do not go from 0 to n-1. This means the indices of the atoms in the graphs and the atoms in the atomic_props dict is different, just the order is the same.
    :param graph: graph to reindex.
    :return: reindexed graph with indices from  0 to n-1 in the same numerical order as before.
    """
    assert all(isinstance(node, int) for node in graph.nodes), f'Graph nodes are not integers: {list(graph.nodes)}.'

    # Skip reindexing if the nodes are already indexed from 0 to n-1
    if sorted(list(graph.nodes)) == list(range(len(graph.nodes))):
        return graph

    old_labels = sorted(list(graph.nodes))
    mapping = {old_label: new_label for new_label, old_label in enumerate(old_labels)}
    reindexed_graph = nx.relabel_nodes(graph, mapping)

    # Create a new graph with sorted nodes
    # Order is important: First add all nodes without edges
    sorted_graph = nx.Graph()
    for node in sorted(reindexed_graph.nodes):
        sorted_graph.add_node(node, **reindexed_graph.nodes[node])
    # Now add edges
    for node in sorted(reindexed_graph.nodes):
        for neighbor, edge_attrs in reindexed_graph[node].items():
            sorted_graph.add_edge(node, neighbor, **edge_attrs)

    nodes = list(sorted_graph.nodes)
    assert nodes == sorted(nodes), f'Nodes are not sorted after reindexing: {nodes}'
    assert list(sorted_graph.nodes) == list(range(len(nodes))), f'Nodes are not indexed from 0 to n-1 after reindexing: {nodes}'
    return sorted_graph

def count_atoms_with_n_bonds(graph: nx.Graph, element: Union[str, None], n_bonds: int, graph_element_label: str='node_label') -> int:
    """
    Count the number of n_ligand_instances of element `element` with exactly `n_bonds` bonds in the given molecular graph.
    :param graph (network.Graph): molecular graph.
    :param element (str, None): specification of the element, e.g. 'C'. If None, all elements are counted.
    :param n_bonds (int): count an atom if it has exactly this number of bonds.
    :param graph_element_label (str): the label of the element string in the graph attributes. Only necessary if `element` is not `None`.
    :return (int): integer count of the n_ligand_instances.
    """
    n = 0
    for atom_idx, atom in graph.nodes(data=True):
        if element is None or atom[graph_element_label] == element:
            n_atom_bonds = len(list(nx.all_neighbors(graph, atom_idx)))
            if n_atom_bonds == n_bonds:
                n += 1

    return n

def count_atoms_with_n_bonds(graph: nx.Graph, element: Union[str, None], n_bonds: int, graph_element_label: str='node_label') -> int:
    """
    Count the number of n_ligand_instances of element `element` with exactly `n_bonds` bonds in the given molecular graph.
    :param graph (network.Graph): molecular graph.
    :param element (str, None): specification of the element, e.g. 'C'. If None, all elements are counted.
    :param n_bonds (int): count an atom if it has exactly this number of bonds.
    :param graph_element_label (str): the label of the element string in the graph attributes. Only necessary if `element` is not `None`.
    :return (int): integer count of the n_ligand_instances.
    """
    n = 0
    for atom_idx, atom in graph.nodes(data=True):
        if element is None or atom[graph_element_label] == element:
            n_atom_bonds = len(list(nx.all_neighbors(graph, atom_idx)))
            if n_atom_bonds == n_bonds:
                n += 1

    return n

def unify_graph(G):
    """
    THis method aims to bring graph from the tmQMG format in their .gml into the format we require for our process, i.e. it assures that the nodes are labeled by integers and that the resulting class is a nx.Graph object rather than a nx.MultiGraph object.
    """
    G = make_graph_labels_integers(G)
    G = nx.Graph(G)

    return G

def rdchem_mol_to_nx(mol) -> nx.Graph:
    """
    convert rdkit.chem Mol object to nx.Graph, as there is nothing built in
    But at least so we have full control over how the graphs should actually look lilke
    :param mol: The atoms as an rdchem atoms object we want to turn into a graph
    """
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   node_label=Element(
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

def graph_to_smiles(graph, element_label='node_label', bond_label='bond_type') -> Union[str,None]:
    """
    Convert networkx graph of a molecule to a SMILES string. This function is currently not deterministic, i.e. it can return different smiles strings for the same molecule.
    :param graph: networkx graph
    :param element_label: element label of the graph
    :param bond_label: bond label of the graph
    :return: SMILES string or None if any bond orders are invalid or unkwown.
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

    # convert element labels to pysmiles format
    for _, atom in graph.nodes(data=True):
        atom['element'] = atom[element_label]

    # convert graph to smiles
    all_smiles = []
    for nodes in nx.connected_components(graph):
        # workaround bug in pysmiles:
        # works only for individual fragments, therefore convert each fragment separately
        single_graph = graph.subgraph(nodes)
        smiles = pysmiles.write_smiles(single_graph)
        all_smiles.append(smiles)
    smiles = '.'.join(sorted(all_smiles))

    # Bugfix of pysmiles: add brackets if smiles is not organic and a single atom
    n_atoms = len(graph)
    organic = 'H B C N O P S F Cl Br I'.split()
    if n_atoms == 1 and not smiles in organic:
        if not smiles.startswith('['):
            smiles = '[' + smiles + ']'  # add brackets if smiles is not organic


    # Doublecheck assumptions
    smiles_graph = pysmiles.read_smiles(smiles, explicit_hydrogen=True)
    # Check that number of heavy atoms in graph matches number of heavy atoms in smiles
    heavy_atoms = sorted(atom[element_label] for _, atom in graph.nodes(data=True) if atom[element_label] != 'H')
    heavy_atoms_smiles = sorted(atom['element'] for _, atom in smiles_graph.nodes(data=True) if atom['element'] != 'H')
    assert heavy_atoms == heavy_atoms_smiles, f'Heavy atoms in graph ({heavy_atoms}) do not match heavy atoms in smiles ({heavy_atoms_smiles}).'
    # Check that the number of hydrogens in graph matches the number of hydrogens in smiles
    n_H = len(graph) - len(heavy_atoms)
    n_H_smiles = len(smiles_graph) - len(heavy_atoms_smiles)
    assert n_H == n_H_smiles, f'Number of hydrogens in graph ({n_H}) does not match number of hydrogens in smiles ({n_H_smiles}).'

    return smiles

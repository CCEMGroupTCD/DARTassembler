from copy import deepcopy
import pysmiles
import warnings
import networkx as nx

bond_order_rdkit_to_pysmiles = {
    0: 0,
    1: 1,
    2: 2,
    3: 3,
    4: 4,
    12: 1.5
    }



def graph_to_smiles(graph, element_label='node_label', bond_label='bond_type') -> str:
    """
    Convert networkx graph of a molecule to smiles string. This function is currently not deterministic, i.e. it can return different smiles strings for the same molecule.
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


# if __name__ == '__main__':
#
#     db_version = '1.6'
#     db_path = f'../data/final_db_versions/complex_db_v{db_version}.json'
#
#
#     # db_path = '/Users/timosommer/PhD/projects/RCA/projects/CreateTMC/data_output/CSD_MM_G_Jsons_test/complex_db.json'
#     df = pd.DataFrame.from_dict(load_complex_db(path=db_path), orient='index')
#     data = load_complex_db(path=db_path, molecule='class')
#     # df = df[df['denticity'] > 0]
#     # lig = data['unq_CSD-FIXVAL-02-a']
#
#     #%%
#     ligands = {}
#     all_bond_orders = []
#     for name, lig in data.items():
#         bond_orders = [edge['bond_type'] for _, _, edge in lig.graph.edges(data=True)]
#         all_bond_orders.extend(bond_orders)
#         ligands[name] = {
#                             'bond_orders': bond_orders,
#                             'good_bond_orders': lig.check_for_good_bond_orders(),
#         }
#     df_ligs = pd.DataFrame.from_dict(ligands, orient='index')
#     df = df.join(df_ligs)
#
#     #%%
#     frac_good_bond_orders = df['good_bond_orders'].sum() / len(df)
#     good_bond_orders = pd.Series(all_bond_orders).value_counts()
#     print(f'Fraction of molecules with good bond orders: {frac_good_bond_orders:.2f}')
#     print(f'Good bond orders: {good_bond_orders}')
#
#     #%% plot bond orders
#     plt.figure()
#     sns.histplot(all_bond_orders)
#     plt.savefig('/Users/timosommer/Downloads/bond_orders.png', dpi=300)
#
#
#
#
#     lig = df.iloc[0]
#
#     bond_orders = {name: [edge['bond_type'] for _, edge in mol.graph.edges(data=True)] for name, mol in data.items()}
#     ligands = {}
#     for name, mol in data.items():
#         _, _, edge_dict = mol.graph.nodes(data=True)
#         bond_orders = edge_dict['bond_type']
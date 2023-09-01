import pandas as pd
from pathlib import Path
from sklearn.metrics import accuracy_score
from pysmiles import read_smiles
from networkx import weisfeiler_lehman_graph_hash as graph_hash
from DARTassembler.src.ligand_extraction.utilities_graph import smiles2nx
from pymatgen.core.composition import Composition


def remove_entries_with_multiple_charges_per_graph_hash_and_remove_duplicate_graph_hashes(df, charge_col):
    charges_per_graph_hash = df.groupby('graph_hash')[charge_col].agg(['unique']).reset_index().rename(columns={'unique': 'all_charges'})
    df = df.merge(charges_per_graph_hash, left_on='graph_hash', right_on='graph_hash', how='left')
    n_charges = df['all_charges'].apply(len)
    df = df[n_charges < 2]
    df = df.drop_duplicates('graph_hash')

    return df

if __name__ == '__main__':

    ms_lig_path = '/Users/timosommer/PhD/projects/RCA/projects/CreateTMC/test/tests_for_paper/data/molSimplify/OctLig.csv'
    ulig_db_version = '1.6'
    only_confident_charges = True
    exclude_unconnected_ligands = True


    df_uligs = pd.read_json(f'../../data/final_db_versions/unique_ligand_db_v{ulig_db_version}.json', orient='index')
    if only_confident_charges:
        df_uligs = df_uligs.query('pred_charge_is_confident')
    if exclude_unconnected_ligands:
        df_uligs = df_uligs.query('denticity > 0')
    df_uligs = remove_entries_with_multiple_charges_per_graph_hash_and_remove_duplicate_graph_hashes(df_uligs, charge_col='pred_charge')

    df_ms = pd.read_csv(ms_lig_path)
    df_ms['graph'] = df_ms['SMILES'].apply(smiles2nx)
    df_ms['graph_hash'] = df_ms['graph'].apply(lambda graph: graph_hash(graph, node_attr='element', iterations=3, digest_size=16))
    df_ms = remove_entries_with_multiple_charges_per_graph_hash_and_remove_duplicate_graph_hashes(df_ms, charge_col='charge')
    assert len(df_ms.drop_duplicates('graph_hash')) == len(df_ms.drop_duplicates(['graph_hash', 'charge'])), 'Graph hash -> charge mapping is not unique.'

    df_same_graph = df_uligs.merge(df_ms, left_on='graph_hash', right_on='graph_hash', how='inner')
    assert len(df_same_graph.drop_duplicates('graph_hash')) == len(df_same_graph.drop_duplicates(['graph_hash', 'charge'])), 'Graph hash -> charge mapping is not unique.'
    df_same_graph
    assert (df_same_graph['stoichiometry'].apply(Composition) == df_same_graph['formula'].apply(Composition)).all(), 'Stoichiometries are different for the same graph hash!'
    df_same_graph['same_catoms'] = df_same_graph['local_elements'].apply(sorted) == df_same_graph['catoms'].apply(lambda s: sorted(s.split(',')))
    df_same_graph = df_same_graph[df_same_graph['same_catoms']]
    assert df_same_graph['diff_catoms'].all()
    assert all(df_same_graph['denticity_x'] == df_same_graph['denticity_y'])
    df_same_graph['diff_charge'] = df_same_graph['pred_charge'] == df_same_graph['charge']

    df_same_graph = df_same_graph[['stoichiometry',  'graph_hash', 'pred_charge', 'charge',  'radicals', 'occurrences', 'diff_charge','all_ligand_names',
        'n_protons', 'denticity_x', 'name_x',
       'ligand_to_metal', 'local_elements', 'was_connected_to_metal',
       'original_metal_position', 'graph_hash_with_metal', 'stats',
       'unique_name',  'same_graph_denticities', 'count_metals',
       'n_same_graph_denticities', 'n_metals', 'n_same_graphs',
       'has_unconnected_ligands',
       'pred_charge_is_confident', 'all_charges_x', 'graph_dict', 'atomic_props', 'global_props',  'name_y', 'SMILES',
       '_id', 'n_heavy_atoms', 'T1', 'TAE', 'max(t1)',
       'C0^2', 'nHOMO_CAS', 'nLUMO_CAS', 'IND_PBE', 'rND_PBE', 'IND_b3lyp',
       'rND_b3lyp', 'B1', 'A25PBE', 'nHOMO_MP2', 'nLUMO_MP2', '%Ecorr[(T)]',
       'lacRACs.f-chi-0-all', 'lacRACs.f-chi-1-all',
       'lacRACs.f-chi-2-all', 'lacRACs.f-chi-3-all', 'lacRACs.f-Z-0-all',
       'lacRACs.f-Z-1-all', 'lacRACs.f-Z-2-all', 'lacRACs.f-Z-3-all',
       'lacRACs.f-I-0-all', 'lacRACs.f-I-1-all', 'lacRACs.f-I-2-all',
       'lacRACs.f-I-3-all', 'lacRACs.f-T-0-all', 'lacRACs.f-T-1-all',
       'lacRACs.f-T-2-all', 'lacRACs.f-T-3-all', 'lacRACs.f-S-0-all',
       'lacRACs.f-S-1-all', 'lacRACs.f-S-2-all', 'lacRACs.f-S-3-all', 'n_bond',
       'nus_bond', 'avrg_bo', 'catoms', 'denticity_y', 'train', 'graph',
       'all_charges_y', 'diff_catoms', 'same_catoms']]

    ms_lig = df_ms.iloc[0,:].to_dict()

    acc = accuracy_score(y_true=df_same_graph['pred_charge'], y_pred=df_same_graph['charge'])







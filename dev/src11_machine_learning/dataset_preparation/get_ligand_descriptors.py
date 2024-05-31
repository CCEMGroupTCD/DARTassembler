"""
This script is for generating a dataset for machine learning of the input complexes.
"""
import numpy as np
import pandas as pd
from copy import deepcopy
import networkx as nx
from DARTassembler.src.ligand_extraction.bond_orders import graph_to_smiles
from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db, load_complex_db, load_full_ligand_db, save_unique_ligand_db, save_full_ligand_db, save_complex_db
from DARTassembler.src.ligand_extraction.utilities import unroll_dict_into_columns
from dev.src11_machine_learning.dataset_preparation.descriptors import SOAP_3D, RDKit_2D
import random

if __name__ == '__main__':

    db_version = '1.7'
    n_molecules = None      # None for all complexes
    dataset = f'../../data/final_db_versions/unique_ligand_db_v{db_version}.json'
    out_csv = f'../../test/ML_detect_missing_H/data/ligand_descriptors_v{db_version}.csv'

    print('Load ligand database.')
    db = load_unique_ligand_db(dataset, 'class', n_max=n_molecules)

    # print('Start creating SOAP descriptors for each complex.')
    # ase_mols = [lig.get_ase_molecule(remove_elements=['H']) for lig in db.values()]
    # r_cut = 10.0
    # n_max= 2
    # l_max= 1
    # crossover= False
    # average= 'inner'
    # weighting = None
    # soap = SOAP_3D(
    #                 ase_molecules=ase_mols,
    #                 r_cut=r_cut,
    #                 n_max=n_max,
    #                 l_max=l_max,
    #                 crossover=crossover,
    #                 average=average,
    #                 weighting=weighting
    #                 )
    # soap_desc = soap.calculate_descriptors()
    # print('Created SOAP descriptors successfully.')

    # smiles_molecules = [graph_to_smiles(lig.graph, element_label='node_label', bond_label='bond_type') for lig in db.values()]

    db_dict = []
    for name, mol in db.items():
        db_dict.append({
                        'name': name,
                        'denticity': mol.denticity,
                        'pred_charge': mol.pred_charge,
                        'pred_charge_is_confident': mol.pred_charge_is_confident,
                        'n_hydrogens': mol.n_hydrogens,
                        'n_missing_H': 0,
                        'graph': deepcopy(mol.graph)
                        })
    df_real = pd.DataFrame(db_dict)

    # Filter out unwanted molecules
    # Filter ligands with not confident charge prediction
    df_real = df_real[df_real['pred_charge_is_confident']]
    # Filter ligands consisting only of H
    df_real = df_real[df_real['graph'].apply(lambda graph: pd.unique([el for _, el in graph.nodes(data='node_label')]).tolist() != ['H'])]


    # Add molecules with missing H
    n_removed_H = 2
    db_dict = []
    for _, mol in df_real.iterrows():
        graph = deepcopy(mol.graph)
        H_nodes = [node for node in graph.nodes if graph.nodes[node]['node_label'] == 'H']
        n_actual_removed_H = n_removed_H if len(H_nodes) >= n_removed_H else len(H_nodes)

        remove_H_nodes = random.sample(H_nodes, n_actual_removed_H)
        graph.remove_nodes_from(remove_H_nodes)

        # Double check that the number of H is correct
        original_atoms_list = [el for _, el in mol.graph.nodes(data='node_label')]
        new_atoms_list = [el for _, el in graph.nodes(data='node_label')]
        assert sorted(original_atoms_list) == sorted(new_atoms_list+['H']*n_actual_removed_H)

        db_dict.append(mol.to_dict())
        db_dict.append({
                        'name': name,
                        'denticity': mol.denticity,
                        'pred_charge': mol.pred_charge,
                        'pred_charge_is_confident': mol.pred_charge_is_confident,
                        'n_hydrogens': mol.n_hydrogens - n_actual_removed_H,
                        'n_missing_H': n_actual_removed_H,
                        'graph': graph
                        })
    df = pd.DataFrame(db_dict)

    # Add smiles from graph
    df['smiles'] = df['graph'].apply(lambda x: graph_to_smiles(x, element_label='node_label', bond_label='bond_type'))
    df = df[df['smiles'].notna()]
    df = df.drop(columns=['graph'])

    # Add RDKit 2D descriptors
    rdkit_desc = RDKit_2D(df['smiles'].values).compute_2Drdkit()
    rdkit_descriptors = rdkit_desc.columns.tolist()
    df = df.merge(rdkit_desc, left_on='smiles', right_index=True, how='left')
    df = df[df[rdkit_descriptors].notna().all(axis=1)]

    # soap_cols = soap.get_descriptor_names()
    # soap_desc = pd.DataFrame(data=soap_desc, columns=soap_cols, index=df.index)
    # df = df.join(soap_desc, validate='1:1')

    print(f'Save csv with molecules, properties and descriptors to {out_csv}.')
    df.to_csv(out_csv, index=False)

    print('Done!')










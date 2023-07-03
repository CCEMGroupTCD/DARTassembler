"""
This script is for generating a dataset for machine learning of the input complexes.
"""
import numpy as np
import pandas as pd
from src01.io_custom import load_unique_ligand_db, load_complex_db, load_full_ligand_db, save_unique_ligand_db, save_full_ligand_db, save_complex_db
from src01.utilities import unroll_dict_into_columns
from src11_machine_learning.dataset_preparation.descriptors import SOAP_3D

if __name__ == '__main__':

    db_version = '1.7'
    n_complexes = 100      # None for all complexes
    tmqm_csv = f'../../../data_input/tmQM/global_mol_properties.csv'  # all complexes in the tmqm as csv
    complex_dataset = f'../../../data/final_db_versions/complex_db_v{db_version}.json' # all complexes
    ligand_dataset = f'../../../data/final_db_versions/unique_ligand_db_v{db_version}.json'  # all complexes



    df = pd.read_csv(tmqm_csv, index_col='CSD_code')

    print('Load complex database.')
    complexes = load_complex_db(complex_dataset, 'class', n_max=n_complexes)
    ligands = load_unique_ligand_db(ligand_dataset, 'class')



    print(f'Save csv with complexes, properties and descriptors to {out_csv}.')
    df.to_csv(out_csv, index=False)

    print('Done!')










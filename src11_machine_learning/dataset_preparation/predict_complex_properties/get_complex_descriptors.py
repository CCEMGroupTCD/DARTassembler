"""
This script is for generating a dataset for machine learning of the input complexes.
"""
import numpy as np
import pandas as pd
from src01.io_custom import load_unique_ligand_db, load_complex_db


if __name__ == '__main__':

    db_version = '1.7'
    n_complexes = None      # None for all complexes
    test_frac = 0.1
    tmqm_csv = f'../../../data_input/tmQM/global_mol_properties.csv'  # all complexes in the tmqm as csv
    complex_dataset = f'../../../data/final_db_versions/complex_db_v{db_version}.json' # all complexes
    ligand_dataset = f'../../../data/final_db_versions/unique_ligand_db_v{db_version}.json'  # all unique ligands
    out_csv = f'../../../test/machine_learning/complex_property_prediction/input_data/complexes_with_descriptors_v{db_version}.csv'




    df_tmqm = pd.read_csv(tmqm_csv, index_col='CSD_code').to_dict(orient='index')

    print('Load complex database.')
    complexes = load_complex_db(complex_dataset, 'class', n_max=n_complexes)
    unique_ligand_names = [[lig.unique_name for lig in c.ligands] for c in complexes.values()]
    unique_ligand_names = pd.unique(np.concatenate(unique_ligand_names)).tolist()
    ligands = load_unique_ligand_db(ligand_dataset, 'class', n_max=unique_ligand_names)

    print('Generate descriptors for complexes.')
    all_complexes_with_descriptors = {}
    skipped_complexes = []
    for c_id, c in complexes.items():
        try:
            properties = df_tmqm[c_id]
        except KeyError:
            continue
        if not c.has_good_bond_orders:
            skipped_complexes.append(c_id)
            continue
        descriptors = c.generate_descriptors_of_complex_graph(only_core_complex=True)
        all_complexes_with_descriptors[c_id] = {**properties, **descriptors}
    df = pd.DataFrame.from_dict(all_complexes_with_descriptors, orient='index')

    # Add test set specification
    df = df.sample(frac=1, random_state=0)
    n_test = int(len(df) * test_frac)
    df['is_test_set'] = [True] * n_test + [False] * (len(df) - n_test)
    df = df.sample(frac=1, random_state=0)

    print(f'Save csv with complexes, properties and descriptors to {out_csv}.')
    df.to_csv(out_csv)

    print('Done!')










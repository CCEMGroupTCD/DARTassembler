"""
This script is for generating a dataset for machine learning of the input complexes.
"""
import pandas as pd
import numpy as np
from src01.DataBase import LigandDB
from src01.io_custom import load_unique_ligand_db, load_complex_db, iterate_complex_db
from tqdm import tqdm

if __name__ == '__main__':

    db_version = '1.7'
    n_complexes = 100      # None for all complexes
    test_frac = 0.1
    tmqm_csv = f'../../../data_input/tmQM/global_mol_properties.csv'  # all complexes in the tmqm as csv
    complex_dataset = f'../../../data/final_db_versions/complex_db_v{db_version}.json' # all complexes
    ligand_dataset = f'../../../data/final_db_versions/unique_ligand_db_v{db_version}.json'  # all unique ligands
    out_csv = f'../../../test/machine_learning/complex_property_prediction/input_data/test_complexes_with_descriptors_v{db_version}.csv'
    only_core_complex = True


    ################# Here the script starts #################

    print('Load tmqm and ligand database.')
    df_tmqm = pd.read_csv(tmqm_csv, index_col='CSD_code').to_dict(orient='index')
    ligands = LigandDB.load_from_json(path=ligand_dataset, n_max=50, only_core_ligands=only_core_complex)
    xtbs = []
    for ulig in tqdm(ligands.db.values(), desc='Generate descriptors for ligands'):
        xtbs.append(ulig.get_xtb_descriptors())
    xtbs = pd.DataFrame(xtbs, index=ligands.db.keys())
    n_nan = xtbs.isna().any(axis=1).sum()
    print(f'Number of NaN values in XTB descriptors: {n_nan} out of {len(xtbs)}')
    # soap_features = ligands.get_SOAP_descriptors(only_from_donors=True).to_dict(orient='index')
    #
    # print('Generate descriptors for complexes.')
    # all_complexes_with_descriptors = {}
    # np.random.seed(0)     # Set random seed for reproducibility of test and train set
    # for c_id, c in iterate_complex_db(complex_dataset, 'class', n_max=n_complexes):
    #     try:
    #         properties = df_tmqm[c_id]
    #     except KeyError:
    #         continue
    #
    #     # Generate test and train set for complexes
    #     properties['is_test_set'] = np.random.rand() < test_frac
    #
    #     # Generate descriptors for complex
    #     descriptors = c.generate_descriptors_of_complex_graph(only_core_complex=only_core_complex)
    #
    #     # Add SOAP descriptors of unique ligands making up the complex
    #     soaps = []
    #     for lig in c.ligands:
    #         uname = lig.unique_name
    #         if lig.denticity > 0 or not only_core_complex:
    #             soap_labels = list(soap_features[uname].keys())
    #             soap_values = list(soap_features[uname].values())
    #             soaps.append(soap_values)
    #     soaps = np.mean(soaps, axis=0)
    #     descriptors.update({label: value for label, value in zip(soap_labels, soaps)})
    #
    #     all_complexes_with_descriptors[c_id] = {**properties, **descriptors}
    # df = pd.DataFrame.from_dict(all_complexes_with_descriptors, orient='index')
    #
    # # Remove descriptors with missing values
    # drop_col = [col for col in descriptors.keys() if df[col].isna().any()]
    # df = df.drop(columns=drop_col)
    # print(f'Removed {len(drop_col)} descriptors with missing values from df: {drop_col}')
    #
    # # Rename columns since we only use the sum of the descriptors on atomistic level
    # df = df.rename(columns={col:col.rstrip('-sum') for col in df.columns if col.startswith('own_graph')})
    #
    # print(f'Save csv with complexes, properties and descriptors to {out_csv}.')
    # df.to_csv(out_csv)
    #
    # # Check if the new csv is equal to the original one
    # if n_complexes == 100:
    #     df_orig = pd.read_csv(f'../../../test/machine_learning/complex_property_prediction/input_data/test_complexes_with_descriptors_v{db_version}_orig.csv', index_col=0)
    #     df_orig = df_orig.drop(columns=['is_test_set'])
    #     common_columns = list(set(df.columns).intersection(set(df_orig.columns)))
    #     pd.testing.assert_frame_equal(df[common_columns], df_orig[common_columns], check_like=True)
    #     print('Dataframe is equal to the original one.')
    #     if not len(common_columns) == len(set(df.columns.tolist()+df_orig.columns.tolist())):
    #         print(f'But there are {len(set(df.columns.tolist()+df_orig.columns.tolist()))-len(common_columns)} not common columns and {len(common_columns)} common columns.')

    print('Done!')










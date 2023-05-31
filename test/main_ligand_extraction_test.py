"""
This is the main script for the extraction of ligands.py from a database.
"""
from copy import deepcopy
from typing import Union

import numpy as np
from sklearn.metrics import accuracy_score, r2_score

from src01.DataBase import LigandDB, ComplexDB
from src01.io_custom import load_json
from src01.main_ligand_extraction import main
from src01.utilities import unroll_dict_into_columns, sort_dict_recursively_inplace
from pathlib import Path
import pandas as pd
from src02_ChargeAssignment.charge_benchmark.merge_benchmark_charges import \
    update_ligands_with_information_from_ligand_db

class ChargeBenchmark:

    def __init__(self, true_charge_name: str):
        self.true_charge_name = true_charge_name

        self.benchmark_charge_dir = '../test/debug/databases/charge_benchmark'
        self.benchmark_charge_filenames = {
            # 'Cian1': 'Cian_already_assigned_ligand_charges.csv',
            # 'Cian2': 'Cian_BenchMark_corrected_170223.csv',
            # 'Manting': 'Manting_Corrected_170223.csv',
            # 'Marconi1': 'Marconi_corrected_by_Cian_020223.csv',
            'Cian3': 'Cian_ligand_charges_CSD_060423.csv'
        }
        self.update_properties = ['unique_name', 'name', 'graph_hash', 'local_elements', 'pred_charge', 'pred_charge_is_confident']

    def merge_benchmark_datasets(self, latest_full_ligand_db_path: str):
        """
        Update the resulting merged dfs with ligand information from the latest full ligand db. Ligands are identified by 'CSD_code' and 'stoichiometry' (hard coded).
        """
        df_all = pd.DataFrame()
        important_columns = ['CSD_code', 'stoichiometry', 'metal', 'charge', 'denticity', 'issue_detected', 'confidence', 'comment', 'author']
        for name, filename in self.benchmark_charge_filenames.items():
            df0 = pd.read_csv(Path(self.benchmark_charge_dir, filename))
            df0['author'] = name
            df0 = df0[important_columns]
            df_all = pd.concat((df_all, df0), axis=0)
        df_all = df_all.reset_index(drop=True)

        self.df_all = df_all[df_all['charge'].notna()]
        assert self.df_all.drop(
            columns=['issue_detected', 'comment']).notna().all().all(), f'Any of the columns {self.df_all.columns} is NaN which should not be NaN'

        self.df_all.loc[self.df_all['issue_detected'] == '-', 'issue_detected'] = np.nan
        self.df_all.loc[self.df_all['comment'] == '-', 'comment'] = np.nan

        # Add ligand information from the ligand db to each ligand in the charge benchmark csv.
        self.df_all = update_ligands_with_information_from_ligand_db(
            df_benchmark=self.df_all,
            latest_full_ligand_db_path=latest_full_ligand_db_path,
            update_properties=self.update_properties
        )
        if 'pred_charge' in self.df_all:
            self.df_all['prediction_error'] = self.df_all['charge'] - self.df_all['pred_charge']

        self.df_all['high_confidence'] = (self.df_all['confidence'] == 3) & self.df_all['issue_detected'].isna() & self.df_all['unique_name'].notna()
        self.df_confident = self.df_all[self.df_all['high_confidence']]
        self.df_confident = self.df_confident.drop(columns='high_confidence')

        # n_duplicates = self.df_all.loc[self.df_all['graph_hash'].notna(), 'graph_hash'].duplicated().sum()
        # if n_duplicates > 0:
        #     print(
        #         f'WARNING: {n_duplicates} duplicates of graph hashes found in the data! They are not excluded as of now.')

        return self.df_all

    def calculate_scores_of_charge_benchmark(self, full_ligand_db: Union[str, Path], expected_charge_interval=(-4, 1)):

        self.merge_benchmark_datasets(latest_full_ligand_db_path=full_ligand_db)

        delimiter = '============================'
        print(f'\nSCORES:\n{delimiter}')

        try:
            frac_not_scale_invariant = sum(~self.df_confident['charge_scale_invariant']) / len(self.df_confident)
            print(f'Frac. not scale invariant: {frac_not_scale_invariant:.2g}')
        except KeyError:
            pass

        y_pred = deepcopy(self.df_confident['pred_charge'])
        y_true = deepcopy(self.df_confident[self.true_charge_name])
        valid = y_true.notna() & y_pred.notna()

        y_true_all, y_pred_all = y_true[valid], y_pred[valid]
        acc = accuracy_score(y_true=y_true_all, y_pred=y_pred_all)
        print(f'Total accuracy (n={len(y_true_all)}): {acc:.2g}')

        confident = valid & self.df_confident['pred_charge_is_confident']
        y_true_conf, y_pred_conf = y_true[confident], y_pred[confident]
        acc = accuracy_score(y_true=y_true_conf, y_pred=y_pred_conf)
        print(f'Confident accuracy (n={len(y_true_conf)}): {acc:.2g}')

        y_true_no_conf, y_pred_no_conf = y_true[~confident], y_pred[~confident]
        acc = accuracy_score(y_true=y_true_no_conf, y_pred=y_pred_no_conf)
        print(f'Non-confident accuracy (n={len(y_true_no_conf)}): {acc:.2g}')

        if len(confident) > 0:
            print(f'Frac confident predictions: {sum(confident) / len(confident):.2g}')
        else:
            print(f'Frac confident predictions: 0.0')

        print(delimiter + '\n')

        return



# TODO in preprocessing CSD:
#   - filters:
#       - if all elements of smiles and xyz match and also with formula in api
#       - recognize valid counter ions/ solvent molecules vs graph errors in ligands.py by setting a cutoff at 3A (or so), if two graph fragments are closer than that it's treated as a graph error
#       - if ligand has C but no H (or maybe if ligand has only C)
#       - more than one metal from d or f block
#   - add properties:
#   - charges:
#
# TODO in pipeline
#   - check data pipeline from initial tmQMg up until input json
#   - add properties
#       - unique_ligand:
#           - minimum sphere
#           - min/ max interatomic distances
#           - graph_origin = graph_creation_method  ???
#   - rename
#   - unique_ligand filter:
#   - unique ligand db
#   - charge assignment:
#   - metrics:




if __name__ == '__main__':
    # example databases, choose between: tmqm, tmqmG, CSD_MM_G
    database_path = '../data_input/CSD_MM_G'  # in github
    data_store_path = '../data_output/CSD_MM_G_Jsons_test'  # directory where we want to store the jsons

    testing = 100  # if we would like to only do a test run. Set to False for full run
    graph_strategy = 'CSD'  # the desired graph strategy: default, ase_cutoff, CSD, pymatgen_NN, molsimplifyGraphs

    overwrite_atomic_properties = False     # if atomic properties json should be overwritten. Only necessary after changing input files.
    use_existing_input_json = True          # if the existing input json should be used. For speeding up test runs.
    store_database_in_memory = True        # if the database should be stored in memory. Only use if you have enough RAM, but can speed up the pipeline by maybe 30%.

    # Input complex filters
    exclude_not_fully_connected_complexes = False   # only keep complexes which are fully connected
    exclude_charged_complexes = False               # Keep only input complexes with charge of 0

    db = main(
        database_path_=database_path,
        data_store_path_=data_store_path,
        overwrite_atomic_properties_=overwrite_atomic_properties,
        use_existing_input_json_=use_existing_input_json,
        exclude_not_fully_connected_complexes_=exclude_not_fully_connected_complexes,
        testing_=testing,
        graph_strat_=graph_strategy,
        exclude_charged_complexes=exclude_charged_complexes,
        store_database_in_memory=store_database_in_memory
    )
    




    ###########################################
    #%%               DEBUGGING
    ###########################################
    if not testing == False:

        print('\nDEBUGGING:')
        reduce_to_intersection_of_rows = []
        original_suffix = f'_original_{testing}.json'

        #%%
        charge_benchmark = ChargeBenchmark(true_charge_name='charge')
        charge_benchmark.calculate_scores_of_charge_benchmark(db.full_ligands_json)
        # charge_benchmark.calculate_scores_of_charge_benchmark('/Users/timosommer/PhD/projects/RCA/projects/CreateTMC/data/final_db_versions/full_ligand_db_v1.6.json')
        #%%

        print('Read in output to look at it.')
        df_unique_ligands = pd.DataFrame.from_dict(load_json(db.unique_ligands_json), orient='index')
        df_unique_ligands = unroll_dict_into_columns(df_unique_ligands, dict_col='global_props', prefix='gbl_', delete_dict=True)
        df_unique_ligands = unroll_dict_into_columns(df_unique_ligands, dict_col='stats', prefix='stats_', delete_dict=True)
        df_full_ligands = pd.DataFrame.from_dict(load_json(db.full_ligands_json), orient='index')
        df_full_ligands = unroll_dict_into_columns(df_full_ligands, dict_col='global_props', prefix='gbl_', delete_dict=True)
        df_full_ligands = unroll_dict_into_columns(df_full_ligands, dict_col='stats', prefix='stats_', delete_dict=True)
        df_complexes = pd.DataFrame.from_dict(load_json(db.output_complexes_json), orient='index')
        df_complexes = unroll_dict_into_columns(df_complexes, dict_col='global_props', prefix='gbl_', delete_dict=True)


        c = df_complexes.iloc[0].to_dict()
        ulig = df_unique_ligands.iloc[0].to_dict()
        lig = df_full_ligands.iloc[0].to_dict()

        print('Double checking if all data is still the same after refactoring:')
        check_db = {
            'tmQM_Ligands_unique': df_unique_ligands,
            'tmQM_Ligands_full': df_full_ligands,
            'complex_db': df_complexes,
        }
        for db_name, df_new in check_db.items():

            old_path = Path(data_store_path, db_name + original_suffix)
            if old_path.exists():
                print(f'Check {db_name}:')
            else:
                print(f'ERROR: Path for {db_name} doesn\'t exist. Cannot doublecheck output.')
                continue

            df_old = pd.DataFrame.from_dict(load_json(old_path), orient='index')
            df_old = unroll_dict_into_columns(df_old, dict_col='global_props', prefix='gbl_', delete_dict=True)
            try:
                df_old = unroll_dict_into_columns(df_old, dict_col='stats', prefix='stats_', delete_dict=True)
            except KeyError:
                pass

            df_old.sort_index(inplace=True)
            df_new.sort_index(inplace=True)

            if reduce_to_intersection_of_rows:
                for col in reduce_to_intersection_of_rows:
                    intersect = set(df_old[col]).intersection(set(df_new[col]))
                    df_old = df_old[df_old[col].isin(intersect)]
                    df_new = df_new[df_new[col].isin(intersect)]

            if db_name == 'tmQMG':
                if not df_old['graph_dict'].equals(df_new['graph_dict']):
                    print('Column `graph_dict` is not equal, sort dictionaries to try to make it equal.')
                    df_old['graph_dict'] = df_old['graph_dict'].apply(sort_dict_recursively_inplace)
                    df_new['graph_dict'] = df_new['graph_dict'].apply(sort_dict_recursively_inplace)

            try:
                pd.testing.assert_frame_equal(df_new, df_old, check_like=True)
                print('Successful refactoring. All data is still the same.')

            except AssertionError:
                drop_cols = list(set(df_new.columns).symmetric_difference(set(df_old.columns)))
                drop_cols = drop_cols + ['ligands']
                print(f'Failed testing whole df. Check again without {drop_cols}.')
                pd.testing.assert_frame_equal(df_new.drop(columns=drop_cols, errors='ignore'),
                                              df_old.drop(columns=drop_cols, errors='ignore'), check_like=True)
                print(f'Mostly successful refactoring. All data is still the same when excluding {drop_cols}.')

        #%% Checks if db is still the same after writing to and re-reading from json.
        # Slow, should only be used with testing <= 100.
        print(f'\nCheck if db is equal after reading in from json:')
        same_ulig_db = db.unique_ligand_db.check_db_equal(db=db.unique_ligands_json)
        if not same_ulig_db:
            print('  Unique ligand db: not the same!')
        else:
            print('  Unique ligand db: good')

        same_full_lig_db = db.full_ligand_db.check_db_equal(db=db.full_ligands_json)
        if not same_full_lig_db:
            print('  Full ligand db: not the same!')
        else:
            print('  Full ligand db: good')

        same_complex_db = db.complex_db.check_db_equal(db=db.output_complexes_json)
        if not same_complex_db:
            print('  Complex db: not the same!')
        else:
            print('  Complex db: good')


    print('Done!')




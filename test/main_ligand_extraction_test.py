"""
This is the main script for the extraction of ligands.py from a database.
"""
from src01.io_custom import load_json, load_full_ligand_db
from src01.main_ligand_extraction import main
from src01.utilities import unroll_dict_into_columns, sort_dict_recursively_inplace
from pathlib import Path
import pandas as pd
from test.Charge_Benchmark import ChargeBenchmark



# TODO in preprocessing CSD:
#   - filters:
#       - if all elements of smiles and xyz match and also with formula in api
#       - recognize valid counter ions/ solvent molecules vs graph errors in ligands.py by setting a cutoff at 3A (or so), if two graph fragments are closer than that it's treated as a graph error
#       - if ligand has C but no H (or maybe if ligand has only C)
#       - more than one metal from d or f block
#   - add properties:
#       - uligs:
#           - has_warnings
#           - warning if ligand has bad bonds
#           - profile:
#               - fraction of homoleptic complexes
#               - list of other donor atoms
#   - charges:
#
# TODO in pipeline
#   - check data pipeline from initial tmQMg up until input json
#   - add properties
#       - unique_ligand:
#           - minimum sphere
#           - min/ max interatomic distances
#           - graph_origin = graph_creation_method  ???
#           - names of originating ligands
#           - IUPAC name of ligand based on complex names
#   - rename
#   - unique_ligand filter:
#   - unique ligand db
#   - charge assignment:
#   - metrics:




if __name__ == '__main__':
    # example databases, choose between: tmqm, tmqmG, CSD_MM_G
    database_path = '../data_input/CSD_MM_G'  # in github
    data_store_path = '../data_output/CSD_MM_G_Jsons_test'  # directory where we want to store the jsons

    testing = 5000  # if we would like to only do a test run. Set to False for full run
    graph_strategy = 'CSD'  # the desired graph strategy: default, ase_cutoff, CSD

    overwrite_atomic_properties = False     # if atomic properties json should be overwritten. Only necessary after changing input files.
    use_existing_input_json = True          # if the existing input json should be used. For speeding up test runs. Not critical
    store_database_in_memory = True         # if the database should be stored in memory. Only use if you have enough RAM, but can speed up the pipeline by maybe 30%.

    # Input complex filters
    exclude_not_fully_connected_complexes = False   # only keep complexes which are fully connected
    exclude_charged_complexes = False               # Keep only input complexes with charge of 0



    # Just for safety
    if testing == False:
        store_database_in_memory = False
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
        charge_benchmark.calculate_scores_of_charge_benchmark(db.output_complexes_json)

        #%%

        print('Read in output to look at it.')
        df_unique_ligands = pd.DataFrame.from_dict(load_json(db.unique_ligands_json), orient='index')
        df_unique_ligands = unroll_dict_into_columns(df_unique_ligands, dict_col='global_props', prefix='gbl_', delete_dict=True)
        df_unique_ligands = unroll_dict_into_columns(df_unique_ligands, dict_col='stats', prefix='stats_', delete_dict=True)
        df_full_ligands = pd.DataFrame.from_dict(load_full_ligand_db(db.output_complexes_json), orient='index')
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
            # 'tmQM_Ligands_full': df_full_ligands,
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
        # print(f'\nCheck if db is equal after reading in from json:')
        # same_ulig_db = db.unique_ligand_db.check_db_equal(db=db.unique_ligands_json)
        # if not same_ulig_db:
        #     print('  Unique ligand db: not the same!')
        # else:
        #     print('  Unique ligand db: good')
        #
        # # same_full_lig_db = db.full_ligand_db.check_db_equal(db=db.full_ligands_json)
        # # if not same_full_lig_db:
        # #     print('  Full ligand db: not the same!')
        # # else:
        # #     print('  Full ligand db: good')
        #
        # same_complex_db = db.complex_db.check_db_equal(db=db.output_complexes_json)
        # if not same_complex_db:
        #     print('  Complex db: not the same!')
        # else:
        #     print('  Complex db: good')


    print('Done!')




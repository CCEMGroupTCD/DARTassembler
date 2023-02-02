"""
This is the main script for the extraction of ligands from a database.
"""
from src01.ligand_extraction import LigandExtraction
from src01.utilities import unroll_dict_into_columns, sort_dict_recursively_inplace
from pathlib import Path
import pandas as pd

# TODO in preprocessing CSD:
#   - filters:
#       - if all elements of smiles and xyz match and also with formula in api
#
# TODO in pipeline
#   - check data pipeline from initial tmQMg up until input json
#       - write function to check that MultiGraphs have only one edge and can be made into simple graphs
#   - add properties
#       - unique_ligand:
#           - global_props
#               - graph_origin = graph_creation_method  ???
#   - rename
#       - complex dict
#           - total_q --> charge
#   - pre-assembly filters:
#       - chosen_denticity_fraction too small
#   - pre-ligand extraction filters:
#       - filter out complexes with H coordinated metals (ask Cian if makes sense)
#
# TODO when refactoring pipeline
#   - BUGFIX:
#       - in the folder `tmQMG_Jsons_fixed_gbl_props_test_full` the newly generated `complex_db.json` has 58429 entries, but the saved `complex_db_original.json` has only 58408 entries (21 less)


if __name__ == '__main__':
    # specify the database path
    database_path = '../database/tmQMg_fixed_gbl_props_cutoffs'         #'../database/tmQMg'
    data_store_path = "../data/tmQMG_Jsons_fixed_gbl_props_cutoffs_test"  # Folder where we want to store the jsons
    calculate_charges = True        # if you want to run charge assignment after ligand extraction
    overwrite_atomic_properties = True
    use_existing_input_json = True
    exclude_not_fully_connected_complexes = False
    get_only_unique_ligand_db_without_charges = False

    db = LigandExtraction(
                            database_path=database_path,
                            data_store_path=data_store_path,
                            exclude_not_fully_connected_complexes=exclude_not_fully_connected_complexes
                            )
    db.run_ligand_extraction(
                                calculate_charges=calculate_charges,
                                overwrite_atomic_properties=overwrite_atomic_properties,
                                use_existing_input_json=use_existing_input_json,
                                get_only_unique_ligand_db_without_charges=get_only_unique_ligand_db_without_charges
                            )


    ###########################################
    #               DEBUGGING
    ###########################################

    print('\nDEBUGGING:')

    print('Read in output to look at it.')
    df_unique_ligands = pd.read_json(Path(data_store_path, 'tmQM_Ligands_unique' + '.json'), orient='index')
    df_full_ligands = pd.read_json(Path(data_store_path, 'tmQM_Ligands_full' + '.json'), orient='index')
    df_complexes = pd.read_json(Path(data_store_path, 'complex_db.json'), orient='index')
    c = df_complexes.iloc[0].to_dict()
    ulig = df_unique_ligands.iloc[0].to_dict()
    lig = df_full_ligands.iloc[0].to_dict()

    print('Double checking if all data is still the same after refactoring:')
    check_db = {
        'tmQM_Ligands_unique': df_unique_ligands,
        'tmQM_Ligands_full': df_full_ligands,
        'complex_db': df_complexes,
        # 'tmQMG': df_tmqmg,
    }
    reduce_to_intersection_of_rows = []
    unroll_global_props = False
    original_suffix = '_original_unconnected.json'
    for db_name, df_new in check_db.items():

        old_path = Path(data_store_path, db_name + original_suffix)
        if old_path.exists():
            print(f'Check {db_name}:')
        else:
            print(f'ERROR: Path for {db_name} doesn\'t exist. Cannot doublecheck output.')
            continue

        df_old = pd.read_json(old_path, orient='index')

        df_old.sort_index(inplace=True)
        df_new.sort_index(inplace=True)

        if unroll_global_props:
            df_old = unroll_dict_into_columns(df_old, dict_col='global_props', prefix='gbl_')
            df_new = unroll_dict_into_columns(df_new, dict_col='global_props', prefix='gbl_')

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
            drop_cols = ['has_unconnected_ligands', 'unique_ligand_information', 'ligands']
            print(f'Failed testing whole df. Check again without {drop_cols}.')
            pd.testing.assert_frame_equal(df_new.drop(columns=drop_cols, errors='ignore'),
                                          df_old.drop(columns=drop_cols, errors='ignore'), check_like=True)
            print(f'Mostly successful refactoring. All data is still the same when excluding {drop_cols}.')

    print('Done!')



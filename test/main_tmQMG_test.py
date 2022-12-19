"""
todo: Write this into a clean class
"""
import pandas as pd
from copy import deepcopy
from src01.DataBase import MoleculeDB, LigandDB
import gc
import numpy as np
from pathlib import Path
from src01.main_tmQMG import unique_ligands_from_Ligand_batch_json_files, update_complex_db_with_ligands


# TODO list
#   - check data pipeline from initial tmQMg up until input json
#       - write function to check that MultiGraphs have only one edge and can be made into simple graphs


if __name__ == '__main__':
    # specify the database path
    database_path = '../database/tmQMg'
    data_store_path = "../data/tmQMG_Jsons_test"  # Folder where we want to store the jsons
    number_of_batches = 10



    # # Establish and safe the Database (in our case tmQM) as object
    # db_dict = DataLoader(database_path_=database_path).data_for_molDB
    # tmQMG_DB = MoleculeDB.from_json(json_=db_dict, type_="Molecule")
    # tmQMG_DB.to_json(path=f'{data_store_path}/tmQMG.json')
    #
    # del db_dict
    # gc.collect()
    #
    # print("safe successfull")
    #
    # #
    # total_keys = list(tmQMG_DB.db.keys())
    # delta = int(len(total_keys) / number_of_batches)
    #
    # ripped_keys = [total_keys[i * delta:(i + 1) * delta] for i in range(number_of_batches - 1)]
    # ripped_keys.append(total_keys[(number_of_batches - 1) * delta:])
    #
    # del tmQMG_DB
    # gc.collect()
    #
    # print("Now the ligand extraction should start running")
    #
    # assert len(ripped_keys) == number_of_batches
    # assert sum([len(key) for key in ripped_keys]) == len(total_keys)


    ripped_keys = [None]                    # todo remove when running from top
    number_of_batches = len(ripped_keys)    # todo remove when running from top

    for i, key_list in enumerate(ripped_keys):
        print(f"Batch {i + 1} of {number_of_batches} running:")
        tmQM_batch = MoleculeDB.from_json(json_=f'{data_store_path}/tmQMG.json', type_="Molecule",
                                          identifier_list=key_list)

        # Normalize graphs
        tmQM_batch.filter_not_fully_connected_molecules()   # important, will be made automatic
        tmQM_batch.remove_node_features_from_molecular_graphs(keep=['node_label'])
        tmQM_batch.remove_edge_features_from_molecular_graphs()
        tmQM_batch.normalize_multigraphs_into_simple_graphs()

        Ligand_batch = LigandDB.from_MoleculeDB(molDB=tmQM_batch)


        Ligand_batch.to_json(f"{data_store_path}/Lig_Batch_{i}.json")
        del Ligand_batch
        del tmQM_batch
        gc.collect()

    # We now have 380k Ligands and could load them batchwise via
    # lig1 = LigandDB.from_json(json_=f"../data/New_DB_jsons/Lig_Batch_{1}.json", type_="Ligand")
    # and then assemble all the ligands in one dict, create a ligand DB out of that and do the duplicant filtering.
    # However, that doesnt work due to memory reasons.

    unique_ligands_from_Ligand_batch_json_files(n=number_of_batches, data_store_path=data_store_path)
    save_complex_db_path = f'{data_store_path}/complex_db.json'
    update_complex_db_with_ligands(
                                    complex_json=f'{data_store_path}/tmQMG.json',
                                    ligand_json=f'{data_store_path}/tmQM_Ligands_full.json',
                                    save_complex_db_path = save_complex_db_path
                                    )

    print("All data established, move to playground")

    print('Read in output to look at it.')
    df_unique_ligands = pd.read_json(Path(data_store_path, 'tmQM_Ligands_unique' + '.json'), orient='index')
    df_full_ligands = pd.read_json(Path(data_store_path, 'tmQM_Ligands_full' + '.json'), orient='index')
    df_complexes = pd.read_json(save_complex_db_path, orient='index')
    c = df_complexes.iloc[0].to_dict()
    ulig = df_unique_ligands.iloc[0].to_dict()


    print('Double checking if all data is still the same after refactoring:')
    check_db = {
                'tmQM_Ligands_unique': df_unique_ligands,
                'tmQM_Ligands_full': df_full_ligands,
                'complex_db': df_complexes
                }
    for db_name, df_new in check_db.items():

        old_path = Path(data_store_path, db_name + '_original.json')
        if old_path.exists():
            print(f'Check {db_name}:')
        else:
            print(f'ERROR: Path for {db_name} doesn\'t exist. Cannot doublecheck output.')
            continue

        df_old = pd.read_json(old_path, orient='index')

        df_old.sort_index(inplace=True)
        df_new.sort_index(inplace=True)

        try:
            pd.testing.assert_frame_equal(df_new, df_old, check_like=True)
            print('Successful refactoring. All data is still the same.')

        except AssertionError:
            drop_cols = []
            print(f'Failed testing whole df. Check again without {drop_cols}.')
            pd.testing.assert_frame_equal(df_new.drop(columns=drop_cols, errors='ignore'), df_old.drop(columns=drop_cols, errors='ignore'), check_like=True)
            print(f'Mostly successful refactoring. All data is still the same when excluding {drop_cols}.')

    print('Done!')

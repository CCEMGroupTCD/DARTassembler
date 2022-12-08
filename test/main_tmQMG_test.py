"""
todo: Write this into a clean class
"""
import pandas as pd

from src01.DataBase import MoleculeDB, LigandDB
from src01.DataLoader import DataLoader
import gc
import json
from src01.Molecule import RCA_Ligand
from tqdm import tqdm
import numpy as np
import random
from pathlib import Path





def unique_ligands_from_Ligand_batch_json_files(n=10):
    gc.collect()

    # first we generate the full ligand dict
    ligand_dict = {}
    # expected_length = 0
    for i in range(n):
        with open(f"{data_store_path}/Lig_Batch_{i}.json", "r") as handle:
            dict_ = json.load(handle)
            ligand_dict.update(dict_)
            # expected_length += len(dict_)

    name_ghash_dict = {}
    # we append the dict piecewise by the graph hash
    for lig_key, lig_dict in tqdm(ligand_dict.items(), desc="Extracting graph hashs"):
        # new_ligand_dict[lig_key] = lig_dict["graph_hash"]
        name_ghash_dict[lig_key] = lig_dict["graph_hash"]

    grouped_ligands_by_hash = {}
    for k, v in tqdm(name_ghash_dict.items(), desc="grouping ligands by hash"):
        grouped_ligands_by_hash.setdefault(v, []).append(k)

    unique_ligand_dict = {}
    for same_ligands in tqdm(grouped_ligands_by_hash.values(), desc="filling unique ligand dict"):
        unique_ligand = same_ligands[0]
        unique_ligand_name = 'unq_' + unique_ligand
        for ligand_name in same_ligands:
            ligand_dict[ligand_name]["unique_ligand_name"] = unique_ligand_name

        unique_ligand_dict[unique_ligand] = ligand_dict[unique_ligand]

    #
    # Now we can dump ligand_dict to json, which is a huge dict containing ALL ligands with their unique name
    # and also unique_ligand_dict which only contains unique ligands
    with open(f"{data_store_path}/tmQM_Ligands_full.json", "w+") as file:
        json.dump(ligand_dict, file)
        print(f"Full ligand database saved to {f'{data_store_path}/tmQM_Ligands_full.json'}")
    with open(f"{data_store_path}/tmQM_Ligands_unique.json", "w+") as file:
        json.dump(unique_ligand_dict, file)
        print(f"Unique ligand database saved to {f'{data_store_path}/tmQM_Ligands_unique.json'}")

    return


if __name__ == '__main__':
    # specify the database path
    database_path = '../database/tmQMg'
    data_store_path = "../data/tmQMG_Jsons_test"  # Folder where we want to store the jsons
    number_of_batches = 10

    np.random.seed(0)
    random.seed(0)


    # Establish and safe the Database (in our case tmQM) as object
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

    ripped_keys = [None]
    number_of_batches = len(ripped_keys)

    for i, key_list in enumerate(ripped_keys):
        print(f"Batch {i + 1} of {number_of_batches} running:")
        tmQM_batch = MoleculeDB.from_json(json_=f'{data_store_path}/tmQMG.json', type_="Molecule",
                                          identifier_list=key_list)

        Ligand_batch = LigandDB.from_MoleculeDB(molDB=tmQM_batch)
        df_ligs = pd.DataFrame.from_dict(Ligand_batch.get_dict_in_json_format(), orient='index')
        df_ligs.to_csv(Path(data_store_path, 'tmQM_Ligands_int.csv'))


        Ligand_batch.to_json(f"{data_store_path}/Lig_Batch_{i}.json")
        del Ligand_batch
        del tmQM_batch
        gc.collect()
    #
    # We now have 380k Ligands and could load them batchwise via
    # lig1 = LigandDB.from_json(json_=f"../data/New_DB_jsons/Lig_Batch_{1}.json", type_="Ligand")
    # and then assemble all the ligands in one dict, create a ligand DB out of that and do the duplicant filtering.
    # However, that doesnt work due to memory reasons.

    unique_ligands_from_Ligand_batch_json_files(n=number_of_batches)

    ligand_filenames = ['tmQM_Ligands_full', 'tmQM_Ligands_unique']
    for filename in ligand_filenames:
        in_filepath = Path(data_store_path, filename + '.json')
        df_ligands = pd.read_json(in_filepath, orient='index')

        out_filepath = Path(data_store_path, filename + '.csv')
        df_ligands.to_csv(out_filepath, index=False)


    print("All data established, move to playground")


    print('Doublechecking if all data is still the same after refactoring.')


    rtol = 1e-4
    atol = 1e-6
    for db_type in ['full']:
        index = 'name' if db_type == 'full' else 'graph_hash'
        df_old = pd.read_csv(f'../data/tmQMG_Jsons_test/tmQM_Ligands_{db_type}_original.csv', index_col=index).sort_index()
        df_new = pd.read_csv(f'../data/tmQMG_Jsons_test/tmQM_Ligands_{db_type}.csv', index_col=index).sort_index()
        try:
            pd.testing.assert_frame_equal(df_new, df_old, check_like=True, rtol=rtol, atol=atol)
            print('Successful refactoring. All data is still the same.')
        except AssertionError:
            drop_cols = ['graph_dict']#['unique_ligand_name', 'atomic_props', 'global_props', 'graph_dict', 'ligand_to_metal', 'denticity', 'graph_hash']
            print(f'Failed testing whole df. Check again without {drop_cols}.')
            pd.testing.assert_frame_equal(df_new.drop(columns=drop_cols), df_old.drop(columns=drop_cols), check_like=True, rtol=rtol, atol=atol)
            print(f'Mostly successful refactoring. All data is still the same when excluding {drop_cols}.')

    check_prop = 'ligand_to_metal'
    df = df_new.join(df_old, rsuffix='_new', lsuffix='_old')
    df = df[[f'{check_prop}_new', f'{check_prop}_old']]
    df['diff'] = df[f'{check_prop}_new'] != df[f'{check_prop}_old']


    # check_jsons = [
    #                 {
    #                     'old': '../data/tmQMG_Jsons_test/tmQM_Ligands_full_original.json',
    #                     'new': '../data/tmQMG_Jsons_test/tmQM_Ligands_full.json',
    #                 },
    #                 {
    #                     'old': '../data/tmQMG_Jsons_test/tmQM_Ligands_unique_original.json',
    #                     'new': '../data/tmQMG_Jsons_test/tmQM_Ligands_unique.json'
    #                 }
    #             ]
    # for old_new_json in check_jsons:
    #     print(f'Checking {old_new_json}.')
    #
    #     with open(old_new_json['old'], 'r') as old_file:
    #         old_d = json.load(old_file)
    #     with open(old_new_json['new'], 'r') as new_file:
    #         new_d = json.load(new_file)
    #
    #     if not len(old_d) == len(new_d):
    #         print(f'WARNING: size of dictionaries differs: {len(old_d), len(new_d)}')
    #
    #     for key, old_value in old_d.items():
    #         if not key in new_d:
    #             print(f'WARNING: old key not in new dict: {key}')
    #         else:
    #             new_value = new_d[key]
    #             if not old_value == new_value:
    #                 print(f'WARNING: Values differ for key {key}: {old_value, new_value}')

    print('Done!')

from src01.DataBase import MoleculeDB, LigandDB
from src01.DataLoader import DataLoader
import gc
from src01.graph_utility import *
from src01.Molecule import RCA_Ligand
from tqdm import tqdm
import json
import pandas as pd


data_store_path = "../data/New_DB_jsons"


def extract_ligands_piecewise(n=10):
    tmQM_DB = MoleculeDB.from_json(json_='../data/New_DB_jsons/tmQM.json',
                                   type_="Molecule",
                                   max_number=None
                                   )

    total_keys = list(tmQM_DB.db.keys())
    delta = int(len(total_keys) / n)

    ripped_keys = [total_keys[i * delta:(i + 1) * delta] for i in range(n-1)]
    ripped_keys.append(total_keys[(n-1) * delta:])

    del tmQM_DB
    gc.collect()

    assert len(ripped_keys) == n
    assert sum([len(key) for key in ripped_keys]) == len(total_keys)

    for i, key_list in enumerate(ripped_keys):
        print(f"Batch {i+1} of {n} running:")
        tmQM_batch = MoleculeDB.from_json(json_='../data/New_DB_jsons/tmQM.json', type_="Molecule", identifier_list=key_list)

        Ligand_batch = LigandDB.from_MoleculeDB(molDB=tmQM_batch)

        Ligand_batch.to_json(f"../data/New_DB_jsons/Lig_Batch_{i}.json")
        del Ligand_batch
        del tmQM_batch
        gc.collect()


def unique_ligands_from_Ligand_batch_json_files(n=10):
    import json

    gc.collect()

    # first we generate the full ligand dict
    ligand_dict = {}
    # expected_length = 0
    for i in range(n):
        with open(f"{data_store_path}/Lig_Batch_{i}.json", "r") as handle:
            dict_ = json.load(handle)
            ligand_dict.update(dict_)
            # expected_length += len(dict_)

    new_ligand_dict = {}
    # we append the dict piecewise by the graph hash
    for lig_key, lig_dict in tqdm(ligand_dict.items(), desc="Generating graph hashs"):
        lig = RCA_Ligand.read_from_mol_dict(lig_dict)
        new_ligand_dict[lig_key] = lig.graph_hash

    # now we the dict appended by graph hashs
    name_ghash_dict = {lig_key: item for lig_key, item in new_ligand_dict.items()}

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
    with open(f"{data_store_path}/tmQM_Ligands_unique.json", "w+") as file:
        json.dump(unique_ligand_dict, file)
        #
        #
    # das jeweilige koennen wir dann im Playground einlesen und sind wieder da wo wir raus
    # kommen wollten

    print("done")

    return


def old():
    ## Establish and safe the Database (in our case tmQM) as object
    """
    database_path = '../database/tmQM/data'
    tmQM_DB = MoleculeDB.from_json(json_=DataLoader(database_path_=database_path).data_for_molDB,
                                   type_="Molecule",
                                   max_number=None
                                   )

    tmQM_DB.to_json(path='../data/New_DB_jsons/tmQM.json')
    """
    # reset memory

    # \etract the ligands batchwise
    # extract_ligands_piecewise(n=10)
    #
    #
    #
    #
    # specify the database path

    # unique_ligands_from_Ligand_batch_json_files(n=10)


def analysis_for_LinEQSolver(unique_ligands_dict, full_ligands_dict, selected_df):

    csd_codes = [selected_df.loc[i, "CSD_code"] for i in selected_df.index]

    unique_ligand_names = [value['unique_ligand_name'] for value in unique_ligands_dict.values()]

    hits, misses = 0, 0

    for unique_ligand in tqdm(unique_ligand_names):

        for lig_name, value in full_ligands_dict.items():
            if value['unique_ligand_name'] == unique_ligand:
                if lig_name.split("-")[1] in csd_codes:
                    hits += 1
                    break

        misses += 1

    return hits, misses


if __name__ == '__main__':

    with open(f"{data_store_path}/tmQM_Ligands_unique.json") as file:
        unique_ligands_dict = json.load(file)

    with open(f"{data_store_path}/tmQM_Ligands_full.json") as file:
        full_ligands_dict = json.load(file)


    df = pd.read_csv("../src04_LigandAnalysis/CSD.csv")
    # lig1 = LigandDB.from_json(json_=f"../data/New_DB_jsons/Lig_Batch_{1}.json", type_="Ligand")

    selected_df = df[~df['metal_nr_if_exists'].isnull()]
    number_equations = len(selected_df)

    hits, misses = analysis_for_LinEQSolver(unique_ligands_dict, full_ligands_dict, selected_df)

    print(f"hits: {hits}, misses: {misses}")

    print("done")

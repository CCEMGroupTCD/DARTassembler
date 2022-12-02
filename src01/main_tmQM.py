"""
todo: Write this into a clean class
"""

from src01.DataBase import MoleculeDB, LigandDB
from src01.DataLoader import DataLoader
import gc
import json
from src01.Molecule import RCA_Ligand
from tqdm import tqdm

#data_store_path = "../data/tmQM_Jsons"      # Folder where we want to store the jsons
data_store_path = "../data/New_DB_jsons"


def extract_ligands_piecewise(n=10):
    tmQM_DB = MoleculeDB.from_json(json_=f'{data_store_path}/tmQM.json',
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
        tmQM_batch = MoleculeDB.from_json(json_=f'{data_store_path}/tmQM.json', type_="Molecule", identifier_list=key_list)

        Ligand_batch = LigandDB.from_MoleculeDB(molDB=tmQM_batch)

        Ligand_batch.to_json(f"{data_store_path}/Lig_Batch_{i}.json")
        del Ligand_batch
        del tmQM_batch
        gc.collect()


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
        #new_ligand_dict[lig_key] = lig_dict["graph_hash"]
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
        print(f"Unique ligand database saved to {f'{data_store_path}/tmQM_Ligands_full.json'}")

    return


if __name__ == '__main__':

    """
    # specify the database path
    database_path = '../database/tmQM/data'
    number_of_batches = 10

    # Establish and safe the Database (in our case tmQM) as object
    tmQM_DB = MoleculeDB.from_json(json_=DataLoader(database_path_=database_path).data_for_molDB, type_="Molecule")
    tmQM_DB.to_json(path='../data/New_DB_jsons/tmQM.json')

    #
    extract_ligands_piecewise(n=number_of_batches)
    #
    # We now have 380k Ligands and could load them batchwise via
    # lig1 = LigandDB.from_json(json_=f"../data/New_DB_jsons/Lig_Batch_{1}.json", type_="Ligand")
    # and then assemble all the ligands in one dict, create a ligand DB out of that and do the duplicant filtering.
    # However, that doesnt work due to memory reasons.
    """
    number_of_batches = 10
    unique_ligands_from_Ligand_batch_json_files(n=number_of_batches)

    print("All data established, move to playground")

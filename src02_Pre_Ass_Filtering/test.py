# only a testfile to crunch some numbers for the group

from src.process import LigandDatabase
from src02_Pre_Ass_Filtering.FilteringStage import FilterHandler
import pickle
import numpy as np


if __name__ == "__main__":
    with open("../data/ligand_db_NO_filtered.pickle", "rb") as handle:
        ligand_db = pickle.load(handle)

    new_ligand_dict = {}

    for dent, ligand_list in ligand_db.NO_filtered_ligand_dict.items():
        print(f"dentcity {dent}")
        print(f"Pre-filtering: {len(ligand_list)}")
        # duplicate filtering
        new_ligand_list = list(set(ligand_list))
        new_ligand_dict[dent] = [lig for lig in new_ligand_list if lig.betaH_check() is False]
        print(f"Post-filtering: {len(new_ligand_dict[dent])}")
        print("\n")

    with open("../data/ligand_dict_NO_and_dup_and_betaH_filtered.pickle", "wb") as handle:
        pickle.dump(new_ligand_dict, handle)

    print("done")


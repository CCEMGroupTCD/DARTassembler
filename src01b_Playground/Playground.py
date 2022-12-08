"""
I always needed a file where, which could simply be run in the debugger
mode to get quick acess to the tmQM and its ligands
in the RCA_Molecule and RCA_Ligand format respecitvely
and here it is
"""
from src01.DataBase import MoleculeDB, LigandDB
import json
from tqdm import tqdm
import matplotlib.pyplot as plt


if __name__ == "__main__":

    data_path = "../data/tmQMG_Jsons"

    # Playground for tmQMG

    # First unique ligands
    #with open(f"{data_store_path}/tmQMG_Ligands_unique.json") as file:
    #    unique_ligands_dict = json.load(file)

    #lig_dict_old_format = {}
    #for k, v in tqdm(unique_ligands_dict.items()):
    #    lig_dict_old_format.setdefault(v['denticity'], []).append(k)


    #
    #
    #plt_dict = {k: len(v) for k, v in lig_dict_old_format.items()}
    #plt.bar(plt_dict.keys(), plt_dict.values(), color="r")
    #plt.ylabel("Number of Ligands")
    #plt.xlabel("Denticity")
    #plt.title("Denticity distribution")
    #plt.xticks(range(1, 11))
    #plt.yticks(range(0, 12000, 1000))
    #plt.grid()
    #plt.show()

    #
    #tmQM_DB = MoleculeDB.from_json(json_=f'{data_path}/tmQMG.json', type_="Molecule")

    # Create the LigandDB from the tmQM
    #tmQM_Ligands = LigandDB.from_json(json_=f'{data_path}/tmQMG_Ligands_full.json', type_="Ligand")

    tmQM_unique_Ligands = LigandDB.from_json(json_='../data/New_DB_jsons/tmQM_ligands_unique.json', type_="Ligand")



    #with open(f"{data_store_path}/tmQMG_Ligands_full.json") as file:
    #    full_ligands_dict = json.load(file)

    print("Playground established")

    print("done")

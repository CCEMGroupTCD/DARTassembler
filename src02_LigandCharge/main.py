from src01.DataBase import MoleculeDB, LigandDB
import json
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd


if __name__ == "__main__":

    data_store_path = "../data/tmQMG_Jsons"

    # Playground for tmQMG

    # First unique ligands
    with open(f"{data_store_path}/tmQMG_Ligands_full.json") as file:
        ligands_dict = json.load(file)

    with open(f"{data_store_path}/tmQMG.json") as file:
        complex_dict = json.load(file)

    #
    df = pd.read_csv("../src02_LigandCharge/CSD.csv")
    #
    # get subdataframe for all metals where MOXS exist
    df = df[df["metal_str_if_exists"].notna()]

    new_dict = {}
    for i in tqdm(df.index):
        csd_code = df.loc[i, "CSD_code"]
        if csd_code in complex_dict:
            d = {}

            metal_sym = df.loc[i, "metal_name"].split(")")[0]
            d["metal"] = metal_sym
            d["metal_oxi_state"] = df.loc[i, "metal_nr_if_exists"]

            d["total_q"] = complex_dict[csd_code]["global_props"]["charge"]

            metal_ind = complex_dict[csd_code]["atomic_props"]["atoms"].index(metal_sym)
            d['Metal_q'] = complex_dict[csd_code]["atomic_props"]["partial_charge"][metal_ind]

            d["ligands"] = []

            for lig_name, lig_dict in ligands_dict.items():
                if lig_name.split("-")[1] == csd_code:
                    ld = {}
                    ld['denticity'] = lig_dict['denticity']
                    ld["name"] = lig_name
                    ld["graph_hash"] = lig_dict["graph_hash"]
                    ld["unique_name"] = lig_dict["unique_ligand_name"]
                    ld["donor_atom_indices"] = lig_dict["ligand_to_metal"]

                    lig_ap = lig_dict["atomic_props"]
                    partial_charge_dict = {i: [a, [pc]] for i, (a, pc) in enumerate(zip(lig_ap["atoms"], lig_ap["partial_charge"]))}
                    coord_dict = {i: [a, [x, y, z]] for i, (a, x, y,z ) in enumerate(zip(lig_ap["atoms"], lig_ap["x"], lig_ap["y"], lig_ap["z"]))}

                    ld["atomic_props"] = {
                        "partial_charge": partial_charge_dict,
                        "coordinates": coord_dict
                    }

                    d['ligands'].append(ld)
                    # maybe for efficiency reasons?
                    #del ligands_dict[lig_name]

            #
            new_dict[csd_code] = d
            # maybe for efficiency reasons?
            #del complex_dict[csd_code]

    with open("../src02_LigandCharge/test.json", "w+") as file:
        json.dump(new_dict, file)

    print("done")

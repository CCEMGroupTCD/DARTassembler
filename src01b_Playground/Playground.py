"""
I always needed a file where, which could simply be run in the debugger
mode to get quick acess to the tmQM and its ligands
in the RCA_Molecule and RCA_Ligand format respecitvely
and here it is
"""
from src01.DataBase import MoleculeDB
from src01.DataLoader import DataLoader


if __name__ == "__main__":

    """
    data_path = "../data_output/tmQMG_Jsons"

   

    #
    tmQM_DB = MoleculeDB.from_json(json_=f'{data_path}/tmQMG.json', type_="Molecule", identifier_list=["NIBTAT"])

    # Create the LigandDB from the tmQM
    tmQM_Ligands = LigandDB.from_json(json_=f'{data_path}/tmQMG_Ligands_full.json', type_="Ligand", identifier_list=["CSD-NIBTAT-03-a"])

    lig = list(tmQM_Ligands.db.values())[0]
    mol = list(tmQM_DB.db.values())[0]

    tmQM_unique_Ligands = LigandDB.from_json(json_='../data/New_DB_jsons/tmQM_ligands_unique.json', type_="Ligand")



    #with open(f"{data_store_path}/tmQMG_Ligands_full.json") as file:
    #    full_ligands_dict = json.load(file)

    print("Playground established")

    print("done")
    """
    database_path = '../database/tmQM'

    tmQM_DB = MoleculeDB.from_json(json_=DataLoader(database_path_=database_path).data_for_molDB,
                                   type_="Molecule",
                                   max_number=100,
                                   graph_strategy="smiles"
                                   )

    print("done")
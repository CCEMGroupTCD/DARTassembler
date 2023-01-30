"""
This is just for me to test some molsimplify functionilaty
"""

from molSimplify.Classes.mol3D import mol3D
from src01.DataBase import MoleculeDB
from src01.DataLoader import DataLoader
from src01.Molecule import RCA_Molecule

database_path = '../database/tmQM/data'


if __name__ == "__main__":

    tmQM_DB = MoleculeDB.from_json(json_=DataLoader(database_path_=database_path).data_for_molDB,
                                   type_="Molecule",
                                   max_number=10
                                   )

    RCA_mol = list(tmQM_DB.db.values())[0]

    molmol = RCA_mol.to_molsimplify_mol()

    print("done")

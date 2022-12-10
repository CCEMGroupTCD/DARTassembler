"""
My File to play arround with the assmebly
"""
from src01.DataBase import LigandDB
from src03_Assembly.Assembly_Test import RandomComplexAssembler

import rdkit


if __name__ == "__main__":

    F = LigandDB.from_json(json_="../data/Filtered_Jsons/filteredLigDB.json", type_="Ligand")

    ligand = list(F.db.values())[0]

    RCA = RandomComplexAssembler(database=F)
    successfull_generated_complexes = {}

    #
    for i in range(10):
        try:
            random_complex = RCA.create_random_TMC(visualize_=True,
                                                   optimize_=False)

        except rdkit.Chem.rdchem.AtomValenceException as ex:
            print(f"The standard Error: {ex} has occured. Solveable by turning sanitize in stk off")
            random_complex = None

        except Exception as e:
            print(f"Something else went wrong: {e}")
            random_complex = None

        if random_complex is None:
            print("None")
        input("press Enter to continue")

    print("done")

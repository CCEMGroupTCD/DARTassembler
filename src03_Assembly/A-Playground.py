"""
My File to play arround with the assmebly
is goingto become the main later
"""
from src01.DataBase import LigandDB
from src03_Assembly.RandomComplexAssembler import RandomComplexAssembler

import rdkit


if __name__ == "__main__":

    F = LigandDB.from_json(json_="../data/Filtered_Jsons/filteredLigDB.json", type_="Ligand")

    ligand = list(F.db.values())[0]

    RCA = RandomComplexAssembler(database=F)
    successfull_generated_complexes = {}

    # c = RCA.create_random_TMC(random_seed=4)
    # for testing, works fine

    for i in range(10):
        try:
            random_complex = RCA.create_random_TMC()


        except rdkit.Chem.rdchem.AtomValenceException as ex:
            print(f"The standard Error: {ex} has occured. Solveable by turning sanitize in stk off")
            random_complex = None

        except Exception as e:
            print(f"Something else went wrong: {e}")
            random_complex = None

        if random_complex is None:
            print("None")
        else:
            print(f"Mr or Ms {random_complex.working_name} got created! Yuhei!\n")
            successfull_generated_complexes[random_complex.name] = random_complex

    print("done")

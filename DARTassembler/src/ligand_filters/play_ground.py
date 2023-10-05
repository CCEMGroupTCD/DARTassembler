from pathlib import Path
from DARTassembler.src.constants.Paths import project_path
from DARTassembler.src.ligand_extraction.DataBase import MoleculeDB, LigandDB
import networkx as nx
from copy import deepcopy
import pandas as pd



class TEST:

    def __init__(self,
                 database: [MoleculeDB, LigandDB],
                 safe_path: [str, Path] = None):
        """
        :param database:
        :param safe_path:   The path we want to safe the filtered databases to. If None, no saving at all
        """
        # The object we are working on. Is hopefully saved, because the filter stage doesnt make copy itself of it
        self.database = database

        self.safe_path = safe_path
        self.safe = True if self.safe_path is not None else False
        self.filter_tracking = {}

    def test(self, denticity: int = None, bonds: list = None, instruction: str = None):
        # instruction = "must_include"
        # instruction = "must_exclude"
        new_db = deepcopy(self.database.db)
        i = 0
        for unq_name, ligand in self.database.db.items():
            if ligand.denticity == denticity:
                print(i)
                print("##########")
                print(ligand.atomic_props["atoms"])
                print(unq_name)
                graph = nx.Graph(ligand.graph)
                detected_bond_tracker = []
                for edge in graph.edges:
                    index_1 = edge[0]
                    index_2 = edge[1]
                    atom_1 = graph.nodes[index_1]["node_label"]
                    atom_2 = graph.nodes[index_2]["node_label"]
                    for bond in bonds:
                        if (atom_1 == bond.split("-")[0] and atom_2 == bond.split("-")[1]) or (atom_1 == bond.split("-")[1] and atom_2 == bond.split("-")[0]):
                            detected_bond_tracker.append(bond)
                        else:
                            pass
                        if sorted(list(set(detected_bond_tracker))) == sorted(list(set(deepcopy(bonds)))):
                            print("we have detected all the bonds in our input")

                            if instruction == "must_include":
                                pass
                            elif instruction == "must_exclude":
                                print("fblsd3")
                                del new_db[unq_name]
                                pass
                            break
                    if sorted(list(set(detected_bond_tracker))) == sorted(list(set(deepcopy(bonds)))):
                        break
                    else:
                        pass

                if (set(detected_bond_tracker).intersection(set(deepcopy(bonds))) and (sorted(list(set(detected_bond_tracker))) == sorted(list(set(deepcopy(bonds)))))) or (
                        detected_bond_tracker == []):
                    print("these bonds have not been found")
                    if instruction == "must_include":
                        print("dfBLIWDEFGLEIGVB")
                        del new_db[unq_name]
                        pass
                    elif instruction == "must_exclude":
                        pass
                i += 1
            else:
                pass

            self.database.db = new_db
            self.filter_tracking[len(self.filter_tracking)] = f"bond filter with bonds {bonds} and instruction {instruction}"

    def loop(self):
        results = []
        i = 0
        print(len(self.database))
        for unq_name, ligand in self.database.db.items():
            if i >= 0:
                ligand.view_3d()
                print(f"CSD_Code: {ligand.CSD_code}")
                print(f"Ligand_Charge: {ligand.pred_charge}")
                input_ = input("Is this Carbene sp2 or sp3\n")
                if (str(input_) == "3") or (str(input_) == "2") or (str(input_) == "c"):
                    results.append({"CSD_Code": str(ligand.CSD_code),
                                    "Ligand_Charge": ligand.pred_charge,
                                    "Carbon_state": str(input_)})
                else:
                    df = pd.DataFrame(results)
                    df.to_csv("file_2.csv")
                    exit()

            else:
                pass
            print(f"Counter: {i}")
            i = i + 1

    def loop_standard(self):
        stoic_string = ""
        for unq_name, ligand in self.database.db.items():
            stoic_string = stoic_string + str(ligand.stoichiometry)

        print(list(set(stoic_string)))
        for character in list(set(stoic_string)):
            try:
                int(character)
            except:
                print(character)


if __name__ == "__main__":
    """
    One easy method which includes all the custom filters one could potentially apply to a ligand
    filtered_db = {unique_lig_name: RCA_Ligand} dict-type, which can then be employed for the ligandAssembly
    """

    tmQM_unique_Ligands = LigandDB.from_json(json_=project_path().extend("data", "Filtered_Jsons", "P_N_with_Benzene_Br", "P_N_080823.json"),
                                             type_="Ligand",
                                             max_number=63053649734659453576395)
    TEST(database=tmQM_unique_Ligands).loop_standard()
    # Filter = TEST(tmQM_unique_Ligands)
    # Filter.test(denticity=2, bonds=["C-C", "C-H"], instruction="must_exclude")
    print("done")

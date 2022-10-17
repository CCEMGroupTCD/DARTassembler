from src.process import LigandDatabase
from src02_Pre_Ass_Filtering.FilteringStage import FilterHandler
import pickle


if __name__ == "__main__":

    with open("../data/ligand_db.pickle", "rb") as handle:
        ligand_db = pickle.load(handle)

    Filter = FilterHandler(ligand_db)

    print("Filter initialized")

    Filter.filter_N_and_O_functional_groups(safe_path="../data/ligand_db_NO_filtered.pickle",
                                            i=1
                                            )

    print("Done and Not Safed")

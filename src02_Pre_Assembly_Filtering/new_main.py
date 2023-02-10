from src02_Pre_Assembly_Filtering.Filter import Filter
from src01.main_ligand_extraction import select_example_database
from src01.DataBase import LigandDB
from pathlib import Path

if __name__ == "__main__":

    selected_db = "CSD_MM_G"  # from which DB the ligands should come
    _, datapath = select_example_database(DB=selected_db)

    F = Filter(datapath,
               max_number=5000
               )

    F.run_filters()
    F.add_constant_ligands()
    F.add_reactant()
    F.safe()

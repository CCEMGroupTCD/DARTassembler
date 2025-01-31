from DARTassembler.filter_ligands import filter_ligands
from pathlib import Path as pt
from DARTassembler.src.constants.Paths import project_path
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from Rotator import LigRotNorm
import pandas as pd
from DARTassembler.src.ligand_extraction.io_custom import load_json
from DARTassembler.src.ligand_extraction.utilities_graph import graph_from_graph_dict


class GenLigands:
    """
    Generate ligands from a given input file.
    i.e. GenLigands(namx=1000, input_file="/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/DART/dev/Assembler_revision_jan_2025/test.yml").get_ligands()
    """

    def __init__(self, namx: int, input_file: str):
        self.nmax = namx
        self.input_file = input_file

    def apply_filters(self) -> None:
        """
        Apply filters to the ligand database.
        :return: None
        """
        filters = filter_ligands(filter_input_path=self.input_file, nmax=self.nmax)


class Utils:
    """
    Utility functions for the project.
    """


if __name__ == '__main__':
    # Initialize the Utils class this class contains various utility functions that do not belong to any class.
    Utilities = Utils()

    # Initialize paths as pathlib objects
    mono_yml = pt("/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/DART/dev/Assembler_revision_jan_2025/ligand_filter_inputs/monodentate.yml")
    bidentate_yml = pt("/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/DART/dev/Assembler_revision_jan_2025/ligand_filter_inputs/bidentate.yml")
    pentadentate_yml = pt("/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/DART/dev/Assembler_revision_jan_2025/ligand_filter_inputs/pentadentate_haptic.yml")

    # Initialize the json database paths
    monodentate_db_path = mono_yml.parent / "monodentate.json"
    bidentate_db_path = bidentate_yml.parent / "bidentate.json"
    pentadentate_db_path = pentadentate_yml.parent / "pentadentate_haptic.json"

    # Filter the ligand database
    max_ligands = 1000
    GenLigands(namx=max_ligands, input_file=str(mono_yml)).apply_filters()
    GenLigands(namx=max_ligands, input_file=str(bidentate_yml)).apply_filters()
    GenLigands(namx=max_ligands, input_file=str(pentadentate_yml)).apply_filters()

    # Generate the db object
    monodentate_db = LigandDB.load_from_json(path=monodentate_db_path, show_progress=True).db
    bidentate_db = LigandDB.load_from_json(path=bidentate_db_path, show_progress=True).db
    pentadentate_db = LigandDB.load_from_json(path=pentadentate_db_path, show_progress=True).db


    # Loop through the ligands
    for mono_lig, bi_lig, penta_lig in zip(monodentate_db.values(), bidentate_db.values(), pentadentate_db.values()):
        LigRotNorm(mono_lig)
        LigRotNorm(bi_lig)
        LigRotNorm(penta_lig)




    # Script finished
    print("Done!")

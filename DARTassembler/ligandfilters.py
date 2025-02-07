"""
This module is a wrapper for the ligand filtering code. It takes in a filter input file and returns a ligand filter object.
"""
from DARTassembler.src.ligand_extraction.io_custom import read_yaml
from typing import Union
from pathlib import Path
from dev.DART_refactoring_to_v1_1_0.refactor_ligandfilters import NewLigandFilters

def ligandfilters(filter_input_path: Union[str, Path], nmax: Union[int, None] = None, delete_output_dir: bool = False) -> NewLigandFilters:
    """
    Filter the full ligand database according to the specified filters. Should be run before assembly to reduce the number of ligands considered in the assembly to the ones that are interesting to the user.
    :param filter_input_path: Path to the filter input file
    :param nmax: Maximum number of ligands to be read in from the initial full ligand database. If None, all ligands are read in. This is useful for testing purposes.
    :param outpath: Path to the output ligand database. If None, the output ligand database is saved in the same directory as the input filter file.
    :param delete_output_dir: If True, the output ligand db file and the output ligand info directory are deleted before the new files are written. This is useful for testing purposes.
    :return: LigandFilters object
    """
    input_dict = read_yaml(filter_input_path)
    input_db_file = input_dict.pop('input_db_file', None)

    filter = NewLigandFilters(input_db_file=input_db_file, n_max=nmax)
    filter.get_filtered_db(**input_dict, pre_delete=delete_output_dir)

    return filter




# Integration test, to check if everything is working and the same as before.
if __name__ == "__main__":
    from DARTassembler.src.constants.Paths import project_path
    from dev.test.Integration_Test import IntegrationTest

    ligand_filter_path = project_path().extend('src05_Assembly_Refactor', 'ligandfilters.yml')
    max_number = 5000

    filter = LigandFilters(filepath=ligand_filter_path, max_number=max_number)
    filter.save_filtered_ligand_db()

    # Check if the new filtered db is the same as the old one
    benchmark_dir = project_path().extend("src14_Assembly_Unit_Test", 'ligandfilters_benchmark')
    if benchmark_dir.exists():
        test = IntegrationTest(new_dir=filter.output_ligand_db_path.parent, old_dir=benchmark_dir)
        test.compare_all()
    else:
        print("No benchmark directory found. Could not perform integration test.")
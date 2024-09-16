"""
This module is a wrapper for the ligand filtering code. It takes in a filter input file and returns a ligand filter object.
"""
from DARTassembler.src.ligand_filters.Ligand_Filters import LigandFilters
from typing import Union
from pathlib import Path
import matplotlib
try:    # Avoid error when running on server
    matplotlib.use('TkAgg')
except ImportError:
    pass

def ligandfilters(filter_input_path: Union[str, Path], nmax: Union[int, None] = None, outpath: Union[str, Path] = None, delete_output_dir: bool = False) -> LigandFilters:
    """
    Filter the full ligand database according to the specified filters. Should be run before assembly to reduce the number of ligands considered in the assembly to the ones that are interesting to the user.
    :param filter_input_path: Path to the filter input file
    :param nmax: Maximum number of ligands to be read in from the initial full ligand database. If None, all ligands are read in. This is useful for testing purposes.
    :param outpath: Path to the output ligand database. If None, the output ligand database is saved in the same directory as the input filter file.
    :param delete_output_dir: If True, the output directory is deleted before the new ligand database is saved. This is useful for testing purposes.
    :return: LigandFilters object
    """
    filter = LigandFilters(filepath=filter_input_path, max_number=nmax, output_ligand_db_path=outpath, delete_output_dir=delete_output_dir)

    filter.save_filtered_ligand_db()

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
"""
This module is a wrapper for the ligand filtering code. It takes in a filter input file and returns a ligand filter object.
"""
from DARTassembler.src.ligand_filters.Ligand_Filters import LigandFilters
from typing import Union
from pathlib import Path

def filter_ligands(filter_input_path: Union[str, Path], nmax: Union[int, None] = None, outpath: Union[str, Path] = None):
    """
    Filter the full ligand database according to the specified filters. Should be run before assembly to reduce the number of ligands considered in the assembly to the the ones that are interesting to the user.
    :param filter_input_path: Path to the filter input file
    :param nmax: Maximum number of ligands to be read in from the initial full ligand database. If None, all ligands are read in. This is useful for testing purposes.
    :return: LigandFilters object
    """
    filter = LigandFilters(filepath=filter_input_path, max_number=nmax)

    # Overwrite the output path if specified
    if outpath is not None:
        filter.output_ligand_db_path = Path(outpath)

    filter.save_filtered_ligand_db()
    # filter.save_filter_tracking()
    # filter.save_ligand_info_csv()
    filter.print_filter_tracking()

    return filter




# Integration test, to check if everything is working and the same as before.
if __name__ == "__main__":
    from DARTassembler.src.constants.Paths import project_path
    from dev.test.Integration_Test import IntegrationTest

    ligand_filter_path = project_path().extend('src05_Assembly_Refactor', 'ligandfilters.yml')
    max_number = 5000

    filter = LigandFilters(filepath=ligand_filter_path, max_number=max_number)
    filter.save_filtered_ligand_db()
    filter.save_filter_tracking()
    filter.print_filter_tracking()

    # Check if the new filtered db is the same as the old one
    benchmark_dir = project_path().extend("src14_Assembly_Unit_Test", 'ligandfilters_benchmark')
    if benchmark_dir.exists():
        test = IntegrationTest(new_dir=filter.output_ligand_db_path.parent, old_dir=benchmark_dir)
        test.compare_all()
    else:
        print("No benchmark directory found. Could not perform integration test.")
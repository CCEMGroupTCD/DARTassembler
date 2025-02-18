"""
This module is a wrapper for the ligand filtering code. It takes in a filter input file and returns a ligand filter object.
"""
from DARTassembler.src.ligand_extraction.io_custom import read_yaml
from typing import Union
from pathlib import Path
from dev.DART_refactoring_to_v1_1_0.refactor_ligandfilters import NewLigandFilters

def ligandfilters(filter_input_path: Union[str, Path], pre_delete: bool = False) -> NewLigandFilters:
    """
    Filter the full ligand database according to the specified filters. Should be run before assembly to reduce the number of ligands considered in the assembly to the ones that are interesting to the user.
    :param filter_input_path: Path to the filter input file
    :param pre_delete: If True, the output ligand db file and the output ligand info directory are deleted before the new files are written. This is useful for testing purposes.
    :return: LigandFilters object
    """
    input_dict = read_yaml(filter_input_path)
    input_db_file = input_dict.pop('input_db_file', None)
    n_max_ligands = input_dict.pop('n_max_ligands', None)

    filter = NewLigandFilters(input_db_file=input_db_file, n_max_ligands=n_max_ligands)
    filter.get_filtered_db(**input_dict, pre_delete=pre_delete)

    return filter




if __name__ == "__main__":

    ligand_filter_path = '/Users/timosommer/PhD/projects/DARTassembler/dev/DART_refactoring_to_v1_1_0/data/assembler/example_ligand_db/ligandfilters_monodentate.yml'

    filter = ligandfilters(ligand_filter_path)

    # # Check if the new filtered db is the same as the old one
    # from DARTassembler.src.constants.Paths import project_path
    # from dev.test.Integration_Test import IntegrationTest
    # benchmark_dir = project_path().extend("src14_Assembly_Unit_Test", 'ligandfilters_benchmark')
    # if benchmark_dir.exists():
    #     test = IntegrationTest(new_dir=filter.output_ligand_db_path.parent, old_dir=benchmark_dir)
    #     test.compare_all()
    # else:
    #     print("No benchmark directory found. Could not perform integration test.")
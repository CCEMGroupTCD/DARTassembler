"""
This file is for testing the DART package.
"""
import DART
from constants.Paths import project_path

if __name__ == '__main__':

    assembly_input = project_path().extend('src05_Assembly_Refactor', 'assembly_input.yml')
    DART.assemble_complexes(assembly_input_path=assembly_input)

    ligand_filter_path = project_path().extend('src05_Assembly_Refactor', 'ligandfilters.yml')
    DART.filter_ligands(filter_input_path=ligand_filter_path, max_number=5000)
"""
Integration test for the assembly of complexes.
"""
from DART.assemble_complexes import assemble_complexes
from pathlib import Path
from constants.Paths import project_path

if __name__ == "__main__":

    USER_INPUT = project_path().extend('src05_Assembly_Refactor', 'assembly_input.yml')   # In this file the user specifies which input they want
    assembly = assemble_complexes(assembly_input_path=USER_INPUT)

    #%% ==============    Doublecheck refactoring    ==================
    from test.Integration_Test import IntegrationTest
    test = IntegrationTest(new_dir=assembly.output_path, old_dir=Path(assembly.output_path.parent, 'output_benchmark'))
    test.compare_all()

    print('Done!')
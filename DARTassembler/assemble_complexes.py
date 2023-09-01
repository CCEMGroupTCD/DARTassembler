"""
This module is a wrapper for the assembly code. It takes in an assembly input file and returns an assembly object.
"""
from DARTassembler.src.assembly.DART_Assembly import DARTAssembly
from pathlib import Path
from typing import Union

def assemble_complexes(assembly_input_path: Union[str, Path]):
    assembly = DARTAssembly(assembly_input_path=assembly_input_path)
    assembly.run_all_batches()

    return assembly



# Integration test, to check if everything is working and the same as before.
if __name__ == "__main__":
    from DARTassembler.src.constants.Paths import project_path

    USER_INPUT = project_path().extend('src05_Assembly_Refactor', 'assembly_input.yml')   # In this file the user specifies which input they want
    assembly = assemble_complexes(assembly_input_path=USER_INPUT)

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    test = IntegrationTest(new_dir=assembly.output_path, old_dir=Path(assembly.output_path.parent, 'output_benchmark'))
    test.compare_all()

    print('Done!')
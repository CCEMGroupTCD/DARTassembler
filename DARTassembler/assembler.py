"""
This module is a wrapper for the assembly code. It takes in an assembly input file and returns an assembly object.
"""
from DARTassembler.src.assembly.DART_Assembly import DARTAssembly
from pathlib import Path
from typing import Union


def assembler(assembly_input_path: Union[str, Path], delete_output_dir: bool = False) -> DARTAssembly:
    assembly = DARTAssembly(assembly_input_path=assembly_input_path, delete_output_dir=delete_output_dir)
    assembly.run_all_batches()

    return assembly


# Integration test, to check if everything is working and the same as before.
if __name__ == "__main__":
    from DARTassembler.src.constants.Paths import project_path

    USER_INPUT = project_path().extend('dev/oerdata/testbatch_data/input_files/oer_test_assembler.yml')  # In this file the user specifies which input they want
    assembly = assembler(assembly_input_path=USER_INPUT)

    # %% ==============    Doublecheck refactoring    ==================
    # from dev.test.Integration_Test import IntegrationTest
    #
    # test = IntegrationTest(new_dir=project_path().extend('examples', 'Pd_Ni_Cross_Coupling', 'output', 'data_before_g16_calcs', 'DART_Example_Pd_Ni_Complexes'),
    #                        old_dir=project_path().extend('Legacy_unit_test_do_not_delete', 'DART_Example_Pd_Ni_Complexes'))
    # test.compare_all()

    print('Done!')

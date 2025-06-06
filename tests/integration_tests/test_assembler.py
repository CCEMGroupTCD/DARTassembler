"""
Integration test for the assembly of complexes.
"""
from pathlib import Path
from DARTassembler.src.constants.paths import project_path
from DARTassembler.src.assembler.assembler import Assembler


def test_assembler():
    assembly_input = project_path().extend('tests', 'integration_tests', 'assembler', 'data_input', 'test_assembly_input.yml')
    assembly = Assembler.run_from_yaml(assembly_input)

    #%% ==============    Doublecheck refactoring    ==================
    from DARTassembler.src.misc.tests import IntegrationTest
    old_dir = Path(assembly.output_directory.parent, 'benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=assembly.output_directory, old_dir=old_dir)
        test.compare_all()
        print('Test for assembly of complexes passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return assembly


if __name__ == "__main__":
    assembly = test_assembler()
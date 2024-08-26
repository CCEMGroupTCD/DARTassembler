"""
Integration test for the assembly of complexes.
"""
from DARTassembler.assembler import assembler
from pathlib import Path
from DARTassembler.src.constants.Paths import project_path

def test_assemble_complexes():
    assembly_input = project_path().extend('testing', 'integration_tests', 'assembly', 'data_input', 'test_assembly_input.yml')   # In this file the user specifies which input they want
    assembly = assembler(assembly_input_path=assembly_input)

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    old_dir = Path(assembly.output_path.parent, 'benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=assembly.output_path, old_dir=old_dir)
        test.compare_all()
        print('Test for assembly of complexes passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return assembly


if __name__ == "__main__":
    assembly = test_assemble_complexes()
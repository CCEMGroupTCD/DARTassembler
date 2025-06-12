"""
Integration test for the assembly of complexes.
"""
from pathlib import Path
from DARTassembler.src.constants.Paths import project_path

def test_assembler(method):
    if method == "A":
        from dev.DART_refactoring_to_v1_1_0.assembler import Assembler
        assembly_input = project_path().extend('tests', 'integration_tests', 'assembler', 'data_input', 'test_assembly_input.yml') #

    elif method == "B":
        from dev.DART_refactoring_to_v1_2_0.assembler import Assembler
        assembly_input = project_path().extend('dev', 'DART_refactoring_to_v1_2_0', 'assembler.yml')
    else:
        raise ValueError(f"Unknown method '{method}'! Use 'A' or 'B'.")



    assembly = Assembler.run_from_yaml(assembly_input)

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    old_dir = Path(assembly.output_directory.parent, 'benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=assembly.output_directory, old_dir=old_dir)
        test.compare_all()
        print('Test for assembly of complexes passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return assembly


if __name__ == "__main__":
    assembly = test_assembler(method="B")
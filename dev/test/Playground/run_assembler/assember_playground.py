from DARTassembler.assemble_complexes import assemble_complexes
from pathlib import Path

if __name__ == '__main__':

    assembly_input = 'input_assembler_playground.yml'
    assembly = assemble_complexes(assembly_input_path=assembly_input)

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    old_dir = Path(assembly.output_path.parent, 'benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=assembly.output_path, old_dir=old_dir)
        test.compare_all()
        print('Test for assembly of complexes passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')
from src05_Assembly_Refactor.DART_Assembly import DARTAssembly
from pathlib import Path


def main(assembly_input_path: str):
    assembly = DARTAssembly(assembly_input_path=assembly_input_path)
    assembly.run_all_batches()

    return assembly



if __name__ == "__main__":

    USER_INPUT = 'EXAMPLE_PN_Pd_Ni.yml'   # In this file the user specifies which input they want
    assembly = main(assembly_input_path=USER_INPUT)



    #%% ==============    Doublecheck refactoring    ==================
    from test.Integration_Test import IntegrationTest
    test = IntegrationTest(new_dir=assembly.output_path, old_dir=Path(assembly.output_path.parent, 'DART_Example_Pd_Ni_Complexes'))
    test.compare_all()


    print('Done!')
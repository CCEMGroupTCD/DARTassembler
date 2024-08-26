"""
Integration test for the filtering of ligands.
"""
from DARTassembler.ligandfilters import ligandfilters
from DARTassembler.src.constants.Paths import project_path

def test_filter_ligands(nmax=3000):
    filter_ligands_path = project_path().extend('testing', 'integration_tests', 'ligandfilters', 'data_input', 'ligandfilters.yml')   # In this file the user specifies which input they want
    filters = ligandfilters(filter_input_path=filter_ligands_path, nmax=nmax)

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    old_dir = project_path().extend('testing', 'integration_tests', 'ligandfilters', 'benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=filters.output_ligand_db_path.parent, old_dir=old_dir)
        test.compare_all()
        print('Test for ligand filters passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return filters


if __name__ == "__main__":
    nmax = 3000
    filters = test_filter_ligands(nmax=nmax)
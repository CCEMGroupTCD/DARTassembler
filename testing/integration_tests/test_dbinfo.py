"""
Integration test for outputting a csv file of the ligand database.
"""
from DARTassembler.make_ligand_db_csv import make_ligand_db_csv
from DARTassembler.src.constants.Paths import project_path
from pathlib import Path


def test_make_ligand_db_csv(nmax=1000):
    output_path = project_path().extend('testing', 'integration_tests', 'make_ligand_db_csv', 'data_output', 'MetaLig_v1.7.csv')
    db = make_ligand_db_csv(input_path='metalig', output_path=output_path, nmax=nmax)

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    old_dir = Path(str(output_path).replace('/data_output/', '/benchmark_data_output/'))
    if old_dir.exists():
        test = IntegrationTest(new_dir=output_path.parent, old_dir=old_dir.parent)
        test.compare_all()
        print('Test for assembly of complexes passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return db


if __name__ == "__main__":
    nmax = 1000

    db = test_make_ligand_db_csv(nmax=nmax)


"""
This file is for testing the DART package.
"""
from time import sleep
from testing.integration_tests.test_assemble_complexes import test_assemble_complexes
from testing.integration_tests.test_ligand_extraction import test_ligand_extraction
from testing.integration_tests.test_filter_ligands import test_filter_ligands


if __name__ == '__main__':
    sleeptime = 5

    c, ulig, lig, df_unique_ligands, df_complexes = test_ligand_extraction()
    sleep(sleeptime)

    filters = test_filter_ligands()
    sleep(sleeptime)

    assembly = test_assemble_complexes()
    sleep(sleeptime)

    print('All tests passed!')

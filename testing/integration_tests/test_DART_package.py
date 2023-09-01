"""
This file is for testing the DART package.
"""
from time import sleep
from testing.integration_tests.test_assemble_complexes import test_assemble_complexes
from testing.integration_tests.test_ligand_extraction import test_ligand_extraction
from testing.integration_tests.test_filter_ligands import test_filter_ligands


if __name__ == '__main__':
    sleeptime = 3

    test_ligand_extraction()
    sleep(sleeptime)

    test_filter_ligands()
    sleep(sleeptime)

    test_assemble_complexes()
    sleep(sleeptime)

    print('All tests passed!')

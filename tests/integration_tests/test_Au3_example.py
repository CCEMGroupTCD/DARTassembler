"""
Integration test for the filtering and assembly of Ni(II) square-planar complexes, the example from the DART paper.
"""
import itertools
import sys
print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/Users/timosommer/PhD/projects/DARTassembler'])
import os
from pathlib import Path
from DARTassembler.src.constants.paths import project_path
from DARTassembler.src.assembler.assembler import Assembler
from DARTassembler import LigandFilters
from shutil import rmtree

def test_Au3():
    outdir = project_path().extend('tests', 'integration_tests', 'Au3_example', 'data_output', 'assembly')
    assembly_input = outdir.parent.parent / 'data_input' / 'Au3_assembly_input.yaml'
    # Remove the output directory if it exists to start fresh
    if outdir.exists():
        rmtree(outdir)
    old_cwd = Path.cwd()    # Save the current working directory to return to it later
    try:
        outdir.mkdir(parents=True, exist_ok=True)
        os.chdir(outdir)
        assembly = Assembler.run_from_yaml(assembly_input)
    finally:
        # Change back to the original working directory
        os.chdir(old_cwd)

    #%% ==============    Doublecheck refactoring    ==================
    from DARTassembler.src.misc.tests import IntegrationTest
    old_dir = Path(outdir.parent.parent, 'benchmark_data_output').resolve()
    if old_dir.exists():
        test = IntegrationTest(new_dir=outdir.parent, old_dir=old_dir)
        test.compare_all()
        print('Test for assembly of complexes passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return assembly


if __name__ == "__main__":
    n_max = None
    outdir = project_path().extend('tests', 'integration_tests', 'Au3_example', 'data_output')

    Cl_filters = [
        {'filter': 'composition', 'elements': 'Cl', 'instruction': 'must_contain_and_only_contain', 'only_donors': False},
        {'filter': 'property', 'name': 'archetype', 'values': ['1-mono']}
    ]
    CN = [
        {'filter': 'property', 'name': 'n_denticities', 'values': [2]},
        {'filter': 'property', 'name': 'n_haptic_groups', 'values': [0]},
        {'filter': 'property', 'name': 'charge', 'values': [-1]},
        {'filter': 'property', 'name': 'archetype', 'values': ['2-cis']},    # This filter was not in the original study, we just manually removed the 2-trans ligand afterwards
        {'filter': 'composition', 'elements': 'CN', 'instruction': 'must_contain_and_only_contain', 'only_donors': True},
    ]


    # Remove the output directory if it exists to start fresh
    if outdir.exists():
        rmtree(outdir)
    old_cwd = Path.cwd()    # Save the current working directory to return to it later
    try:
        outdir.mkdir(parents=True, exist_ok=True)
        os.chdir(outdir)
        filter = LigandFilters(db='metalig', n=n_max)
        Cl_ligandnames = filter.run(Cl_filters, outpath=Path('ligands')/'Cl.jsonlines', metal=True, dbinfo=False)
        CN_ligandnames = filter.run(CN, outpath=Path('ligands')/'CN.jsonlines', metal=True, dbinfo=False)
    finally:
        # Change back to the original working directory
        os.chdir(old_cwd)

    test_Au3()


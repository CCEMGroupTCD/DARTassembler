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

def test_sqpl_Ni():
    outdir = project_path().extend('tests', 'integration_tests', 'assembler', 'data_output')
    assembly_input = outdir.parent / 'data_input' / 'test_assembly_input.yaml'
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
    old_dir = Path(assembly.output_directory.parent, 'benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=assembly.output_directory, old_dir=old_dir)
        test.compare_all()
        print('Test for assembly of complexes passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return assembly


if __name__ == "__main__":
    n_max = 1000
    outdir = project_path().extend('tests', 'other_tests', 'sqpl_Ni', 'data_output')

    Br_filters = [
        {'filter': 'composition', 'elements': 'Br', 'instruction': 'must_contain_and_only_contain', 'only_donors': False},
        {'filter': 'property', 'name': 'archetype', 'values': ['1-mono']}
    ]
    NN_filters = [
        {'filter': 'composition', 'elements': 'N2', 'instruction': 'must_contain_and_only_contain', 'only_donors': True},
        {'filter': 'composition', 'elements': 'CHNF', 'instruction': 'must_only_contain_in_any_amount', 'only_donors': False},
        {'filter': 'property', 'name': 'archetype', 'values': ['2-cis']},
    ]
    Hap_filters = [
        {'filter': 'composition', 'elements': 'C', 'instruction': 'must_only_contain_in_any_amount', 'only_donors': True},
        {'filter': 'property', 'name': 'n_haptic_groups', 'values': [1]},
        {'filter': 'property', 'name': 'archetype', 'values': ['1-mono']},
    ]




    # Remove the output directory if it exists to start fresh
    if outdir.exists():
        rmtree(outdir)
    old_cwd = Path.cwd()    # Save the current working directory to return to it later
    try:
        outdir.mkdir(parents=True, exist_ok=True)
        os.chdir(outdir)
        filter = LigandFilters(db='metalig', n=n_max)
        Br_ligandnames = filter.run(Br_filters, outpath='Br.jsonlines', metal=False, dbinfo=False)
        NN_ligandnames = filter.run(NN_filters, outpath='NN.jsonlines', metal=False, dbinfo=False)
        Hap_ligandnames = filter.run(Hap_filters, outpath='Hap.jsonlines', metal=False, dbinfo=False)
    finally:
        # Change back to the original working directory
        os.chdir(old_cwd)

    # For each set, get a histogram of charges for each ligand
    print('NN ligand charges:')
    import pandas as pd
    NN_charges = pd.Series([filter.db.db[name].charge for name in NN_ligandnames]).value_counts().sort_index()
    print(NN_charges)
    print('Br ligand charges:')
    Br_charges = pd.Series([filter.db.db[name].charge for name in Br_ligandnames]).value_counts().sort_index()
    print(Br_charges)
    print('Hap ligand charges:')
    Hap_charges = pd.Series([filter.db.db[name].charge for name in Hap_ligandnames]).value_counts().sort_index()
    print(Hap_charges)
    # Get all combinations of charges that give overall -1 for NN and Hap
    import itertools
    combinations = list(itertools.product(NN_charges.index, Hap_charges.index, Br_charges.index))
    total_n = 0
    for NN_charge, Hap_charge, Br_charge in combinations:
        if Hap_charge + NN_charge + Br_charge == -1:
            n = Hap_charges[Hap_charge] * NN_charges[NN_charge]
            print(f'Hap: {Hap_charge}, NN: {NN_charge}, Br: {Br_charge} n: {n}')
            total_n += n
    print('Total combinations with correct charge (NN + Hap = -1):', total_n)
    print('Done!')


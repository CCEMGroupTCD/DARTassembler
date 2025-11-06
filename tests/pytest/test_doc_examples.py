"""
Integration test for the assembly of complexes.
"""
import sys; print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/Users/timosommer/PhD/projects/DARTassembler'])
import os
from pathlib import Path
from DARTassembler.src.constants.paths import project_path
from DARTassembler import Assembler, LigandFilters
from shutil import rmtree

def test_doc_examples():
    outdir = project_path().extend('tests', 'pytest', 'doc_examples', 'data_output')

    C6H6_filters_path = outdir.parent / 'data_input' / 'C6H6_ligandfilters.yml'
    HMDS_filters_path = outdir.parent / 'data_input' / 'amide_ligandfilters.yml'
    assembly_input = outdir.parent / 'data_input' / 'test_doc_examples.yml'
    # Remove the output directory if it exists to start fresh
    if outdir.exists():
        rmtree(outdir)
    old_cwd = Path.cwd()    # Save the current working directory to return to it later
    try:
        outdir.mkdir(parents=True, exist_ok=True)
        os.chdir(outdir)
        LigandFilters.run_from_yaml(C6H6_filters_path)
        LigandFilters.run_from_yaml(HMDS_filters_path)
        assembly = Assembler.run_from_yaml(assembly_input)
    finally:
        # Change back to the original working directory
        os.chdir(old_cwd)

    #%% ==============    Doublecheck refactoring    ==================
    from DARTassembler.src.misc.tests import IntegrationTest
    old_dir = Path(outdir.parent, 'benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=outdir, old_dir=old_dir)
        test.compare_all()
        print('Test for assembly of complexes passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return assembly


if __name__ == "__main__":
    assembly = test_doc_examples()
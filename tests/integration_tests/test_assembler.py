"""
Integration test for the assembly of complexes.
"""
import sys; print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/Users/timosommer/PhD/projects/DARTassembler'])
import os
from pathlib import Path
from DARTassembler.src.constants.paths import project_path
from DARTassembler.src.assembler.assembler import Assembler
from shutil import rmtree

def test_assembler():
    outdir = project_path().extend('tests', 'integration_tests', 'assembler', 'data_output')
    assembly_input = outdir.parent / 'data_input' / 'test_assembly_input.yml'
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
    assembly = test_assembler()
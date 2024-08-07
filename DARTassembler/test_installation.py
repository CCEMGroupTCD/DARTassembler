"""
This file tests the installation of the DARTassembler package. It executes all the modules and checks if any errors occur. It does not check if the output is correct, only if the code runs without errors.
"""
import os
from pathlib import Path
import tempfile
from typing import Union
from DARTassembler import filter_ligands, assemble_complexes, save_dbinfo
from DARTassembler.src.constants.Paths import test_installation_dirpath


def test_module_in_directory(function, input: Union[Path, str], outdir: [Path, str, None]) -> None:
    """
    Tests a DART module without leaving any files saved. This is done by copying the test files to a temporary directory and running the module there.
    :param function: DART module to test, given as a python function
    :param input: Input of the module
    :param outdir: Output directory of the module. If None, a temporary directory is used and no files are saved.
    """
    prev_cwd = os.getcwd()  # Save the current working directory to return to it after the test

    if outdir is None:  # Leave no files saved
        print(f'==========  Testing module {function.__name__} in temporary directory  ==========')
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            function(input)
    else:               # Save the files in the specified directory
        print(f'==========  Testing module {function.__name__} in directory {Path(outdir).name}  ==========')
        os.chdir(outdir)
        function(input)

    os.chdir(prev_cwd)  # Reset working directory

def test_assembler_installation(outdir: Union[Path, str, None]) -> None:
    """
    Test the assembler module in a temporary directory without leaving any files saved.
    :return: None
    """
    input_filepath = Path(test_installation_dirpath, 'assembler.yml').resolve()
    test_module_in_directory(function=assemble_complexes, input=input_filepath, outdir=outdir)

def test_ligandfilters_installation(outdir: Union[Path, str, None]) -> None:
    """
    Test the ligandfilters module in a temporary directory without leaving any files saved.
    :return: None
    """
    input_filepath = Path(test_installation_dirpath, 'ligandfilters.yml').resolve()
    test_module_in_directory(function=filter_ligands, input=input_filepath, outdir=outdir)

def test_dbinfo_installation(outdir: Union[Path, str, None]) -> None:
    """
    Test the save_dbinfo module in a temporary directory without leaving any files saved.
    :return: None
    """
    test_module_in_directory(function=save_dbinfo, input='test_metalig', outdir=outdir)



def run_installation_test(outdir: Union[Path, str, None] = None) -> None:

    test_assembler_installation(outdir=outdir)

    test_ligandfilters_installation(outdir=outdir)

    test_dbinfo_installation(outdir=outdir)

    print('Done! All modules tested successfully.')

if __name__ == '__main__':
    run_installation_test()




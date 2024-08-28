"""
This file tests the installation of the DARTassembler package. It executes all the modules and checks if any errors occur. It does not check if the output is correct, only if the code runs without errors.
"""
import os
import io
import sys
from pathlib import Path
import tempfile
from typing import Union
from DARTassembler import ligandfilters, assembler, dbinfo
from DARTassembler.src.constants.Paths import test_installation_dirpath
from contextlib import redirect_stdout


def test_module_in_directory(function, input: Union[Path, str], outdir: [Path, str, None], suppress_stdout: bool = True) -> None:
    """
    Tests a DART module without leaving any files saved. This is done by copying the test files to a temporary directory and running the module there.
    :param function: DART module to test, given as a python function
    :param input: Input of the module
    :param outdir: Output directory of the module. If None, a temporary directory is used and no files are saved.
    :param suppress_stdout: If True, the stdout of the module is suppressed. If False, the stdout is printed.
    """
    prev_cwd = os.getcwd()  # Save the current working directory to return to it after the test

    # Suppress stdout if suppress_stdout is True
    stdout = io.StringIO() if suppress_stdout else sys.stdout
    with redirect_stdout(stdout):

        if outdir is None:  # Leave no files saved
            testing_title = f'  Testing module {function.__name__} in temporary directory  '
            print(f'{testing_title:#^80}')
            with tempfile.TemporaryDirectory() as tmpdir:
                os.chdir(tmpdir)
                function(input)
        else:               # Save the files in the specified directory
            testing_title = f'  Testing module {function.__name__} in directory {Path(outdir).name}  '
            print(f'{testing_title:#^80}')
            outdir = Path(outdir)
            outdir.mkdir(parents=True, exist_ok=True)
            os.chdir(outdir)
            function(input)

    os.chdir(prev_cwd)  # Reset working directory

def test_assembler_installation(outdir: Union[Path, str, None], suppress_stdout: bool = True) -> None:
    """
    Test the assembler module in a temporary directory without leaving any files saved.
    :return: None
    """
    input_filepath = Path(test_installation_dirpath, 'assembler.yml').resolve()
    test_module_in_directory(function=assembler, input=input_filepath, outdir=outdir, suppress_stdout=suppress_stdout)

def test_ligandfilters_installation(outdir: Union[Path, str, None], suppress_stdout: bool = True) -> None:
    """
    Test the ligandfilters module in a temporary directory without leaving any files saved.
    :return: None
    """
    input_filepath = Path(test_installation_dirpath, 'ligandfilters.yml').resolve()
    test_module_in_directory(function=ligandfilters, input=input_filepath, outdir=outdir, suppress_stdout=suppress_stdout)

def test_dbinfo_installation(outdir: Union[Path, str, None], suppress_stdout: bool = True) -> None:
    """
    Test the save_dbinfo module in a temporary directory without leaving any files saved.
    :return: None
    """
    test_module_in_directory(function=dbinfo, input='test_metalig', outdir=outdir, suppress_stdout=suppress_stdout)



def installtest(outdir: Union[Path, str, None] = None, suppress_stdout: bool = True) -> None:

    print(f"{'  Testing DART installation...  ':=^80}")

    print('Testing DART Assembler Module...')
    test_assembler_installation(outdir=outdir, suppress_stdout=suppress_stdout)
    print('Successful.')

    print('Testing DART Ligand Filters Module...')
    test_ligandfilters_installation(outdir=outdir, suppress_stdout=suppress_stdout)
    print('Successful.')

    print('Testing DART DBinfo Module...')
    test_dbinfo_installation(outdir=outdir, suppress_stdout=suppress_stdout)
    print('Successful.')

    print('Done! All modules tested successfully. Exiting DART InstallTest Module.')

if __name__ == '__main__':
    installtest()




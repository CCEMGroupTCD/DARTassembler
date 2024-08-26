"""
This module returns the default configuration files the ligand filters and the assembler. These files can be used as templates for the user to adapt to their specific needs.
"""
from pathlib import Path
from typing import Union
from DARTassembler.src.constants.Paths import default_assembler_yml_path, default_ligandfilters_yml_path
from shutil import copyfile

def configs(outdir: str) -> None:
    """
    Get the default configuration files for the ligand filters and the assembler and save them to the specified output path.
    :return: None
    """
    print(f"Getting default configuration files...")

    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    # Copy assembler.yml
    filename = default_assembler_yml_path.name
    print(f'\t- {filename}')
    dest = Path(outpath, filename)
    copyfile(default_assembler_yml_path, dest)

    # Copy ligandfilters.yml
    filename = default_ligandfilters_yml_path.name
    print(f'\t- {filename}')
    dest = Path(outpath, filename)
    copyfile(default_ligandfilters_yml_path, dest)

    print(f"Saved in current directory.")
    print('Done!')

    return





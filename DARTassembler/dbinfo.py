"""
This module reads in a ligand db from file and saves a .csv file with an overview of the ligands.
"""
from typing import Union
from pathlib import Path

from DARTassembler.src.metalig.metalig_utils import get_correct_ligand_db_path_from_input
from DARTassembler.src.constants.Paths import default_ligand_db_path

import DARTassembler.src.constants.Paths
from DARTassembler.src.ligand_extraction.DataBase import LigandDB

def get_ligand_csv_output_path(output_path: Union[str, Path], input_path: Union[str, Path]):
    """
    Save to csv. If no output path is specified, save to the same directory as the input file with the same name but with the .csv extension.
    """
    if output_path is None:
        if str(input_path).endswith('.csv'):
            raise ValueError(f"Input path {input_path} ends with .csv. Please specify an output path explicitly.")
        else:
            current_dir = Path.cwd()
            output_path = current_dir.joinpath(Path(input_path).with_suffix('.csv').name)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    return output_path

def dbinfo(input_path: Union[str, Path], output_path: Union[str, Path, None] = None, nmax: Union[int, None] = None):
    """
    Reads in the given ligand database and saves a .csv file and a concatenated .xyz file with an overview of the ligands.
    :param input_path: Path to the ligand database
    :param output_path: Path to the output .csv file. If None, the output file will be saved in the same directory as the input file with the same name as the input file but with the .csv extension.
    :param nmax: Maximum number of ligands to be read in from the initial full ligand database. If None, all ligands are read in. This is useful for testing purposes.
    :return: LigandDB object
    """
    input_path = get_correct_ligand_db_path_from_input(input_path)
    if input_path is None:
        raise ValueError(f"Invalid ligand database path.")

    print(f"Starting DART DBinfo Module.")
    print(f'Input ligand database: {input_path.name}')
    db = LigandDB.load_from_json(input_path, n_max=nmax)

    print('Saving ligand info and structures...')

    # Save to csv
    output_path = get_ligand_csv_output_path(output_path, input_path)
    db.save_reduced_csv(output_path)
    print(f'  - Saved .csv to {output_path.name}')

    # Save to concatenated xyz file
    xyz_filename = str(Path(f'concat_{output_path.with_suffix("").name}.xyz'))
    xyz_output_path = output_path.parent.joinpath(xyz_filename)
    db.save_concat_xyz(xyz_output_path, with_metal=True)
    print(f'  - Saved .xyz to {xyz_filename}')

    print(f"Done! All info files saved. Exiting DART DBinfo Module.")

    return db
"""
Utility functions for input and output.
"""
import bz2
import json
import sys
import yaml

from DARTassembler.src.metalig.metalig_utils import get_correct_ligand_db_path_from_input
from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand, RCA_Complex
import numpy as np
from datetime import datetime, date, timedelta
from DARTassembler.src.ligand_extraction.utilities import get_duration_string
from typing import Union
from pathlib import Path, PurePath
import jsonlines
from tqdm import tqdm
import functools
import ase
import zipfile
import shutil
import os
from DARTassembler.src.constants.Paths import default_ligand_db_path, test_ligand_db_path

def check_if_MetaLig_exists_else_uncompress_from_zip(delete_zip=False):
    """
    Checks if the MetaLig database exists as uncompressed file, and if not uncompresses it.
    """
    files = [default_ligand_db_path, test_ligand_db_path]

    for unzipped_file in files:
        zip_file = Path(str(unzipped_file) + '.bz2')

        if not unzipped_file.exists():
            name = Path(zip_file).name
            if not Path(zip_file).exists():
                raise FileNotFoundError(f"DART Error: Could not find MetaLig database file {name} at {Path(zip_file).resolve()}.")

            db_dir = Path(zip_file).parent
            try:
                print(f"Uncompressing MetaLig database file {name}...")
                uncompress_file(zip_file, db_dir)
            except Exception as e:
                raise Exception(f"DART Error: Could not uncompress MetaLig database file {name} at {Path(zip_file).resolve()}. Error message: {e}")

            if delete_zip:
                Path(zip_file).unlink()

    return

def compress_file(file_path, output_zip_path=None, compression_level=6):
    """
    Compress a file into a zip file.
    """
    if output_zip_path is None:
        output_zip_path = str(file_path) + '.zip'

        with zipfile.ZipFile(output_zip_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=compression_level) as zipf:
            zipf.write(file_path, os.path.basename(file_path))

    return

def uncompress_file(zip_file_path, output_dir = None):
    """
    Uncompress a zip file into a file.
    """
    zip_file_path = Path(zip_file_path)

    if output_dir is None:
        output_dir = Path(zip_file_path).parent

    # Create a temporary directory to extract the files
    temp_dir = Path(output_dir, 'temp_extract')
    temp_dir.mkdir(parents=True, exist_ok=True)

    # Extract all files to the temporary directory
    if zip_file_path.suffix == '.zip':
        with zipfile.ZipFile(zip_file_path, 'r') as zipf:
            zipf.extractall(temp_dir)

    elif zip_file_path.suffix == '.bz2':
        newfile = zip_file_path.with_suffix('')
        with open(newfile, 'wb') as new_file, bz2.BZ2File(zip_file_path, 'rb') as file:
            for data in iter(lambda: file.read(100 * 1024), b''):
                new_file.write(data)

    # Move the contents of the temporary directory to the output directory
    for filename in os.listdir(temp_dir):
        shutil.move(Path(temp_dir, filename), Path(output_dir, filename))

    # Remove the temporary directory
    Path(temp_dir).rmdir()

    return


class NumpyEncoder(json.JSONEncoder):
    """Special json encoder. This is important to use in json.dump so that if json encounters a np.array, it converts it to a list automatically, otherwise errors arise. Use like this:
    dumped = json.dump(dic, cls=NumpyEncoder)
    """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, np.str_):
            return str(obj)
        elif isinstance(obj, np.string_):
            return str(obj)
        elif isinstance(obj, Path):
            return str(obj)
        elif isinstance(obj, (datetime, date)):
            return obj.isoformat()
        elif isinstance(obj, timedelta):
            return str(obj)
        return json.JSONEncoder.default(self, obj)




def check_if_return_entry(i: int, n_max: Union[int, list]=None) -> bool:
    """
    Check if the entry should be returned. It will not be returned if the index i is larger than n_max.
    """
    # Accept False, None or np.inf to disable n_max
    if n_max is None or n_max == False or n_max is np.inf:
        return True

    # If n_max is an integer, check if the index is smaller than n_max
    is_good = i < n_max

    return is_good

def ensure_path_exists(path: Union[str, Path]) -> Path:
    """
    Check if the path exists. If not, raise an error.
    """
    try:
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f'Could not find file {path}')
    except TypeError:
        raise TypeError(f'Expected path to be a string or Path, but got {type(path)} for file {path}')

    return path

def load_json(path: Union[str, Path], n_max: int=None, show_progress: bool=True) -> dict:
    """
    Load a JSON or JSON Lines file. If the file is a JSON Lines file, it is converted to a dictionary.
    :param path: Path to the JSON or JSON Lines file
    :return: Dictionary with the contents of the file
    """
    db = {key: value for key, value in iterate_over_json(path, n_max=n_max, show_progress=show_progress)}

    return db

def iterate_over_json(path: Union[str, Path], n_max: int=None, show_progress: bool=True) -> tuple[str, dict]:
    """
    Iterate over a JSON or JSON Lines file and yield the key and value of each entry.
    :param path: Path to the JSON or JSON Lines file
    :return: Tuple with the key and value of each entry
    """
    path = ensure_path_exists(path)
    try:
        # Try to load as normal JSON file first
        with open(path, 'r') as file:
            db = json.load(file)

            if set(db.keys()) == {'key', 'value'}:
                # If the file is a JSON Lines file, go into correct mode
                raise UserWarning

            for i, (key, value) in enumerate(db.items()):
                if check_if_return_entry(i, n_max):
                    yield key, value
                else:
                    return

    except (json.JSONDecodeError, UserWarning):
        # If normal JSON fails, try to load as JSON Lines
        with jsonlines.open(path, 'r') as reader:
            for i, line in tqdm(enumerate(reader), disable=not show_progress, desc='Load json'):
                key, value = line['key'], line['value']
                if check_if_return_entry(i, n_max):
                    yield key, value
                else:
                    return

    return

def iterate_over_jsonlines(path: Union[str, Path], n_max: int=None, show_progress: bool=True) -> tuple[str, dict]:
    """
    Iterate over a JSON Lines file and yield the key and value of each entry.
    :param path: Path to the JSON Lines file
    :return: Tuple with the key and value of each entry
    """
    path = ensure_path_exists(path)
    with jsonlines.open(path, 'r') as reader:
        for i, line in tqdm(enumerate(reader), disable=not show_progress, desc='Load jsonlines'):
            key, value = line['key'], line['value']
            if check_if_return_entry(i, n_max):
                yield key, value
            else:
                return

    return

def get_n_entries_of_json_db(path: Union[str, Path]) -> int:
    """
    Get the number of entries in a JSON or JSON Lines file.
    :param path: Path to the JSON or JSON Lines file
    :return: Number of entries in the file
    """
    n_entries = 0
    for _ in iterate_over_json(path):
        n_entries += 1

    return n_entries

def save_json(db: dict, path: Union[str, Path], **kwargs):
    with open(path, 'w') as file:
        json.dump(db, file, cls=NumpyEncoder, **kwargs)

    return

def save_jsonlines(db: dict, path: Union[str, Path]):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with jsonlines.open(path, mode='w', dumps=functools.partial(json.dumps, cls=NumpyEncoder)) as writer:
        for key, value in db.items():
            data = {'key': key, 'value': value}
            writer.write(data)

    return

def load_jsonlines(path: Union[str, Path], n_max: int=None, show_progress: bool=True) -> dict:
    db = {key: value for key, value in iterate_over_jsonlines(path, n_max=n_max, show_progress=show_progress)}
    return db

def check_molecule_value(output: str):
    possible_values = ['dict', 'class']
    if not output in possible_values:
        raise ValueError(f'Unknown value for `output`: {output}')

    return

def iterate_complex_db(path: Union[str, Path], molecule: str='dict', n_max=None, show_progress: bool=True) -> dict:
    check_molecule_value(molecule)  # Check if the molecule value is valid
    for name, mol in tqdm(iterate_over_json(path, n_max=n_max, show_progress=False), disable=not show_progress, desc='Load complex db'):
        if molecule == 'class':
            mol = RCA_Complex.read_from_mol_dict(mol)
        yield name, mol


def load_complex_db(path: Union[str, Path], molecule: str='dict', n_max=None, show_progress: bool=True) -> dict:
    db = {name: mol for name, mol in iterate_complex_db(path=path, molecule=molecule, n_max=n_max, show_progress=show_progress)}
    return db

def load_full_ligand_db(path: Union[str, Path], molecule: str='dict') -> dict:
    start = datetime.now()

    check_molecule_value(molecule)

    if 'complex' in path.stem:
        # If the complex database is provided, the ligand database is a subset of it
        db = {}
        for c_id, c in iterate_over_json(path):
            for lig in c['ligands']:
                db[lig['name']] = lig
    else:
        db = load_json(path)

    if molecule == 'class':
        db = {name: RCA_Ligand.read_from_mol_dict(mol) for name, mol in db.items()}

    duration = get_duration_string(start=start)
    print(f'Loaded full ligand db. Time: {duration}. ')
    return db

def iterate_unique_ligand_db(path: Union[str, Path], molecule: str= 'dict', n_max=None, show_progress: bool=False) -> dict:
    check_molecule_value(molecule)  # Check if the molecule value is valid
    db_path = get_correct_ligand_db_path_from_input(path)
    if db_path is None:
        raise ValueError(f'Invalid ligand database path specified: {path}')

    filename = Path(db_path).name
    for name, mol_dict in tqdm(iterate_over_json(db_path, n_max=n_max, show_progress=False), disable=not show_progress, desc=f'Load ligand db `{filename}`', file=sys.stdout, unit=' ligands'):
        if molecule == 'class':
            try:
                mol = RCA_Ligand.read_from_mol_dict(mol_dict)
            except Exception as e:
                raise ValueError(f"Error: the provided file '{db_path}' seems to not be a valid ligand database file. Internal error message:\n{e}.")
        else:
            mol = mol_dict
        yield name, mol

def load_unique_ligand_db(path: Union[str, Path], molecule: str='dict', n_max=None, show_progress: bool=True) -> dict:
    db = {name: mol for name, mol in iterate_unique_ligand_db(path=path, molecule=molecule, n_max=n_max, show_progress=show_progress)}
    return db

def save_complex_db(db: dict, path: Union[str, Path]):
    start = datetime.now()

    save_json(db, path)

    duration = get_duration_string(start=start)
    print(f"Complex database saved to {path}. Time: {duration}.")

    return

def save_full_ligand_db(db: dict, path: Union[str, Path]):
    start = datetime.now()

    save_json(db, path)

    duration = get_duration_string(start=start)
    print(f"Full ligand database saved to {path}. Time: {duration}.")

    return

def save_unique_ligand_db(db: dict, path: Union[str, Path]):
    start = datetime.now()

    save_json(db, path)

    duration = get_duration_string(start=start)
    print(f"Unique ligand database saved to {path}. Time: {duration}.")

    return

def write_yaml(path: Union[str, Path], data: dict) -> None:
    with open(path, 'w') as file:
        yaml.dump(data, file)

    return

def read_yaml(path: Union[str, Path]) -> dict:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'Filepath does not exist: {path}.')

    with open(path, 'r') as f:
        txt = f.read()

        try:
            data = yaml.load(txt, yaml.SafeLoader)

        except yaml.YAMLError as exc:
            # Print a nice error message pointing to the line where the error is
            error_start = f"\n\n-->Error while parsing the input YAML file '{path}':"
            error_end = 'The error probably is either on this line or the line before. A common error is wrong indentation. Please correct input file and retry.'
            if hasattr(exc, 'problem_mark'):
                if exc.context != None:
                    raise Exception(f'{error_start}\n\tThe issue seems to be here:\n\t{exc.problem_mark}\n\tYAML error message: {exc.problem}{exc.context}\n\t{error_end}')
                else:
                    raise Exception(f'{error_start}\n\tThe issue seems to be here:\n\t{exc.problem_mark}\n\tYAML error message: {exc.problem}\n\t{error_end}')

            else:
                raise Exception(f'There was an error while reading the YAML file {path}. Please make sure the file is a proper yaml file. For example, a common error is that the indentation might be wrong. This is the error message from yaml, please google it if you don\'t know how to fix it: {exc}')

    return data

def safe_read_yaml(data: Union[str, Path, list, dict]) -> Union[dict,list]:
    """
    Safely read a YAML file or return the data if it is already a dict/list.
    """
    try:
        data = read_yaml(data)
    except TypeError:
        # If data is already a dict/list, just return it
        if not isinstance(data, (dict, list)):
            raise TypeError(f'Expected data to be a dict or list, but got {type(data)}')
        data = data

    return data

def read_xyz(path: str):
    """
    Read an xyz file and return the atoms and coordinates.
    :param path: Path to the xyz file
    :return: Atoms and coordinates
    """
    atoms = ase.io.read(path)
    coords = atoms.get_positions()
    elements = atoms.get_chemical_symbols()
    return elements, coords


if __name__ == '__main__':

    # metalig = '/Users/timosommer/PhD/projects/RCA/projects/DART/DARTassembler/data/metalig/MetaLigDB_v1.0.0.jsonlines'
    test_metalig = '/Users/timosommer/PhD/projects/RCA/projects/DART/DARTassembler/data/metalig/test1000_MetaLigDB_v1.0.0.jsonlines'
    # compress_file(metalig)
    # uncompress_file(metalig)

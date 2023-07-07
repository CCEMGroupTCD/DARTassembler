"""
Utility functions for input and output.
"""
import json
from src01.Molecule import RCA_Ligand, RCA_Complex
import numpy as np
from datetime import datetime
from src01.utilities import get_duration_string
from typing import Union
from pathlib import Path
import jsonlines
from tqdm import tqdm

class NumpyEncoder(json.JSONEncoder):
    """Special json encoder for numpy types. This is important to use in json.dump so that if json encounters a np.array, it converts it to a list automatically, otherwise errors arise. Use like this:
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
        return json.JSONEncoder.default(self, obj)


def load_json(path: Union[str, Path], n_max: int=None, show_progress: bool=True) -> dict:
    """
    Load a JSON or JSON Lines file. If the file is a JSON Lines file, it is converted to a dictionary.
    :param path: Path to the JSON or JSON Lines file
    :return: Dictionary with the contents of the file
    """
    db = {key: value for key, value in iterate_over_json(path, n_max=n_max, show_progress=show_progress)}

    return db

def check_if_return_entry(i: int, n_max: Union[int, list]=None) -> bool:
    """
    Check if the entry should be returned. It will not be returned if the index i is larger than n_max.
    """
    # Accept False, None or np.inf to disable n_max
    if n_max is None or n_max is False or n_max is np.inf:
        return True

    # If n_max is an integer, check if the index is smaller than n_max
    is_good = i < n_max

    return is_good

def iterate_over_json(path: Union[str, Path], n_max: int=None, show_progress: bool=True) -> tuple[str, dict]:
    """
    Iterate over a JSON or JSON Lines file and yield the key and value of each entry.
    :param path: Path to the JSON or JSON Lines file
    :return: Tuple with the key and value of each entry
    """
    try:
        # Try to load as normal JSON file first
        with open(path, 'r') as file:
            db = json.load(file)
            for i, (key, value) in enumerate(db.items()):
                if check_if_return_entry(i, n_max):
                    yield key, value
                else:
                    return

    except json.JSONDecodeError:
        # If normal JSON fails, try to load as JSON Lines
        with jsonlines.open(path, 'r') as reader:
            for i, line in tqdm(enumerate(reader), disable=not show_progress, desc='Load json'):
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

def load_unique_ligand_db_iteratively(path: Union[str, Path], molecule: str='dict', n_max=None, show_progress: bool=False) -> dict:
    check_molecule_value(molecule)  # Check if the molecule value is valid
    for name, mol in tqdm(iterate_over_json(path, n_max=n_max, show_progress=False), disable=not show_progress, desc='Load unique ligand db'):
        if molecule == 'class':
            mol = RCA_Ligand.read_from_mol_dict(mol)
            yield name, mol

def load_unique_ligand_db(path: Union[str, Path], molecule: str='dict', n_max=None, show_progress: bool=True) -> dict:
    db = {name: mol for name, mol in load_unique_ligand_db_iteratively(path=path, molecule=molecule, n_max=n_max, show_progress=show_progress)}
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

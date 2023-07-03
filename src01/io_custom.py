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


def load_json(path: Union[str, Path], n_max: int=None) -> dict:
    """
    Load a JSON or JSON Lines file. If the file is a JSON Lines file, it is converted to a dictionary.
    :param path: Path to the JSON or JSON Lines file
    :return: Dictionary with the contents of the file
    """
    # Accept False, None or np.inf to disable n_max
    if n_max is None or n_max is False:
        n_max = np.inf

    db = {}
    for i, (key, value) in enumerate(iterate_over_json(path)):
        if i >= n_max:
            break
        db[key] = value

    return db

def iterate_over_json(path: Union[str, Path]) -> tuple[str, dict]:
    """
    Iterate over a JSON or JSON Lines file and yield the key and value of each entry.
    :param path: Path to the JSON or JSON Lines file
    :return: Tuple with the key and value of each entry
    """
    try:
        # Try to load as normal JSON file first
        with open(path, 'r') as file:
            db = json.load(file)
            for key, value in db.items():
                yield key, value
    except json.JSONDecodeError:
        # If normal JSON fails, try to load as JSON Lines
        with jsonlines.open(path, 'r') as reader:
            for line in reader:
                # Since 'line' is a dictionary, no need to use json.loads
                yield line['key'], line['value']

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

def load_complex_db(path: Union[str, Path], molecule: str='dict', n_max=None) -> dict:
    start = datetime.now()

    check_molecule_value(molecule)
    db = load_json(path, n_max=n_max)

    if molecule == 'class':
        db = {name: RCA_Complex.read_from_mol_dict(mol) for name, mol in db.items()}

    duration = get_duration_string(start=start)
    print(f'Loaded complex db. Time: {duration}. ')
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

def load_unique_ligand_db(path: Union[str, Path], molecule: str='dict', n_max=None) -> dict:
    start = datetime.now()

    db = load_json(path, n_max=n_max)
    if molecule == 'class':
        db = {name: RCA_Ligand.read_from_mol_dict(mol) for name, mol in db.items()}

    duration = get_duration_string(start=start)
    print(f'Loaded unique ligand db. Time: {duration}. ')
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

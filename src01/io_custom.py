"""
Utility functions for input and output.
"""
import json
from src01.Molecule import RCA_Molecule, RCA_Ligand
from tqdm import tqdm


def load_json(path) -> dict:
    with open(path, 'r') as file:
        db = json.load(file)

    return db

def save_json(db: dict, path: str):
    with open(path, 'w') as file:
        json.dump(db, file)

    return

def check_molecule_value(output: str):
    possible_values = ['dict', 'class']
    if not output in possible_values:
        raise ValueError(f'Unknown value for `output`: {output}')

    return
def load_complex_db(path: str, molecule: str='dict') -> dict:
    check_molecule_value(molecule)
    db = load_json(path)
    if molecule == 'class':
        db = {name: RCA_Molecule.read_from_mol_dict(mol) for name, mol in db.items()}
    return db

def load_full_ligand_db(path: str, molecule: str='dict') -> dict:
    check_molecule_value(molecule)
    db = load_json(path)
    if molecule == 'class':
        db = {name: RCA_Ligand.read_from_mol_dict(mol) for name, mol in db.items()}
    return db

def load_unique_ligand_db(path: str) -> dict:
    return load_json(path)

def save_complex_db(db: dict, path: str):
    save_json(db, path)

def save_full_ligand_db(db: dict, path: str):
    save_json(db, path)

def save_unique_ligand_db(db: dict, path: str):
    save_json(db, path)

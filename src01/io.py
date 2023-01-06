"""
Utility functions for input and output.
"""
import json

def load_json(path) -> dict:
    with open(path, 'r') as file:
        db = json.load(file)

    return db

def save_json(db: dict, path: str):
    with open(path, 'w') as file:
        json.dump(db, file)

    return

def load_complex_db(path: str) -> dict:
    return load_json(path)

def load_ligand_db(path: str) -> dict:
    return load_json(path)

def load_unique_ligand_db(path: str) -> dict:
    return load_json(path)

def save_complex_db(db: dict, path: str):
    save_json(db, path)

def save_ligand_db(db: dict, path: str):
    save_json(db, path)

def save_unique_ligand_db(db: dict, path: str):
    save_json(db, path)

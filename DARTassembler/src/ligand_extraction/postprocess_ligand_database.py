"""
This file is for postprocessing the ligand database as it comes out of the ligand extraction.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
from pathlib import Path
from DARTassembler.src.constants.Paths import project_path
from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand
from DARTassembler.src.ligand_extraction.io_custom import iterate_unique_ligand_db

from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db

from DARTassembler.src.ligand_extraction.DataBase import LigandDB

# def create_ligand_ID(original_complex_id: str, denticity: int, other_IDs) -> str:
#     name = original_complex_id.lower().capitalize()
#     name += f'{denticity:.2d}'



if __name__ == '__main__':

    db_version = '1.7'
    db_path = project_path().extend(*f'data/final_db_versions/unique_ligand_db_v{db_version}.json'.split('/'))
    out_db_name  = 'MetaLigDB_v1.0.0.jsonlines'
    nmax = False

    db = {}
    for uname, moldict in iterate_unique_ligand_db(path=db_path, molecule='dict', n_max=nmax, show_progress=True):
        mol = RCA_Ligand.read_from_mol_dict(moldict)
        if not mol.pred_charge_is_confident or mol.denticity <= 0:
            continue
        db[uname] = mol

    db = LigandDB(db)
    outpath = Path(db_path.parent, out_db_name)
    db.to_json(path=outpath, json_lines=True)
    print(f'Saved ligand database with {len(db.db)} ligands to {outpath}.')















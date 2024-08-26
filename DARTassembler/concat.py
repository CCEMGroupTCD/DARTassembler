"""
This module concatenates multiple ligand databases into one.
"""
import warnings
from pathlib import Path
from typing import Union
from DARTassembler.src.ligand_extraction.DataBase import LigandDB


def concat(paths: list[str]) -> None:
    """
    Concatenate multiple ligand databases into one. The output will be saved in the current working directory as `concat_ligand_db.jsonlines`.
    :param paths: Paths to the ligand databases.
    """
    print(f"===============     CONCAT MODULE     =================")
    print(f"This module concatenates multiple ligand databases into one. The output ligand db will be saved in the current working directory as `concat_ligand_db.jsonlines`.")

    # Check if all paths are valid
    for path in paths:
        if not Path(path).exists():
            raise FileNotFoundError(f"File `{path}` not found.")

    # Load all ligand databases
    ligand_dbs = [LigandDB.load_from_json(path) for path in paths]

    # Print number of ligands in each database
    for i, db in enumerate(ligand_dbs):
        print(f"Ligand database {i+1} contains {len(db.db)} ligands.")

    # Concatenate ligand databases
    full_db = {}
    for i, db in enumerate(ligand_dbs):
        # Check if any ligands are duplicated
        if len(set(db.db.keys()).intersection(full_db.keys())) > 0:
            warnings.warn(f"Ligand database Nr {i} contains ligands already existing in previous databases. Duplicates will be overwritten. This should probably be ok since all ligands in the MetaLig are unique. If you want to avoid this warning, please remove duplicates from the ligand databases before concatenating them.")
        full_db.update(db.db)
    full_db = LigandDB(full_db)

    # Print number of ligands in the concatenated database
    print(f"The concatenated ligand database contains {len(full_db.db)} ligands.")

    # Save concatenated ligand database
    outpath = 'concat_ligand_db.jsonlines'
    full_db.to_json(outpath)

    print(f"Concatenated ligand databases saved to `{outpath}`.")
    print(f"Done! Exiting concat module.")


# Integration test, to check if everything is working and the same as before.
if __name__ == "__main__":
    from DARTassembler.src.constants.Paths import project_path

    print('Done!')

from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.metalig.mol import Ligand
from tqdm import tqdm
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['svg.fonttype'] = 'none'
import numpy as np


# Goals
# 1. show the variety of number of donor atoms for 1-mono
# 2. show a small chart of the number of unique donor atom configs for 2-cis

def filter_ligands_by_geometry(db: LigandDB, geometry: str):
    """
    Filter ligands by their geometry.
    """
    return {k: v for k, v in db.items() if v.geometry == geometry}


def filter_ligands_by_donor_atom_config(db: LigandDB):
    """
    Filter ligands by their donor atom configuration.
    :param db: The ligand database.
    :return: a dictionary with key unique donor atom configurations like 'N-N-C'
             and values as lists of Ligand objects.
    """
    donor_atom_configs = {}
    for ligand in db.values():
        donor_atoms = ligand.donor_elements  # e.g., ["N", "O", "S"]
        donor_atoms.sort()  # Sort for consistent representation
        donor_atoms_string = "-".join(donor_atoms)  # e.g., "C-N-N"

        if donor_atoms_string not in donor_atom_configs:
            donor_atom_configs[donor_atoms_string] = []
        donor_atom_configs[donor_atoms_string].append(ligand)

    return donor_atom_configs



if __name__ == "__main__":
    db = LigandDB.from_json(n_max=500).db

    # Filter ligands by geometry
    cis_bidentate_ligands = filter_ligands_by_geometry(db, "2-cis")

    print("done")

from typing import Tuple, List

import numpy as np
from openbabel import openbabel as ob

from DARTassembler.src.constants.Periodic_Table import DART_Element


def get_coordinates_and_elements_from_OpenBabel_mol(mol: ob.OBMol) -> Tuple[np.ndarray, List[str]]:
    """
    Returns the 3D coordinates of each atom in the molecule and the corresponding list of chemical elements.

    Args:
        mol (ob.OBMol): An Open Babel molecule object.

    Returns:
        Tuple[np.ndarray, List[str]]: A tuple containing a N x 3 numpy array of xyz coordinates and a list of chemical elements for each atom.
    """
    n_atoms = mol.NumAtoms()
    coords = np.empty((n_atoms, 3))
    elements = []

    for idx, atom in enumerate(ob.OBMolAtomIter(mol)):
        coords[idx, 0] = atom.GetX()
        coords[idx, 1] = atom.GetY()
        coords[idx, 2] = atom.GetZ()
        atomic_number = atom.GetAtomicNum()
        elements.append(DART_Element(atomic_number).symbol)

    return coords, elements

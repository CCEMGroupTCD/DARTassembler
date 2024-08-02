from typing import Tuple, List
import numpy as np
from openbabel import openbabel as ob
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from DARTassembler.src.constants.Periodic_Table import DART_Element
from openbabel import pybel
pybel.ob.obErrorLog.StopLogging()   # Remove Openbabel warnings


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


class ForceField(object):

    def __init__(self):
        """
        Initialize the UFF force field object.

        Args:
            backend (str, optional): The backend to use for the UFF force field. Defaults to 'openbabel'.
        """
        pass

    def singlepoint(self, mol):
        return self._singlepoint_with_openbabel(mol)

    def _singlepoint_with_openbabel(self, mol: ob.OBMol) -> float:
        """
        Calculate the energy of the molecule using the UFF force field with Open Babel.
        """
        if isinstance(mol, str):
            mol = pybel.readstring("mol", mol)
        obmol = mol.OBMol
        ff = ob._openbabel.OBForceField_FindType("uff")
        assert (ff.Setup(obmol))
        kj_to_kcal = 1.0 / 4.184
        ff.SetCoordinates(mol.OBMol)
        uffE = ff.Energy(False) * kj_to_kcal

        return uffE

    def optimize(self, mol_rdkit, fixed_atom_indices, nsteps=100):
        return self._optimize_with_openbabel(mol_rdkit, fixed_atom_indices, nsteps)

    def _optimize_with_openbabel(self, mol_rdkit: Chem.Mol, fixed_atom_indices: List[int], nsteps: int):

        # setup conversion
        xyz_string = rdmolfiles.MolToXYZBlock(mol_rdkit)
        conv = ob.OBConversion()
        conv.SetInAndOutFormats('xyz', 'xyz')
        mol = ob.OBMol()
        conv.ReadString(mol, xyz_string)

        # Define constraints
        constraints = ob.OBFFConstraints()
        for idx in fixed_atom_indices:
            constraints.AddAtomConstraint(1 + idx)  # The one is to account for open babel indexing starting at 1

        # Set up the force field with the constraints
        forcefield = ob.OBForceField.FindForceField("Uff")
        forcefield.Setup(mol, constraints)
        forcefield.SetConstraints(constraints)

        # Optimize the molecule coordinates using the force field with constrained atoms.
        optimized_coords = []
        optimized_elements = []
        forcefield.ConjugateGradientsInitialize(nsteps)
        while forcefield.ConjugateGradientsTakeNSteps(1):
            forcefield.GetCoordinates(mol)
            coords, elements = get_coordinates_and_elements_from_OpenBabel_mol(mol)
            optimized_coords.append(coords)
            optimized_elements.append(elements)
        forcefield.GetCoordinates(mol)
        xyz_string_output = conv.WriteString(mol)

        return xyz_string_output, optimized_coords, optimized_elements





from typing import Tuple, List
import numpy as np
import rdkit
import stk
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from DARTassembler.src.constants.Periodic_Table import DART_Element


def get_coordinates_and_elements_from_OpenBabel_mol(mol) -> Tuple[np.ndarray, List[str]]:
    """
    Returns the 3D coordinates of each atom in the molecule and the corresponding list of chemical elements.

    Args:
        mol (ob.OBMol): An Open Babel molecule object.

    Returns:
        Tuple[np.ndarray, List[str]]: A tuple containing a N x 3 numpy array of xyz coordinates and a list of chemical elements for each atom.
    """
    from openbabel import openbabel as ob
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

    def __init__(self, backend='openbabel'):
        """
        Initialize the force field object. This object is used to calculate the energy of a molecule using the UFF force field.
        """
        self.backend = backend


    def singlepoint(self, mol):
        if self.backend == 'openbabel':
            return self._singlepoint_with_openbabel(mol)
        elif self.backend == 'rdkit':
            return self._singlepoint_with_rdkit(mol)
        else:
            raise ValueError(f"Unsupported backend: {self.backend}")

    def _singlepoint_with_openbabel(self, mol) -> float:
        """
        Calculate the energy of the molecule using the UFF force field with Open Babel.
        :param mol: openbabel.OBMol
        :return: Energy in kcal/mol
        """
        from openbabel import openbabel as ob
        from openbabel import pybel
        pybel.ob.obErrorLog.StopLogging()  # Remove Openbabel warnings

        if isinstance(mol, str):
            mol = pybel.readstring("mol", mol)
        elif isinstance(mol, stk.BuildingBlock):
            string = stk.MolWriter().to_string(mol)
            mol = pybel.readstring("mol", string)
        else:
            raise ValueError(f"Unsupported input type: {type(mol)}")

        obmol = mol.OBMol
        ff = ob._openbabel.OBForceField_FindType("uff")
        assert (ff.Setup(obmol))
        kj_to_kcal = 1.0 / 4.184
        ff.SetCoordinates(mol.OBMol)
        uffE = ff.Energy(False) * kj_to_kcal

        return uffE

    def _singlepoint_with_rdkit(self, mol):
        """
        Calculate the energy of the molecule using the UFF force field with RDKit.
        :param mol: rdkit.Chem.Mol
        :return: Energy in kcal/mol
        """
        if isinstance(mol, stk.BuildingBlock):
            string = stk.MolWriter().to_string(mol)
            rdkit_mol = Chem.MolFromMolBlock(string, sanitize=False, removeHs=False, strictParsing=False)
        else:
            raise ValueError(f"Unsupported input type: {type(mol)}")

        # rdkit_mol = Chem.AddHs(rdkit_mol)
        # ff = AllChem.UFFGetMoleculeForceField(rdkit_mol)
        # ff.Initialize()

        # Adjust the valence of problematic atoms temporarily
        for atom in rdkit_mol.GetAtoms():
            if atom.GetExplicitValence() > atom.GetTotalValence():
                atom.SetNumExplicitHs(atom.GetTotalValence() - atom.GetExplicitValence())

        # Initialize the UFF force field
        ff = AllChem.UFFGetMoleculeForceField(rdkit_mol)
        ff.Initialize()

        # Revert valence adjustments
        for atom in rdkit_mol.GetAtoms():
            atom.SetNoImplicit(True)

        energy = ff.CalcEnergy()
        assert isinstance(energy, float)
        return energy

    def optimize(self, mol_rdkit, fixed_atom_indices, nsteps):
        return self._optimize_with_openbabel(mol_rdkit, fixed_atom_indices, nsteps)

    def _optimize_with_openbabel(self, mol_rdkit: Chem.Mol, fixed_atom_indices: List[int], nsteps: int):
        from openbabel import openbabel as ob
        # Read in the molecule and convert it to an Open Babel molecule
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





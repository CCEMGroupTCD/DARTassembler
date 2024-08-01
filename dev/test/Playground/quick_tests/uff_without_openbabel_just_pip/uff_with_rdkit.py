"""
A small script to test if the UFF force field is usable with rdkit when installing rdkit via pip, without openbabel.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from copy import deepcopy
import openbabel as ob
from rdkit.Chem import rdmolfiles

from DARTassembler.src.assembly.forcefields import get_coordinates_and_elements_from_OpenBabel_mol


def ob_uff_optimize(complex_mol, constrained_indices, nsteps=1000):
    # stk to xyz string
    xyz_string = rdmolfiles.MolToXYZBlock(complex_mol)

    # setup conversion
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('xyz', 'xyz')
    mol = ob.OBMol()
    conv.ReadString(mol, xyz_string)

    # Define constraints
    # OPEN BABEL INDEXING STARTS AT ONE !!!!!!!!!!!!!!!!!!!!!!!!!!!
    constraints = ob.OBFFConstraints()
    constraints.AddAtomConstraint(1)  # Here we lock the metal

    # Add constraints
    for atom_index in constrained_indices:
        constraints.AddAtomConstraint(1 + atom_index)  # The one is to account for open babel indexing starting at 1 and to account for the metal

    # Set up the force field with the constraints
    forcefield = ob.OBForceField.FindForceField("Uff")
    forcefield.Setup(mol, constraints)
    forcefield.SetConstraints(constraints)

    # Optimize the molecule coordinates using the force field with constrained atoms.
    forcefield.ConjugateGradientsInitialize(nsteps)
    while forcefield.ConjugateGradientsTakeNSteps(1):
        forcefield.GetCoordinates(mol)
    forcefield.GetCoordinates(mol)
    coords, elements = get_coordinates_and_elements_from_OpenBabel_mol(mol)

    return coords, elements


def get_positions(mol):
    mol_positions = []
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        mol_positions.append(positions)
    mol_positions = np.array(mol_positions)

    return mol_positions

if __name__ == '__main__':

    # Load SMILES as mol obj
    mol = Chem.MolFromSmiles('[Zn]([NH])([NH])([NH])N1C=C(Cl)N=C1')
    mol = AllChem.AddHs(mol) # make sure to add explicit hydrogens
    constrained_indices = (0, 1)
    nsteps = 1000

    AllChem.EmbedMolecule(mol)
    mol.GetConformer()
    mol_init = deepcopy(mol)

    # Get initial positions
    mol_init_positions = deepcopy(get_positions(mol))

    # Perform UFF optimization
    ff = AllChem.UFFGetMoleculeForceField(mol)
    ff.Initialize()
    for i in constrained_indices:
        ff.AddFixedPoint(i)
    ff.Minimize(maxIts=nsteps)
    optxyz_rdkit = Chem.AllChem.MolToXYZBlock(mol)
    # Get optimized positions
    optcoords_rdkit = get_positions(mol)

    # Perform UFF optimization with openbabel
    optcoords_ob, ob_elements = ob_uff_optimize(deepcopy(mol_init), constrained_indices, nsteps)

    # Check that the positions of fixed atoms are the same
    for i in constrained_indices:
        assert np.allclose(mol_init_positions[i], optcoords_rdkit[i])
    # Check that the positions of the other atoms have changed
    for i in range(len(mol_init_positions)):
        if i not in constrained_indices:
            assert not np.allclose(mol_init_positions[i], optcoords_rdkit[i])
    # Check that ob and rdkit uff optimization give the same results
    assert np.allclose(optcoords_rdkit, optcoords_ob)

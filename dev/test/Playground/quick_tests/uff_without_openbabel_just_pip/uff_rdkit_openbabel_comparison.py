from copy import deepcopy

from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
import pandas as pd
import numpy as np
from DARTassembler.src.assembly.forcefields import get_coordinates_and_elements_from_OpenBabel_mol

# Fix random and numpy seed for reproducibility
np.random.seed(0)
import random
random.seed(0)

energy_tol = 1e-6
force_tol = 1e-4

# Function for UFF optimization with RDKit
def optimize_with_rdkit(mol, nsteps, fixed_atom_indices):
    # mol = Chem.AddHs(mol)

    # Create the UFF force field and fix the specified atoms
    ff = AllChem.UFFGetMoleculeForceField(mol)
    ff.Initialize()
    for idx in fixed_atom_indices:
        ff.AddFixedPoint(idx)

    # Perform the optimization
    result = ff.Minimize(energyTol=energy_tol, forceTol=force_tol, maxIts=nsteps)
    if result != 0:
        print(f"RDKit UFF optimization did not converge. Result code: {result}")

    return mol


# Function for UFF optimization with Open Babel
def optimize_with_openbabel(mol, nsteps, fixed_atom_indices):
    obconversion = openbabel.OBConversion()
    obconversion.SetInAndOutFormats("mol", "xyz")

    obmol = openbabel.OBMol()
    # input_coords, _ = get_coordinates_and_elements_from_OpenBabel_mol(obmol)

    obconversion.ReadString(obmol, Chem.MolToMolBlock(mol))
    # obmol.AddHydrogens()

    ff = openbabel.OBForceField.FindForceField("uff")

    # Fix the specified atoms
    constraints = openbabel.OBFFConstraints()
    for idx in fixed_atom_indices:
        constraints.AddAtomConstraint(1 + idx)

    ff.Setup(obmol, constraints)
    ff.SetConstraints(constraints)

    # Perform the optimization
    ff.ConjugateGradients(nsteps, energy_tol)
    ff.GetCoordinates(obmol)

    # opt_coords, _ = get_coordinates_and_elements_from_OpenBabel_mol(obmol)
    # assert np.isclose(input_coords[fixed_atom_indices], opt_coords[fixed_atom_indices], atol=1e-4).all(), "The positions of the fixed atoms have changed!"
    # assert not np.isclose(input_coords, opt_coords, atol=1e-4).all(), "The positions of the other atoms have not changed!"


    return obmol


# Convert Open Babel molecule to XYZ string
def obmol_to_xyz(mol):
    obconversion = openbabel.OBConversion()
    obconversion.SetOutFormat("xyz")
    return obconversion.WriteString(mol)


# Save RDKit molecule to XYZ file
def save_rdkit_xyz(mol, filename):
    Chem.MolToXYZFile(mol, filename)


# Save Open Babel molecule to XYZ file
def save_openbabel_xyz(mol, filename):
    xyz_str = obmol_to_xyz(mol)
    with open(filename, 'w') as f:
        f.write(xyz_str)


# Compare the optimized geometries
def get_coordinates_from_rdkit(mol):
    conf = mol.GetConformer()
    return np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])


def get_coordinates_from_openbabel(mol):
    coords = []
    for atom in openbabel.OBMolAtomIter(mol):
        coords.append([atom.GetX(), atom.GetY(), atom.GetZ()])
    return np.array(coords)


if __name__ == '__main__':
    smiles = '[Zn]([NH])([NH])([NH])N1C=C(Cl)N=C1'
    nsteps = 85
    xyz_tol = 1e-2
    fixed_atom_indices = [0, 1, 2, 3]

    # Create the initial RDKit molecule with 3D coordinates
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    # Save the initial geometry to an XYZ file
    save_rdkit_xyz(mol, "initial.xyz")

    # Perform the optimizations
    rdkit_mol = optimize_with_rdkit(mol, nsteps, fixed_atom_indices)
    obabel_mol = optimize_with_openbabel(mol, nsteps, fixed_atom_indices)

    # Save the optimized geometries to XYZ files
    save_rdkit_xyz(rdkit_mol, "optimized_rdkit.xyz")
    save_openbabel_xyz(obabel_mol, "optimized_openbabel.xyz")

    coords_rdkit = get_coordinates_from_rdkit(rdkit_mol)
    coords_obabel = get_coordinates_from_openbabel(obabel_mol)

    # Assert that the optimized geometries are similar
    close = np.isclose(coords_rdkit, coords_obabel, atol=xyz_tol).all(axis=1)

    # Make pandas df to compare the coordinates
    df = pd.DataFrame(close, columns=['is_close'])
    for i in range(3):
        df[f'rdkit{i}'] = coords_rdkit[:, i]
        df[f'obabel{i}'] = coords_obabel[:, i]
        df[f'diff{i}'] = np.abs(coords_rdkit[:, i] - coords_obabel[:, i])

    assert close.all(), "The optimized geometries are not similar!"


    print("The optimized geometries are similar!")
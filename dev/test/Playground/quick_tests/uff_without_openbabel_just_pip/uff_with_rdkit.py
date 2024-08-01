"""
A small script to test if the UFF force field is usable with rdkit when installing rdkit via pip, without openbabel.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from copy import deepcopy


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

    AllChem.EmbedMolecule(mol)
    mol.GetConformer()

    # Get initial positions
    mol_init_positions = deepcopy(get_positions(mol))

    # Perform UFF optimization
    ff = AllChem.UFFGetMoleculeForceField(mol)
    ff.Initialize()
    for i in constrained_indices:
        ff.AddFixedPoint(i)
    ff.Minimize(energyTol=1e-7,maxIts=100000)

    xyz = Chem.AllChem.MolToXYZBlock(mol)

    # Get optimized positions
    mol_opt_positions = get_positions(mol)

    # Check that the positions of fixed atoms are the same
    for i in constrained_indices:
        assert np.allclose(mol_init_positions[i], mol_opt_positions[i])
    # Check that the positions of the other atoms have changed
    for i in range(len(mol_init_positions)):
        if i not in constrained_indices:
            assert not np.allclose(mol_init_positions[i], mol_opt_positions[i])

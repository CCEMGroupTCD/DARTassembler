"""
A small script to test if the UFF force field is usable with rdkit when installing rdkit via pip, without openbabel.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# Load SMILES as mol obj
mol = Chem.MolFromSmiles('[Zn]([NH])([NH])([NH])N1C=C(Cl)N=C1')
mol = AllChem.AddHs(mol) # make sure to add explicit hydrogens

# Generate 3D coordinates
AllChem.EmbedMolecule(mol)
xyz = Chem.AllChem.MolToXYZBlock(mol) # initial coords

# Perform UFF optimization
ff = AllChem.UFFGetMoleculeForceField(mol)
ff.Initialize()
ff.Minimize(energyTol=1e-7,maxIts=100000)

xyz = Chem.AllChem.MolToXYZBlock(mol)

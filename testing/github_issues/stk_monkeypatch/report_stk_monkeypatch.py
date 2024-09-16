

if __name__ == '__main__':

    import rdkit
    import stk
    testmol = stk.BuildingBlock.init_from_file('testmol.mol')
    rdkit_mol = testmol.to_rdkit_mol()
    rdkit.Chem.AllChem.SanitizeMol(rdkit_mol)
from src.process import LigandDatabase

if __name__ == '__main__':
    """
    Runs the full extraction
    TestSize: Decide if we want to run the programm on a testset (1/100) of the actual size
    Metals_of_interest: MetalAtom in TMC we want to extract the ligands from in the TMQM DB
    denticity_numbers: Denticities we are interested in
    
    """
    Testing = False

    ligand_db = LigandDatabase(TestSize=Testing)

    ligand_db.extract_ligands(denticity_numbers_of_interest=[2, 3, 4, 5],
                              metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"]
                              )

    ligand_db.add_monodentate_ligands()

    ligand_db.filter_duplicates()

    if Testing is True:
        pickle.dump(ligand_db, open("../data/ligand_db_test.pickle", "wb"))
    else:
        pickle.dump(ligand_db, open("../data/ligand_db.pickle", "wb"))

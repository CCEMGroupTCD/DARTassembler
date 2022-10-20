from src.LigandDatabase import LigandDatabase
import pickle
import pandas as pd
import numpy as np

if __name__ == '__main__':
    """
    Runs the full extraction
    
    TestSize: Decide if we want to run the programm on a testset (1/100) of the actual size
    Metals_of_interest: MetalAtom in TMC we want to extract the ligands from in the TMQM DB
    denticity_numbers: Denticities we are interested in
    
    """
    Testing = True

    ligand_db = LigandDatabase(TestSize=Testing)

    ligand_db.extract_ligands(denticity_numbers_of_interest=[2, 3, 4, 5],
                              metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"]
                              )

    ligand_db.add_monodentate_ligands()

    if Testing is True:
        pickle.dump(ligand_db, open("../Old/ligand_db_test.pickle", "wb"))
    else:
        pickle.dump(ligand_db, open("../data/LigandDatabases/ligand_db.pickle", "wb"))

    print("Ligand DB fully established")
    
    
    df = ligand_db.get_df_of_all_ligands().sort_values('name').reset_index(drop=True)
    df.to_csv('../data/ligand_db_test.csv', index=False)
    
    old_df = pd.read_csv('../data/Felix_original_ligand_db_test.csv').sort_values('name').reset_index(drop=True)
    pd.testing.assert_frame_equal(df[old_df.columns], old_df, check_like=True)
    print('Ligand database same as old database.')
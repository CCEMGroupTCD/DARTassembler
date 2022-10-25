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
    database_path = '../database/tmQM/data'
    id_col = 'CSD_code'
    atomic_properties_json = '../database/tmQM/data/atomic_properties/atomic_properties.json'

    ligand_db = LigandDatabase(data_path=database_path, id_col=id_col, TestSize=Testing)
    ligand_db.load_atomic_properties(atomic_properties_json)

    ligand_db.extract_ligands(denticity_numbers_of_interest=[1, 2, 3, 4, 5],
                              metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"]
                              )

    # ligand_db.add_monodentate_ligands()

    if Testing is True:
        pickle.dump(ligand_db, open("../data/LigandDatabases/ligand_db_test.pickle", "wb"))
    else:
        pickle.dump(ligand_db, open("../data/LigandDatabases/ligand_db.pickle", "wb"))

    print("Ligand DB fully established")
    
    df_extr_mols = ligand_db.get_Extracted_Molecules_global_information()
    
    df = ligand_db.get_df_of_all_ligands().sort_values('name').reset_index(drop=True)
    outpath = '../data/ligand_db_test.csv' if Testing else '../data/ligand_db.csv'
    df.to_csv(outpath, index=False)
    
    if Testing:
        old_outpath = '../data/221025_original_ligand_db_test_with_monodentates.csv'
        old_df = pd.read_csv(old_outpath).sort_values('name').reset_index(drop=True)
        pd.testing.assert_frame_equal(df[old_df.columns], old_df, check_like=True)
        print('Ligand database same as old database.')
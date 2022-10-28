from src.LigandDatabase import LigandDatabase
import pickle
import pandas as pd
import numpy as np
from pathlib import Path

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
    atomic_properties_json = Path('../database/tmQM/data/atomic_properties/atomic_properties.json')

    ligand_db = LigandDatabase(data_path=database_path, id_col=id_col, TestSize=Testing)
    
    if atomic_properties_json.exists():
        ligand_db.load_atomic_properties(atomic_properties_json)
    else:
        print('No json file of atomic properties found, so they will be read in from the .xyz files. You might want to save the atomic properties to json for faster future runs using `ligand_db.save_atomic_properties()`.')

    ligand_db.extract_ligands(denticity_numbers_of_interest=[1, 2, 3, 4, 5, 6],
                              metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"]
                              )
    ligand_db.filter_duplicates()
    
    Extracted_Molecules_path = '../data/LigandDatabases/all_Extracted_Molecules_test.json' if Testing else None
    ligand_db.save_Extracted_Molecules_to_json(Extracted_Molecules_path)

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
        old_outpath = '../data/221026_ligand_db_test_with_6-dentates.csv'
        old_df = pd.read_csv(old_outpath).sort_values('name').reset_index(drop=True)
        pd.testing.assert_frame_equal(df[old_df.columns], old_df, check_like=True)
        print('Ligand database same as old database.')
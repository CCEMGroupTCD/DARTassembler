from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from assembly import AssemblyComplex
from utilities import save_to_xyz
from ast import literal_eval
import pandas as pd
from ase import visualize




if __name__ == '__main__':

    # Define paths to the ligand databases and the OER complex csv
    oer_ligand_db_path = '/Users/cianclarke/Documents/PhD/OER_DATABASE/OERdatabase/data/testbatch/ligand_db/oer_all_ligands.jsonlines'
    OH_ligand_db_path = '/Users/cianclarke/Documents/PhD/OER_DATABASE/OERdatabase/data/testbatch/ligand_db/oer_OH.jsonlines'
    oer_complexes_csv = '/Users/cianclarke/Documents/PhD/OER_DATABASE/OERdatabase/data/testbatch/info_table.csv' # todo: CIAN UPDATE WITH TIMO'S CSV FOR RU UPDATE
    save_concat_xyz = 'concat_complexes.xyz'
    n_max = 1000000

    # Load ligand databases and append OH ligands to OER ligands to only have one database to search in. (all OER ligands)
    db = LigandDB.load_from_json(path=oer_ligand_db_path)
    OH_db = LigandDB.load_from_json(path=OH_ligand_db_path)
    db.db.update(OH_db.db)

    # Read in OER complex csv and extract the ligand names of the successfully assembled complexes.
    df = pd.read_csv(oer_complexes_csv)
    df = df[df['success']]    # todo this line needs to be uncommented with Timo's updated csv
    df['ligand names'] = df['ligand names'].apply(literal_eval)
    df = df.set_index('complex name')
    complexes = df.to_dict(orient='index')

    comments = []
    structures = []
    i = 1
    for complex_name, row in complexes.items():     # Loop through all the OER complexes
        # Extract the ligands from the database as old Ligand objects
        tridentate_name, bidentate_name, monodentate_name = row['ligand names']
        tridentate, bidentate, monodentate = db.db[tridentate_name], db.db[bidentate_name], db.db[monodentate_name]



        # Assemble the complex
        complex = AssemblyComplex(ligands=[monodentate, bidentate, tridentate], geometry='mer-3-2-1', metal='Ru')

        # Get the assembled complex as ASE Atoms object
        atoms = complex.assemble_complex()
        structures.append(atoms)
        comments.append(complex_name)
        #visualize.view(atoms) # Uncomment to visualize each assembled complex in ASE
        i = i + 1
        print(f"Complex index {i}")
        if len(structures) == n_max:
            break

    # Save all assembled complexes as concatenated xyz file so that they can be visualized in VESTA or Mercury
    save_to_xyz(outpath=save_concat_xyz, structures=structures, comments=comments)

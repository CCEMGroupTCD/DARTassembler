from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db
from pathlib import Path
import pandas as pd
from DARTassembler.src.ligand_extraction.utilities import unroll_dict_into_columns
from DARTassembler.src.constants.Paths import project_path

if __name__ == '__main__':

    denticity = 1
    coordinating_atoms = ['C']  # order doesn't matter
    charge = -1
    save_xyz_path = '../data/C_ligands_for_Manting'
    db_version = '1.7'
    ligand_db_path = project_path().extend(*f'data/final_db_versions/unique_ligand_db_v{db_version}.json'.split('/'))
    n_max = False

    db = LigandDB.from_json(ligand_db_path, molecule='class', max_number=n_max).db

    # Filter correct ligands
    correct_ligands = {}
    for name, mol in db.items():
        correct = (mol.denticity == denticity) and (sorted(mol.local_elements) == sorted(coordinating_atoms)) and (mol.pred_charge == charge) and mol.pred_charge_is_confident
        if correct:
            correct_ligands[name] = mol

    # Save all wanted ligands as .xyz
    save_xyz_path = Path(save_xyz_path).resolve()
    save_xyz_path.mkdir(parents=True, exist_ok=True)
    print(f'Saving {len(correct_ligands)} ligands to {save_xyz_path} as individual .xyz files.')
    for i, (name, mol) in enumerate(correct_ligands.items()):
        save_file = Path(save_xyz_path, f'{i}.xyz')
        xyz_string = mol.get_xyz_file_format_string(comment=name)
        with open(save_file, 'w') as file:
            file.write(xyz_string)

    # Save csv with these ligands and some information.
    df = pd.DataFrame.from_dict(
                                {name: mol.write_to_mol_dict() for name, mol in correct_ligands.items()},
                                orient='index')
    df = df.drop(columns=['graph_dict', 'atomic_props'])
    df = unroll_dict_into_columns(df, dict_col='global_props', prefix='gbl_', delete_dict=True)
    df = unroll_dict_into_columns(df, dict_col='stats', prefix='stats_', delete_dict=True)
    df.to_csv(Path(save_xyz_path, 'C_ligands_for_Manting.csv'))
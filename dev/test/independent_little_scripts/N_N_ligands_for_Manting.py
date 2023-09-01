from src01.io_custom import load_unique_ligand_db
from pathlib import Path
import pandas as pd
from src01.utilities import unroll_dict_into_columns

if __name__ == '__main__':

    denticity = 2
    coordinating_atoms = ['N', 'N']  # order doesn't matter
    only_allow_elements = ['C', 'H', 'N', 'Si']
    must_have_elements = ['C', 'H']
    charge = -1
    only_confident_charges = True
    save_xyz_path = '../data/N_N_charge_1_ligands_for_Manting'
    ligand_db_version = '1.5'

    db = load_unique_ligand_db(f'../../data/final_db_versions/unique_ligand_db_v{ligand_db_version}.json', molecule='class')
    print('Read in ligands.')

    #%% Filter correct ligands
    correct_ligands = {}
    for name, mol in db.items():
        correct_denticity = mol.denticity == denticity
        correct_coordinating_atoms = sorted(mol.local_elements) == sorted(coordinating_atoms)
        correct_charge = (mol.pred_charge == charge) and (mol.pred_charge_is_confident if only_confident_charges else True)
        correct_elements = (not any([el not in only_allow_elements for el in mol.atomic_props['atoms']])) and all([el in mol.atomic_props['atoms'] for el in must_have_elements])
        correct = correct_denticity and correct_coordinating_atoms and correct_charge and correct_elements
        if correct:
            correct_ligands[name] = mol

    # Save all wanted ligands as .xyz
    save_xyz_path = Path(save_xyz_path)
    save_xyz_path.mkdir(parents=True, exist_ok=True)
    for i, (name, mol) in enumerate(correct_ligands.items()):
        save_file = Path(save_xyz_path, f'{i}.xyz')
        # Sort coordinating atoms to be the first ones in atomic props
        mol.sort_atomic_props_to_have_coordinating_atoms_first()
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
    df.to_csv(Path(save_xyz_path, 'Ligands_for_Manting.csv'))

    print('Done!')
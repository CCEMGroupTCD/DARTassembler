from DARTassembler.src.metalig.db import LigandDB
from pathlib import Path
import pandas as pd
from DARTassembler.src.metalig.utils import unroll_dict_into_columns
from DARTassembler.src.constants.paths import project_path

if __name__ == '__main__':

    denticity = 1
    composition_only = ['C', 'H']  # order doesn't matter
    charge = -1
    save_xyz_path = '../data/CH_ligands_for_Ana'
    db_version = '1.7'
    ligand_db_path = project_path().extend(*f'data/final_db_versions/unique_ligand_db_v{db_version}.json'.split('/'))
    n_max = False

    db = LigandDB.from_json(ligand_db_path, molecule='class', max_number=n_max).db
    db = {key: mol for key, mol in db.items() if mol.n_donors > 0 and mol.has_confident_charge}

    # Filter correct ligands
    correct_ligands = {}
    filtered_out = []
    for name, mol in db.items():
        # if atoms.n_donors <= 0 or (not atoms.has_confident_charge):
        #     continue
        correct_denticity = mol.n_donors == denticity
        correct_composition = set(mol.atomic_props['atoms']) == set(composition_only)
        correct_charge = mol.charge == charge and mol.has_confident_charge
        correct_molweight = mol.global_props['molecular_weight'] > 400
        correct_n_Carbons = True#atoms.atomic_props['atoms'].count('C') == 40
        if not correct_composition:
            filtered_out.append('composition')
            continue
        if not correct_denticity:
            filtered_out.append('n_donors')
            continue
        if not correct_charge:
            filtered_out.append('charge')
            continue
        if not correct_molweight:
            filtered_out.append('molweight')
            continue

        correct_ligands[name] = mol

    filtered_out = pd.value_counts(filtered_out)
    print(f'Filtered out ligands from the initial {len(db)} ligands:')
    print(filtered_out)

    # Save all wanted ligands as .xyz
    save_xyz_path = Path(save_xyz_path).resolve()
    save_xyz_path.mkdir(parents=True, exist_ok=True)
    print(f'Saving {len(correct_ligands)} ligands to {save_xyz_path} as individual .xyz files.')
    for i, (name, mol) in enumerate(correct_ligands.items()):
        save_file = Path(save_xyz_path, f'{i}.xyz')
        xyz_string = mol.get_xyz_string(comment=name)
        with open(save_file, 'w') as file:
            file.write(xyz_string)

    # Save csv with these ligands and some information.
    df = pd.DataFrame.from_dict(
                                {name: mol.to_dict() for name, mol in correct_ligands.items()},
                                orient='index')
    df = df.drop(columns=['graph_dict', 'atomic_props'])
    df = unroll_dict_into_columns(df, dict_col='global_props', prefix='gbl_', delete_dict=True)
    df = unroll_dict_into_columns(df, dict_col='stats', prefix='stats_', delete_dict=True)
    df.to_csv(Path(save_xyz_path, 'CH_ligands_for_Ana.csv'))
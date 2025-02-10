from copy import deepcopy
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from DARTassembler.src.ligand_extraction.DataBase import LigandDB

# todo refactor metalig and DART:
# MetaLig:
# DART:


if __name__ == '__main__':
    n_max = None

    outpath = 'data/ligand_db_v1_1_0.jsonlines'
    db = LigandDB.load_from_json(path='metalig', n_max=n_max)
    data = []
    for uname, ligand in tqdm(db.db.items(), desc='Refactoring ligands'):
        eff_atoms, isomers_eff_idc = ligand.get_isomers_effective_ligand_atoms_with_effective_donor_indices()
        for i, eff_idc in enumerate(isomers_eff_idc):
            donor_positions = list(eff_atoms[eff_idc].get_positions())
            data.append([f'{uname}-{i}', ligand.geometry, ligand.has_neighboring_coordinating_atoms, ligand.n_eff_denticities, ligand.n_haptic_atoms, ligand.n_denticities, ligand.denticity, eff_atoms.get_chemical_formula(), eff_idc, donor_positions])
    df = pd.DataFrame(data, columns=['uname', 'geometry', 'haptic', 'n_eff_denticities', 'n_haptic_atoms', 'n_denticities', 'denticity', 'stoichiometry', 'eff_idc', 'donor_positions'])
    df_haptic = df[df['haptic']]
    df_haptic.sort_values(by=['geometry', 'uname'], inplace=True)

    db.save_to_file(outpath)
    df0 = db.save_reduced_csv(Path(outpath).with_suffix('.csv'))
    df1 = db.get_ligand_output_df()

    # Read in again to check if the data is saved correctly
    if n_max is not None and n_max < 500:
        db2 = LigandDB.load_from_json(path=outpath)
        df2 = db2.get_ligand_output_df()
        pd.testing.assert_frame_equal(df1, df2)
        import filecmp
        new_outpath = str(Path(outpath)).replace('.jsonlines', '_new.jsonlines')
        db2.save_to_file(new_outpath)
        assert filecmp.cmp(outpath, new_outpath)
        print('Data is saved correctly.')
    # cols = set(df1.columns).union(set(df2.columns))
    # for col in cols:
    #     if not pd.testing.assert_series_equal(df1[col], df2[col]):
    #         print(f'Column {col} is not equal!')

    # from DARTassembler.src.ligand_extraction.io_custom import load_jsonlines
    # db_dict = load_jsonlines(default_ligand_db_path, n_max=n_max)
    # ligand_old = db_dict['unq_CSD-OZIYON-02-a']
    # ligand_new = refactor_metalig_entry_from_v1_0_0_to_v1_1_0(ligand_old)

    print('Done')

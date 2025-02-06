from copy import deepcopy
import pandas as pd
from DARTassembler.src.ligand_extraction.DataBase import LigandDB

# todo refactor metalig and DART:
# MetaLig:
# - remove denticity and instead use kappa, eta, elcn, n_donors etc.
# - rename original_complex_indices in atomic_props to parent_complex_idc
# - improve stoichiometry string. Do this by deleting the stoichiometry property and re-calculating it from the atomic_props
# - add properties
#   - kappa
#   - eta
#   - elcn
#   - hapdent_idc
#   - isomers_hapdent_idc
#   - ligand_planarity
#   - donor_planarity
# DART
# refactor ligandinfo output to include only relevant output


if __name__ == '__main__':
    n_max = 10000

    outpath = 'ligand_db_v1_1_0.jsonlines'
    db = LigandDB.load_from_json(path='metalig', n_max=n_max)
    data = []
    for uname, ligand in db.db.items():
        eff_atoms, isomers_eff_idc = ligand.get_isomers_effective_ligand_atoms_with_effective_donor_indices()
        for i, eff_idc in enumerate(isomers_eff_idc):
            donor_positions = list(eff_atoms[eff_idc].get_positions())
            data.append([f'{uname}-{i}', ligand.geometry, ligand.has_neighboring_coordinating_atoms, ligand.elcn, ligand.eta, ligand.kappa, ligand.denticity, eff_atoms.get_chemical_formula(), eff_idc, donor_positions])
    df = pd.DataFrame(data, columns=['uname', 'geometry', 'haptic', 'elcn', 'eta', 'kappa', 'denticity', 'stoichiometry', 'eff_idc', 'donor_positions'])
    df_haptic = df[df['haptic']]
    df_haptic.sort_values(by=['geometry', 'uname'], inplace=True)

    # db.save_to_file(outpath)

    # from DARTassembler.src.ligand_extraction.io_custom import load_jsonlines
    # db_dict = load_jsonlines(default_ligand_db_path, n_max=n_max)
    # ligand_old = db_dict['unq_CSD-OZIYON-02-a']
    # ligand_new = refactor_metalig_entry_from_v1_0_0_to_v1_1_0(ligand_old)

    print('Done')

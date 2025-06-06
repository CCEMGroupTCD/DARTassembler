import pandas as pd
from DARTassembler.src.constants.paths import project_path
from DARTassembler.src.misc.io import load_unique_ligand_db


if __name__ == '__main__':

    db_version = '1.7'
    db_path = project_path().extend(*f'data/final_db_versions/unique_ligand_db_v{db_version}.json'.split('/'))
    exclude_unconnected_ligands = True
    exclude_uncertain_charges = True
    nmax = 1000

    important_cols = {'stoichiometry': 'Stoichiometry', 'donor_elements': 'Donors', 'n_donors': 'Denticity', 'charge': 'Formal Charge',  'n_atoms': 'Num. Atoms', 'n_electrons': 'Num. Electrons', 'parent_complex_id': 'CSD Complex ID', 'original_metal_symbol': 'CSD Metal', 'original_metal_os': 'CSD Metal OS', 'has_betaH': 'Beta Hydrogen', 'is_haptic': 'Haptic',  'n_ligand_instances': 'CSD Occurrences'}

    df = pd.DataFrame.from_dict(load_unique_ligand_db(path=db_path ,n_max=nmax), orient='index')
    df = df.query('has_confident_charge == True and n_donors > 0')
    df = df[important_cols.keys()]
    df = df.rename(columns=important_cols)

    df['Donors'] = df['Donors'].apply(lambda x: ', '.join(x))
    df['Formal Charge'] = df['Formal Charge'].astype(int)


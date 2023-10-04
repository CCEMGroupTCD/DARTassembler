import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
from DARTassembler.src.constants.Paths import project_path
from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db


if __name__ == '__main__':

    db_version = '1.7'
    db_path = project_path().extend(*f'data/final_db_versions/unique_ligand_db_v{db_version}.json'.split('/'))
    exclude_unconnected_ligands = True
    exclude_uncertain_charges = True
    nmax = 1000

    important_cols = {'stoichiometry': 'Stoichiometry', 'local_elements': 'Donors', 'denticity': 'Denticity', 'pred_charge': 'Formal Charge',  'n_atoms': 'Num. Atoms', 'n_electrons': 'Num. Electrons', 'original_complex_id': 'CSD Complex ID', 'original_metal_symbol': 'CSD Metal', 'original_metal_os': 'CSD Metal OS', 'has_betaH': 'Beta Hydrogen', 'has_neighboring_coordinating_atoms': 'Haptic',  'occurrences': 'CSD Occurrences'}

    df = pd.DataFrame.from_dict(load_unique_ligand_db(path=db_path ,n_max=nmax), orient='index')
    df = df.query('pred_charge_is_confident == True and denticity > 0')
    df = df[important_cols.keys()]
    df = df.rename(columns=important_cols)

    df['Donors'] = df['Donors'].apply(lambda x: ', '.join(x))
    df['Formal Charge'] = df['Formal Charge'].astype(int)


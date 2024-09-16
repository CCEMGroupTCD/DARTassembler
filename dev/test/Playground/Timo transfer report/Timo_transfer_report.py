import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
from DARTassembler.src.constants.Paths import project_path, default_ligand_db_path
from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db
try:    # Avoid error when running on server
    matplotlib.use('TkAgg')
except ImportError:
    pass
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none' # for correct text rendering in some programs
plt.rcParams['savefig.facecolor'] = 'white'

from pathlib import Path
sns.set_theme(style='ticks')






if __name__ == '__main__':

    db_version = '1.7'
    db_path = project_path().extend(*f'data/final_db_versions/unique_ligand_db_v{db_version}.json'.split('/'))
    exclude_unconnected_ligands = True
    exclude_uncertain_charges = True
    nmax = False


    # ligands = LigandDB.load_from_json(db_path, n_max=nmax)
    # df_ligands = pd.DataFrame.from_dict({name: lig.get_ligand_output_info(add_confident_charge=True) for name, lig in ligands.db.items()}, orient='index'
    metalig_csv_path = Path(default_ligand_db_path).with_suffix('.csv')
    df_metalig = pd.read_csv(metalig_csv_path, index_col=0)


    ####                         Full ligand database                        ####
    # ========================================================================= #
    df_all_ligands = pd.read_csv(db_path.with_suffix('.csv'), index_col=0)
    df_ligands = df_all_ligands.query('Denticity > 0')
    df_ligands_with_charge = df_ligands[df_ligands['Formal Charge'].notna()]
    df_all_ligands_with_charge = df_all_ligands[df_all_ligands['Formal Charge'].notna()]

    # print output stats for all ligands, including solvent molecules and counter ions
    print('All ligands including solvents:')
    n_denticities = df_all_ligands['Denticity'].value_counts()
    print(f'Denticities:\n{n_denticities}')
    print(f'Total number of ligands: {len(df_all_ligands)}')
    n_solvent = (df_all_ligands['Denticity'] <= 0).sum()
    print(f'Number of solvent molecules/counter ions: {n_solvent}')
    n_no_charge = df_all_ligands['Formal Charge'].isna().sum()
    print(f'Number of ligands with no charge: {n_no_charge}')
    n_no_confident_charge = (~df_all_ligands_with_charge['Charge Confident']).sum()
    print(f'Number of ligands with no confident charge: {n_no_confident_charge}')
    print()

    # print output stats for ligands only
    print('Only inner-sphere ligands:')
    n_no_charge = df_ligands['Formal Charge'].isna().sum()
    print(f'Number of ligands with no charge: {n_no_charge}')
    n_no_confident_charge = (~df_ligands_with_charge['Charge Confident']).sum()
    print(f'Number of ligands with no confident charge: {n_no_confident_charge}')

    # METALIG

    # histogram of denticity
    plt.figure()
    sns.histplot(data=df_metalig, x='Denticity', discrete=True)
    plt.yscale('log')
    plt.savefig('figures/metalig_denticity_hist.png', dpi=300)

    # histogram of formal charge
    plt.figure()
    sns.histplot(data=df_metalig, x='Formal Charge', discrete=True)
    plt.savefig('figures/metalig_charge_hist.png', dpi=300)






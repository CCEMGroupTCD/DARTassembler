"""
Playground which reads in the unique ligand db and lets you play with it.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
from DARTassembler.src.constants.Paths import project_path

from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db

matplotlib.use('TkAgg')
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
sns.set_theme()



def calculate_n_321_complexes_with_OH_fixed(charges_2: np.array, charges_3: np.array, metals: list[str],
                                            filter_factor: float) -> int:
    """
    Calculate the number of possible complexes for the 321 topology with one OH fixed.
    @param charges_2: 1D array of charges of 2-dentate ligands
    @param charges_3: 1D array of charges of 3-dentate ligands
    @param metals: list of metals to take into account
    @param filter_factor: factor to take into account which reduces the number of complexes due to pre and post filtering in the assembly process
    @return: number of possible complexes
    """
    if not charges_2.ndim == 1:
        raise ValueError('Charges of 2-dentate ligands must be 1D array.')
    if not charges_3.ndim == 1:
        raise ValueError('Charges of 3-dentate ligands must be 1D array.')

    charge_OH = -1
    charges_2 = charges_2.reshape((-1, 1))
    charges_3 = charges_3.reshape((1, -1))
    charge_combinations = charges_2 + charges_3 + charge_OH

    print('Calculate number of possible complexes for 321 topology with one OH fixed:')
    n_complexes = 0
    for metal in metals:
        for os in DART_Element(metal).common_oxidation_states:
            n_possible_complexes = (charge_combinations == - os).sum()
            n_complexes += n_possible_complexes
            print(f'\t{metal} {os}+: {n_possible_complexes:e}')

    n_complexes *= filter_factor
    print(f'Total number of possible 321 complexes with filter factor {filter_factor:.1g}: {n_complexes:e}')
    return n_complexes


if __name__ == '__main__':

    db_version = '1.7'
    db_path = project_path().extend(*f'data/final_db_versions/unique_ligand_db_v{db_version}.json'.split('/'))
    exclude_unconnected_ligands = True
    exclude_uncertain_charges = True
    nmax = 10000

    ligands = LigandDB.from_json(db_path, max_number=nmax)
    df_ligands = pd.DataFrame.from_dict(load_unique_ligand_db(path=db_path ,n_max=nmax), orient='index')
    df_ligands = df_ligands[df_ligands['denticity'] == 2]
    df_ligands = df_ligands[df_ligands['has_neighboring_coordinating_atoms']]

    planarities = {}
    for name, lig in ligands.db.items():
        if lig.denticity == 2 and not lig.has_neighboring_coordinating_atoms:
            is_planar, max_dist = lig.check_bidentate_planarity(return_tol=True)
            planarities[name] = {'is_planar': is_planar, 'max_dist': max_dist, 'stoi': lig.stoichiometry}
    df_planarities = pd.DataFrame.from_dict(planarities, orient='index')
    df_ligands = df_planarities.join(df_ligands)


    # ligands = load_unique_ligand_db(path=db_path, molecule='class')
    # #%%
    # df = pd.DataFrame.from_dict(load_unique_ligand_db(path=db_path), orient='index')
    #
    # if exclude_unconnected_ligands:
    #     df = df[df['denticity'] > 0]
    # if exclude_uncertain_charges:
    #     df = df[df['pred_charge_is_confident']]
    # df = unroll_dict_into_columns(df, dict_col='global_props', prefix='gbl_', delete_dict=True)
    # df = unroll_dict_into_columns(df, dict_col='stats', prefix='stats_', delete_dict=True)

    #%%
    # add_props = pd.DataFrame.from_dict({name: {'has_BetaH': lig.has_betaH, 'only_CHNO': len(set(['C', 'H', 'N', 'O'] + lig.atomic_props['atoms'])) <= 4, 'only_ON_donors': len(set(['N', 'O']+lig.local_elements))<=2} for name, lig in ligands.items()}, orient='index')
    # df = df.join(add_props)
    # df['mono_and_small'] = (df['denticity'] == 1) & (df['gbl_n_atoms'] <= 20)
    #
    # df = df.query('not has_BetaH and only_CHNO and only_ON_donors')
    # df = df[((df['denticity'] > 1) & (df['denticity'] < 6)) | df['mono_and_small']]
    #
    # groups = df.groupby(['denticity', 'pred_charge']).size().reset_index().rename(columns={0: 'count'})


    #%% Start playing here

    # Count possible complexes by combining denticity and charge
    # This assumes that we have 10 metals in oxidation state 3, which is already connected to one OH ligand.
    # The other two ligands are assumed to be a bidentate and a tridentate which account together for a charge of -2.
    # metals = ['Cr', 'Mn', 'Fe', 'Ru', 'Co', 'Ni']
    # filter_factor = 0.5         # Factor to account for the pre and post filtering in the assembly
    # charges_2 = df.query('denticity == 2')['pred_charge'].to_numpy(dtype=int)
    # charges_3 = df.query('denticity == 3')['pred_charge'].to_numpy(dtype=int)
    # n_complexes = calculate_n_321_complexes_with_OH_fixed(charges_2, charges_3, metals, filter_factor)


    # n_metals = 10               # Each metal of first row of transition metals, assuming os = 3
    # n_monodentates_1 = 3        # Number of self-defined monodentates with charge = -1
    # n_monodentates_0 = 3        # Number of self-defined monodentates with charge = 0
    # N_3_2 = len(df.query('denticity == 3 and pred_charge == -2'))
    # N_2_0 = len(df.query('denticity == 2 and pred_charge == 0'))
    # N_3_1 = len(df.query('denticity == 3 and pred_charge == -1'))
    # N_2_1 = len(df.query('denticity == 2 and pred_charge == -1'))
    # N_3_0 = len(df.query('denticity == 3 and pred_charge == 0'))
    # N_2_2 = len(df.query('denticity == 2 and pred_charge == -2'))
    # N_top321_os3 = N_3_2*N_2_0 + N_3_1*N_2_1 + N_3_0*N_2_2
    # N_top2211_os3 = N_2_1 * N_2_0 * n_monodentates_1 + N_2_2 * N_2_0 * n_monodentates_0
    # n_complexes = n_metals * (N_top321_os3 + N_top2211_os3) * filter_factor
    # print(f'Number of possible complexes: {n_complexes:e}')

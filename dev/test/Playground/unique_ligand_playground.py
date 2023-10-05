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
    nmax = 100

    ligands = LigandDB.from_json(db_path, max_number=nmax)
    for lig in ligands.db.values():
        lig.check_bidentate_planarity()
    # df_ligands = pd.DataFrame.from_dict(load_unique_ligand_db(path=db_path ,n_max=nmax), orient='index')
    # df_ligands = df_ligands[(df_ligands['denticity'] > 0) & df_ligands['pred_charge_is_confident']]
    # # df_ligands = df_ligands.groupby(['heavy_atoms_graph_hash_with_metal', 'pred_charge']).agg(lambda subdf: )
    # z = df_ligands[df_ligands['atomic_props'].apply(lambda props: 'C' in props['atoms'] and not 'H' in props['atoms'])]
    # z_all = df_ligands[df_ligands['atomic_props'].apply(lambda props: len(pd.unique(props['atoms'])) == 1) & (df_ligands['n_atoms'] >= 3)]
    # df = df_ligands[df_ligands['stoichiometry'] == 'C6']

    # # Get number of possible ligands:
    # ligands.filter_exclude_unconnected_ligands()
    # ligands.filter_not_charge_confident_ligands()
    # ligands.filter_ligands_with_neighboring_coordinating_atoms()
    # ligands.filter_non_centrosymmetric_monodentates()
    # ligands.filter_n_atoms(max_n_atoms=15, denticities=[1])
    # df_num_possible_complexes = ligands.calc_number_of_possible_complexes()
    # df = ligands.get_df_of_all_ligands()




    print('Done')

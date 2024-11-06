"""
Playground which reads in the unique ligand db and lets you play with it.
"""
from pathlib import Path
import pandas as pd
import numpy as np
# from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db, load_jsonlines, load_json
from tqdm import tqdm



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
    from DARTassembler.src.constants.Periodic_Table import DART_Element

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

    # db_version = '1.7'
    # db_path = project_path().extend(*f'data/final_db_versions/unique_ligand_db_v{db_version}.json'.split('/'))
    # exclude_unconnected_ligands = True
    # exclude_uncertain_charges = True
    # nmax = False
    #
    # # for thesis
    # df_ligands = pd.read_csv(db_path.with_suffix('.csv'), index_col=0)
    # n_denticities = df_ligands['Denticity'].value_counts()
    # print(f'Denticities:\n{n_denticities}')
    # print(f'Total number of ligands: {len(df_ligands)}')
    # n_solvent = (df_ligands['Denticity'] <= 0).sum()
    # print(f'Number of solvent molecules/counter ions: {n_solvent}')


    # ===== Code for MetaLig documentation =====
    from DARTassembler.src.ligand_extraction.DataBase import LigandDB

    test_path = '/Users/timosommer/PhD/projects/DARTassembler/testing/github_issues/#4 metalig jsonlines cannot be read in on windows/MetaLigDB_v1.0.0.jsonlines'

    # Load the first 1000 out of 41,018 ligands in the MetaLig database.
    metalig = LigandDB.load_from_json(path=test_path)

    # # Get an overview of all tridentate ligands that were used in projects to try to not make them different in planarity.
    # used_ligands_paths = {
    #     'OER': '/Users/timosommer/PhD/projects/OERdatabase/data/testbatch/ligand_db/oer_all_ligands.jsonlines'
    # }
    # used_ligands = {}
    # for project in used_ligands_paths.keys():
    #     db_path = used_ligands_paths[project]
    #     db_used = LigandDB.load_from_json(path=db_path)
    #     for uname, ligand in db_used.db.items():
    #         if ligand.denticity == 3:
    #             used_ligands[uname] = project
    # df_all_used_ligands = pd.DataFrame.from_dict(used_ligands, orient='index', columns=['where'])
    #
    #
    # data = {}
    # ligand_names = metalig.db.keys()
    # for uname in tqdm(ligand_names, desc='Calculating donor planarity'):
    #     ligand = metalig.db[uname]
    #     new_planar = None
    #     if ligand.denticity == 3:
    #         new_planar = ligand.calculate_donors_planarity(with_metal=True)
    #     # if ligand.denticity == 4:
    #     #     new_planar = ligand.calculate_donors_planarity(with_metal=True)
    #
    #     if new_planar is not None:
    #         data[uname] = {
    #         'denticity': ligand.denticity,
    #         # 'old_planar': ligand.planar_check(),
    #         'new_planar': new_planar,
    #         }
    # df_planarity = pd.DataFrame.from_dict(data, orient='index')
    #
    # df = metalig.get_ligand_output_df()
    # df = df_planarity.join(df, how='inner')

    # # Save to .jsonlines file
    # outfile = Path('/Users/timosommer/Downloads/test_DART/speedup_metalig/data_output/metalig_3000.jsonlines')
    # metalig.save_to_file(outfile)
    #
    # #%% ==============    Doublecheck refactoring    ==================
    # from dev.test.Integration_Test import IntegrationTest
    # old_dir = Path(outfile.parent.parent, 'benchmark_data_output')
    # if old_dir.exists():
    #     test = IntegrationTest(new_dir=outfile.parent, old_dir=old_dir)
    #     test.compare_all()
    #     print('Test for assembly of complexes passed!')
    # else:
    #     print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')


    # dict_metalig = load_jsonlines(default_ligand_db_path)

    #
    # # Set some criteria to filter ligands
    # keep_denticity = 2
    # keep_charge = -1
    # max_n_atoms = 50
    #
    # ligands_to_keep = []
    # for ligand_name, ligand in metalig.db.items():
    #     correct_denticity = ligand.denticity == keep_denticity
    #     correct_charge = ligand.pred_charge == keep_charge
    #     correct_n_atoms = ligand.n_atoms <= max_n_atoms
    #     if correct_denticity and correct_charge and correct_n_atoms:
    #         ligands_to_keep.append(ligand_name)
    #
    # # Reduce MetaLig database to only keep ligands which adhere to the above criteria
    # filtered_metalig_dict = {ligand_name: ligand for ligand_name, ligand in metalig.db.items() if ligand_name in ligands_to_keep}
    # filtered_metalig = LigandDB(filtered_metalig_dict)
    #
    # # Save filtered MetaLig database as .jsonlines file.
    # filtered_metalig.save_to_file('filtered_metalig.jsonlines')
    #
    # # Save an overview table of the filtered ligand database as csv file.
    # filtered_metalig.save_reduced_csv('filtered_metalig.csv')










    # df_ligands['has_metal_neighbors'] = df_ligands['graph_dict'].apply(lambda graph_dict: 'metal_neighbor' in str(graph_dict))
    #
    # planarity = [ligand.calculate_planarity() for ligand in ligands.db.values()]
    # df_ligands['planarity'] = planarity

    # has_identical_ligand_info = [hasattr(ligand, 'identical_ligand_info') for ligand in ligands.db.values()]
    # df_ligands['has_identical_ligand_info'] = has_identical_ligand_info
    # df_ligands = df_ligands[['has_identical_ligand_info', 'pred_charge', 'pred_charge_is_confident', 'denticity', 'n_atoms', 'stoichiometry', 'heavy_atoms_graph_hash_with_metal']]

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

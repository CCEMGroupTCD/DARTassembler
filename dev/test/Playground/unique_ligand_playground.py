"""
Playground which reads in the unique ligand db and lets you play with it.
"""
from pathlib import Path
import pandas as pd
import numpy as np
# from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db, load_jsonlines, load_json
from tqdm import tqdm
from copy import deepcopy
from DARTassembler.src.ligand_extraction.DataBase import LigandDB

import ase
from ase.visualize import view
from DARTassembler.src.constants.Paths import default_ligand_db_path, test_ligand_db_path
from DARTassembler.src.ligand_extraction.utilities_graph import get_graph_hash


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

def rotate_ligand_around_donors_inplace(ligand: ase.Atoms, donors_positions: np.array, metal_position: np.array, angle: float=180) -> None:
    """
    Returns the rotated ligand as an ase.Atoms object.
    """
    assert np.allclose(metal_position, (0, 0, 0)), f'Metal position is not at (0, 0, 0), but at {metal_position}. This is covered by the code, but it\'s unexpected and possibly a bug.'
    # The ideal vector to rotate is the average vector from the metal to all donors.
    rot_vector = (donors_positions - metal_position).sum(axis=0) / len(donors_positions)
    orig_ligand_positions = deepcopy(ligand.get_positions())
    ligand.rotate(angle, rot_vector, center=metal_position)
    if angle == 180:    # Check that rotation by further 180 degrees brings the ligand back to its original position
        check_back_rotation_ligand = deepcopy(ligand)
        check_back_rotation_ligand.rotate(angle, rot_vector, center=metal_position)
        assert np.allclose(check_back_rotation_ligand.get_positions(), orig_ligand_positions)

    return None


if __name__ == '__main__':

    n_max = 5000
    denticities = [2, 3]



    # ===== Code for MetaLig documentation =====

    # Load the first 1000 out of 41,018 ligands in the MetaLig database.
    metalig = LigandDB.load_from_json(path=default_ligand_db_path, n_max=n_max)
    metalig = metalig.get_db_with_only_certain_denticities(denticities=denticities)

    # Find symmetrical ligands, first for bidentates.
    # 3D symmetry with SOAP is here outcommented because I want to focus on graph symmetry first.
    # from dscribe.descriptors import SOAP
    # from dscribe.kernels import AverageKernel, REMatchKernel
    # desc = SOAP(species=list(range(1, 100)), r_cut=100.0, n_max=2, l_max=2, sigma=0.2, compression={"mode": "crossover"})
    # re = REMatchKernel()
    run = 'DART'
    data = []
    for uname, ligand in tqdm(metalig.db.items()):
        if ligand.denticity in denticities:
            # Check if ligand graph is symmetrical between donors
            graph, metal_idx = ligand.get_graph_with_metal(metal_symbol='Hg', return_metal_index=True)
            # Make two new graphs, each where one bond connected to the metal is removed
            donor_graphs = []
            for donor_idx in graph.neighbors(metal_idx):
                donor_graph = graph.copy()
                donor_graph.remove_edge(metal_idx, donor_idx)
                donor_graphs.append(donor_graph)
            # Calculate graph hashes of these two graphs. If they are identical, the ligand is symmetrical.
            graph_hashes = [get_graph_hash(donor_graph) for donor_graph in donor_graphs]
            symmetrical = len(set(graph_hashes)) < len(graph_hashes)
            data.append({'uname': uname, 'dent': ligand.denticity, 'formula': ligand.stoichiometry, 'symm_graph': symmetrical})

            ### Detect 3D symmetrical ligands. ###
            ### Quite involved, therefore first I focused on just 2D symmetry. ###
            # atoms_original = ligand.mol
            # atoms_flipped = deepcopy(atoms_original)
            # rotate_ligand_around_donors_inplace(atoms_flipped, donors_positions=ligand.get_donor_positions(), metal_position=ligand.original_metal_position, angle=180)
            # donor_positions = ligand.get_donor_positions()
            # # view(atoms_original)
            # # view(atoms_flipped)
            # # Choose a point that is close to many atoms so that changes in the ligand structure are captured.
            # com = ligand.mol.get_center_of_mass()
            # features1 = desc.create(atoms_original, centers=[com])
            # features2 = desc.create(atoms_flipped, centers=[com])
            # re_kernel = re.create([features1, features2])
            # similarity = re_kernel[0, 1]
            # dissimilarity = (1 - similarity) * 10000
            # data[-1]['diss'] = dissimilarity

    df = pd.DataFrame(data)
    # df = df.sort_values(['symm_graph', 'diss'], ascending=[False, False])
    # Save as concatenated .xyz file
    outpath = Path('/Users/timosommer/Downloads/metalig1000.xyz')
    outpath.unlink(missing_ok=True) # Remove if already exists
    with open(outpath, 'w') as f:
        for name in df['uname'].values:
            lig = metalig.db[name]
            f.write(lig.get_xyz_file_format_string(comment=name, with_metal=True))





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

"""
This script is for benchmarking methods to identify the middle donor atom of a mer-tridentate ligand.
"""
import itertools
from io import StringIO

from DARTassembler.src.ligand_extraction.DataBase import LigandDB
import numpy as np
import ase
import pandas as pd
import networkx as nx
from DARTassembler.src.ligand_extraction.io_custom import save_to_xyz
from DARTassembler.src.ligand_extraction.utilities_Molecule import get_planarity

from oerdev.timo.playground.try_out_individual_problems.new_dart_assembler.new_assembler import align_donor_atoms

def get_sorted_donor_idc_for_square_bipyramidal(atoms, donor_idc) -> list[int]:
    """
    Get the order of the donor atoms of a square bipyramidal ligand so that the top atom is the first one and then come the indices sorted as neighbors in the square.
    :param atoms: ASE atoms object of the ligand
    :param donor_idc: List of donor atom indices
    :return: List of donor atom indices in the correct order
    """
    # Calculate the planarity of each four of the five atoms
    planarities = []
    for top_idx in donor_idc:
        planar_atoms_idc = [i for i in donor_idc if i != top_idx] # Try these as the planar atoms
        planar_atoms_positions = atoms[planar_atoms_idc].get_positions()
        planarity = get_planarity(planar_atoms_positions)
        planarities.append(planarity)
    # Find the top idx where the other atoms are most planar
    top_idx = donor_idc[np.argmax(planarities)]
    planar_atoms_idc = [i for i in donor_idc if i != top_idx]
    # Now sort the planar atoms in the same way as for the mer-tetradentate ligand
    sorted_planar_idc = get_sorted_donor_idc_for_mer_tetradentate(atoms, planar_atoms_idc)
    sorted_donor_idc = [top_idx] + sorted_planar_idc

    return sorted_donor_idc


def get_mer_tridentate_central_donor_idx_by_largest_distance(atoms, donor_idc) -> int:
    """
    Get the index of the central donor atom of a mer-tridentate ligand by the largest distance between the donor atoms.
    """
    donor_positions = atoms[donor_idc].get_positions()
    distances = np.linalg.norm(donor_positions[:, None] - donor_positions, axis=-1)
    np.fill_diagonal(distances, -np.inf)
    largest_distance_idc = list(np.unravel_index(np.argmax(distances), distances.shape))
    pseudo_bidentate_idc = [donor_idc[i] for i in largest_distance_idc]
    central_idx = [i for i in donor_idc if i not in pseudo_bidentate_idc][0]

    return central_idx

def get_mer_tridentate_central_donor_idx_by_graph_closeness(graph, donor_idc) -> int:
    """
    Get the index of the central donor atom of a mer-tridentate ligand by the closeness centrality of the donor atoms.
    """
    closeness_centrality = nx.closeness_centrality(graph)
    central_idx = max(donor_idc, key=lambda idx: closeness_centrality[idx])

    return central_idx

def get_sorted_donor_idc_for_mer_tridentate(atoms, donor_idc=None) -> list[int]:
    """
    Return the indices of the donor atoms of a mer-tridentate ligand so that indices are in the following order: outer atom, central atom, other outer atom
    :param atoms: ASE atoms object of the ligand
    :param donor_idc: List of donor atom indices
    :return: List of donor atom indices in the correct order
    """
    if donor_idc is None:
        donor_idc = list(range(len(atoms)))
    # Define the target vectors so that the central donor atom is the second one
    target_vectors = [[0, 1, 0], [1, 0, 0], [0, -1, 0]]
    all_rssds = []
    for central_idx in donor_idc:
        outer_idc = [i for i in donor_idc if i != central_idx]
        try_idc_order = [outer_idc[0], central_idx, outer_idc[1]]

        # Try to align the donor atoms and get the resulting rssd
        _, rssd = align_donor_atoms(atoms, try_idc_order, target_vectors, return_rssd=True)
        all_rssds.append(rssd)

    # Find the donor indices with the smallest rssd
    central_idx = donor_idc[np.argmin(all_rssds)]
    # Order the donor atoms so that the central atom is the second one
    other_idc = [i for i in donor_idc if i != central_idx]
    sorted_donor_idc = [other_idc[0], central_idx, other_idc[1]]

    return sorted_donor_idc

def get_sorted_donor_idc_for_mer_tetradentate(atoms, donor_idc) -> list[int]:
    """
    Get the order of the donor atoms of a mer-tetradentate ligand so that indices which follow each other are neighbors in the square of the donor atoms. The final order of indices is normalized by sorting each subset of indices which are in opposite corners of each other.
    :param atoms: ASE atoms object of the ligand
    :param donor_idc: List of donor atom indices
    :return: List of donor atom indices in the correct order
    """
    # Define the target vectors so that the central donor atom is the second one
    target_vectors = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    donor_idc = sorted(donor_idc)   # sort for normalization

    # Keep randomly one atom and check which atom is on the opposite side of the square
    corner_idx = donor_idc[0]
    data = {'rssd': [], 'opposite_corner_idx': []}
    for opposite_corner_idx in donor_idc[1:]:
        other_idc = [i for i in donor_idc if i not in [corner_idx, opposite_corner_idx]]
        assert len(other_idc) == 2
        # Sort the indices so that we assume that atoms which are on opposite sides of the square are not neighbors
        try_idc_order = [corner_idx, other_idc[0], opposite_corner_idx, other_idc[1]]
        # Try to align the donor atoms and get the resulting rssd
        _, rssd = align_donor_atoms(atoms, try_idc_order, target_vectors, return_rssd=True)
        data['rssd'].append(rssd)
        data['opposite_corner_idx'].append(opposite_corner_idx)

    # Find the donor indices with the smallest rssd
    best_opposite_idc = data['opposite_corner_idx'][np.argmin(data['rssd'])]
    other_idc = sorted([i for i in donor_idc if i not in [corner_idx, best_opposite_idc]])  # Sort for normalization
    sorted_donor_idc = [corner_idx, other_idc[0], best_opposite_idc, other_idc[1]]

    return sorted_donor_idc


def get_sorted_donor_idc_for_mer_tetradentate_by_largest_distance(atoms, donor_idc) -> list[int]:
    """
    Get the order of the donor atoms of a mer-tetradentate ligand so that indices which follow each other are neighbors in the square of the donor atoms. Use the largest distance between the donor atoms to determine the order.
    :param atoms: ASE atoms object of the ligand
    :param donor_idc: List of donor atom indices
    :return: List of donor atom indices in the correct order
    """
    donor_positions = atoms[donor_idc].get_positions()
    distances = np.linalg.norm(donor_positions[:, None] - donor_positions, axis=-1)
    np.fill_diagonal(distances, -np.inf)
    largest_distance_idc = list(np.unravel_index(np.argmax(distances), distances.shape))
    corner_idc = sorted([donor_idc[i] for i in largest_distance_idc])   # Sort for normalization
    other_idc = sorted([i for i in donor_idc if i not in corner_idc])   # Sort for normalization
    # Normalize so that the lowest index comes first
    if other_idc[0] < corner_idc[0]:
        sorted_donor_idc = [other_idc[0], corner_idc[0], other_idc[1], corner_idc[1]]
    else:
        sorted_donor_idc = [corner_idc[0], other_idc[0], corner_idc[1], other_idc[1]]

    return sorted_donor_idc

if __name__ == '__main__':

    n_max = 1000

    db = LigandDB.load_from_json(n_max=n_max)
    db.db = {name: ligand for name, ligand in db.db.items() if ligand.has_neighboring_coordinating_atoms}
    # Calculate hapticities
    for name, ligand in db.db.items():
        ligand.get_denticities_and_hapticities_idc()

    #
    # data = db.get_ligand_output_df().to_dict(orient='index')
    # for name, ligand in db.db.items():
    #     atoms = ligand.mol
    #     donor_idc = ligand.ligand_to_metal
    #     graph = ligand.get_reindexed_graph()
    #     data[name]['idx_distance'] = get_sorted_donor_idc_for_mer_tetradentate_by_largest_distance(atoms, donor_idc)
    #     data[name]['idx_alignment'] = get_sorted_donor_idc_for_mer_tetradentate(atoms, donor_idc)
    #
    # df = pd.DataFrame.from_dict(data, orient='index')
    # new_cols = ['idx_distance', 'idx_alignment']
    # # Make new columns that are any two of the three columns are the same or different
    # for col1, col2 in itertools.combinations(new_cols, 2):
    #     name1, name2 = col1.removeprefix('idx_'), col2.removeprefix('idx_')
    #     df[f'{name1}={name2}'] = df[col1] == df[col2]
    #     new_cols.append(f'{name1}={name2}')
    #
    # # Shift the new columns to the front
    # df = df[new_cols + [col for col in df.columns if col not in new_cols]]
    #
    # # Print accuracy of distance and graph method compared to alignment method
    # print(f"Accuracy of distance method compared to alignment method: {df['distance=alignment'].mean():.2f}")
    #
    # structures = []
    # for ligand_name, ligand in db.db.items():
    #     top_idx = get_sorted_donor_idc_for_square_bipyramidal(ligand.mol, ligand.ligand_to_metal)[0]
    #     # Change the element of the top atom to Cu
    #     atoms = ligand.get_ase_molecule_with_metal('Ir')
    #     atoms[top_idx].symbol = 'Cu'
    #     structures.append(atoms)
    #
    # # Concatenate the xyz files of the ligands which do not agree
    # ligands = list(db.db.keys())
    # save_to_xyz(outpath='pentadentates.xyz', structures=structures, comments=ligands)

    # todo: Implement a function that finds the effective denticity of ligands and returns a list of lists of donor atom indices.






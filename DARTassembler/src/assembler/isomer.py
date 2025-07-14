"""
This file contains the classes and methods that are used to process the input data and generate the assembled transition metal complex isomers.
"""
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union
from ase.visualize import view
from ase import Atoms
import numpy as np
import itertools
import ase
import networkx as nx
from copy import deepcopy
import logging
from scipy.spatial.distance import cdist, pdist
from scipy.spatial.transform import Rotation as R
from scipy.optimize import linear_sum_assignment, differential_evolution
import pandas as pd

from scipy.optimize import brute

from DARTassembler.src.assembler.utils import are_atoms_equal, assign_ligands_to_vectors
from DARTassembler.src.metalig.geometry import try_all_geometrical_isomer_possibilities, all_geometries, align_vectors, \
    align_donor_atoms
from DARTassembler.src.constants.chem import Element
from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.metalig.mol import BaseMolecule, Ligand
from DARTassembler.src.metalig.utils import check_equal
from DARTassembler.src.metalig.utils_molecule import hapdent_idc_to_donor_idc, get_atomic_props_from_ase_atoms
from DARTassembler.src.misc.io import load_json
from DARTassembler.src.metalig.utils_graph import graph_from_graph_dict

from functools import lru_cache
from logging import getLogger, Logger

# Keep track of 10 different messages and then warn again
@lru_cache(1)
def logging_warn_once(msg: str):
    logging.warning(msg)

class Isomer(BaseMolecule):

    def __init__(self,
                    atomic_props: Union[ase.Atoms, Dict[str, Any]],
                    graph: nx.Graph,
                    metal_idc: List[int],
                    donor_idc: List[List[int]],
                    ligand_idc: List[List[int]],
                    ligand_info: Dict[str, Any] = None,
                    global_props: Dict[str, Any] = None,
                    validity_check: bool = False,
                    warning: str = '',
                    ):
        if global_props is None:
            global_props = {}
        if ligand_info is None:
            ligand_info = {}

        super().__init__(
                         atomic_props=atomic_props,
                         global_props=global_props,
                         graph=graph,
                         validity_check=False   # will be performed later if required
                         )

        self.metal_idc = metal_idc
        self.donor_idc = donor_idc
        self.ligand_idc = ligand_idc
        self.ligand_info = ligand_info
        self.warning = warning
        self.metals = [self.atoms[idx].symbol for idx in self.metal_idc]

        # A little expensive, but small compared to the mono-axial optimization.
        self.ligands = self._get_ligands()  # Ligand() objects for convenient access to ligands. Won't be saved to disk in the DART workflow.

        if validity_check:
            self._tmc_validity_checks()

    def _get_ligands(self, validity_check: bool = True) -> List[Ligand]:
        ligands = []
        for idx, (isomer_donor_idc, isomer_ligand_indices) in enumerate(zip(self.donor_idc, self.ligand_idc)):
            ligand = Ligand(
                atomic_props=self.atoms[isomer_ligand_indices],
                donor_idc=self.ligand_info['donor_idcs'][idx],
                global_props={'geometry': self.ligand_info['geometries'][idx]},
                graph=self.graph.subgraph(isomer_ligand_indices),
                unique_name=self.ligand_info['unique_names'][idx],
                charge=self.ligand_info['charges'],
                hapdent_idc=self.ligand_info['hapdent_idcs'][idx],
                geometric_isomers_hapdent_idc=self.ligand_info['geometric_isomers_hapdent_idcs'][idx],
                validity_check=validity_check,
            )
            ligands.append(ligand)

        return ligands

    def get_metal_symbols(self) -> List[str]:
        """
        Get the symbols of the metal centers in this isomer.
        :return: List of metal symbols
        """
        return self.atoms[self.metal_idc].get_chemical_symbols()

    def _tmc_validity_checks(self) -> None:
        """Some short checks specifically for transition metal complexes."""
        self._check_if_molecule_valid() # Checks basic molecular properties like atomic_props and graph
        # Doublecheck if all the metals are really metals. Don't raise an error in case it's intentional.
        all_metals = all(Element(metal).is_metal for metal in self.metals)
        if not all_metals:
            warnings.warn(f"Any of the metal centers {self.metals} in Isomer() is not a metal. Providing a chemical element as `metal center` that is not a metal is not a problem, this is just to make you aware.")

        return

    def to_dict(self):
        """
        Converts the Isomer object to a dictionary.
        :return: Dictionary representation of the Isomer object
        """
        d = super().to_dict()   # Base class takes care of atomic_props, global_props, and graph
        d.update({
            "metal_idc": self.metal_idc,
            "donor_idc": self.donor_idc,
            "ligand_idc": self.ligand_idc,
            "ligand_info": self.ligand_info,
        })

        return d

    def get_metal_center_atoms(self) -> ase.Atoms:
        """
        Get the metal atoms of the complex.
        :return: ASE Atoms object containing the metal atoms
        """
        atoms = ase.Atoms()
        for metal_idx in self.metal_idc:
            atoms += self.atoms[metal_idx]
        return atoms

    @classmethod
    def from_json(cls, filepath) -> 'Isomer':
        """
        Loads an Isomer object from a .json file.
        :param filepath: Path to the .json file
        :return: Isomer object
        """
        data = load_json(filepath)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Isomer':
        """
        Creates an Isomer object from a dictionary in the correct format.
        :param data: Dictionary containing the Isomer data
        :return: Isomer object
        """
        data['graph'] = graph_from_graph_dict(data['graph'])
        data['charge'] = 0  # todo
        return cls(**data)


class IsomerFactory(object):

    def __init__(
                    self,
                    ligands: List[Ligand],
                    target_vectors: list[list[list[float]]],
                    metal_centers: Union[List[List[Union[str, List[float]]]], str],
                    ligand_origins: List[List[float]] = None,
                    filter_duplicate_isomers: bool = True,
                    filter_clashing_structures: bool = True,
                    filter_clashing_structures_cov_radii_buffer: float = -0.3,
                    check_metal_clashes: bool = True,
                    filter_duplicate_isomers_method: str = "fingerprint",
                    filter_duplicate_isomers_grid_size: int = 9,
                    isomer_comparison_mode: str = "max_diff",
                    isomer_comparison_grouping_mode: str = "cluster",
                    isomer_comparison_grouping_cutoff: float = 1.0,
                    swap_groups: Optional[List[int]] = None,
                    opt_mono_rot: Optional[bool] = True,
                    all_combinations: bool = False,
                    ):
        """
        Generates isomers from a list of ligands, target vectors and metal centers. The orientation of the ligands relative to its metal center is determined by the target vectors. The ligand_origins can be used to shift the ligand with respect to the metal center. If `metal_centers` is a chemical element such as 'Ru', it is assumed to be a mono-metallic complex at the origin.
        :param ligands: List of Ligand objects from the MetaLig database.
        :param target_vectors: List of target vectors for each ligand.
        :param metal_centers: List of tuple with element and position for each metal center. If a string is provided, it is assumed to be the chemical element of a mono-metallic complex at the origin.
        :param ligand_origins: List of the origin for each ligand.

        Example usage for assembling a mono-metallic square-planar Pd complex with 2 cis bidentate ligands in the xy-plane:
            factory = IsomerFactory(
                                    ligands=..., (list of two bidentate Ligand() objects, usually from the MetaLig database)
                                    metal_centers='Pd',
                                    target_vectors=[
                                                    [[1, 0, 0], [0, 1, 0]],     # the donor atoms of the first bidentate ligand are oriented in (+x,+y) direction
                                                    [[-1, 0, 0], [0, -1, 0]],   # the donor atoms of the second bidentate ligand are oriented in (-x,-y) direction
                                                    ],
                                    )
            isomers = factory.generate_isomers()

        Example usage for assembling a bi-metallic complex with three monodentate ligands, one of them bridging:
            ligands = ... (list of three monodentate Ligand() objects, usually from the MetaLig database)
            ru = ['Ru', [1, 0, 0]]
            fe = ['Fe', [-1, 0, 0]]
            metal_centers = [
                                [ru],       # metal center for the first ligand
                                [ru, fe],   # metal centers for the second, bridging ligand
                                [fe]        # metal center for the third ligand
                            ]
            target_vectors = [
                                [[1, 0, 0]],
                                [[0, 0, 1]],
                                [[-1, 0, 0]],
                             ]
            factory = IsomerFactory(
                                        ligands=ligands,
                                        target_vectors=target_vectors,
                                        metal_centers=metal_centers
                                        )
            isomers = factory.generate_isomers()
        """
        metal_centers, ligand_origins, ligands, target_vectors = self._check_input_and_handle_defaults(metal_centers, ligand_origins, ligands, target_vectors)

        self.ligands = ligands
        self.target_vectors = target_vectors
        self.ligand_origins = ligand_origins
        self.metal_centers = metal_centers

        self.filter_duplicate_isomers = filter_duplicate_isomers
        self.filter_clashing_structures = filter_clashing_structures
        self.filter_clashing_structures_cov_radii_buffer = filter_clashing_structures_cov_radii_buffer
        self.check_metal_clashes = check_metal_clashes
        self.filter_duplicate_isomers_method = filter_duplicate_isomers_method
        self.filter_duplicate_isomers_grid_size = filter_duplicate_isomers_grid_size
        self.isomer_comparison_mode = isomer_comparison_mode
        self.isomer_comparison_grouping_mode = isomer_comparison_grouping_mode
        self.isomer_comparison_grouping_cutoff = isomer_comparison_grouping_cutoff
        self.swap_groups = swap_groups
        self.opt_mono_rot = opt_mono_rot
        self.all_combinations = all_combinations

    def _check_input_and_handle_defaults(self, metal_centers, ligand_origins, ligands, target_vectors):
        if isinstance(metal_centers, str):
            # If the metal center is provided as a chemical element, it's a mono-metallic complex at the origin
            metal_centers = [[ase.Atom(symbol=metal_centers, position=[0, 0, 0])] for _ in ligands]
        else:
            # If the metal center is provided as a list of elements and positions, convert to ASE Atoms objects
            metal_centers = [[ase.Atom(symbol=metal[0], position=metal[1]) for metal in metal_list] for metal_list in
                             metal_centers]

        if ligand_origins is None:
            # The `ligand_origins` defaults for ligands that are connected to only one metal to its metal center coordinates, and for ligands that are connected to multiple metal centers to the average of the metal center coordinates (i.e. in the center of them).
            ligand_origins = []
            for lig, centers in zip(ligands, metal_centers):
                # If the ligand is connected to multiple metal centers, use the average position of the metal centers
                avg_position = np.mean([metal.position for metal in centers], axis=0)
                ligand_origins.append(avg_position.tolist())

        # Check input format of lists
        input_lengths = (len(ligands), len(target_vectors), len(ligand_origins), len(metal_centers))
        all_same_length = all(length == input_lengths[0] for length in input_lengths)
        if not all_same_length:
            raise ValueError(f'Input lists must all have the same length. Got lengths: {input_lengths}')

        # Use abbreviations for target vectors
        target_vector_abbreviations = {'x': [1, 0, 0], 'y': [0, 1, 0], 'z': [0, 0, 1], '-x': [-1, 0, 0],
                                       '-y': [0, -1, 0], '-z': [0, 0, -1]}
        try:
            target_vectors = [[target_vector_abbreviations.get(v, v) for v in target_vector_list] for target_vector_list in target_vectors]
        except:
            pass
        # Check target vectors format
        try:
            for target_vector_list in target_vectors:
                try:
                    array = np.array(target_vector_list)
                except ValueError as e:
                    raise ValueError(
                        f"Target vector is not a list of lists of floats: {target_vector_list}. Error: {e}")
                if array.ndim != 2:
                    raise ValueError(
                        f"Target vector must have 2 dimensions (list of list), got {array.shape}: {target_vector_list}")
                elif array.shape[1] != 3:
                    raise ValueError(
                        f"Target vector must have 3 elements (list of list of floats), got {array.shape[1]}: {target_vector_list}")
                elif not np.issubdtype(array.dtype, np.floating) and not np.issubdtype(array.dtype, np.integer):
                    raise ValueError(
                        f"Target vector must be a list of lists of floats, got {array.dtype}: {target_vector_list}")
        except TypeError as e:
            raise TypeError(f"Target vectors must be a list of lists of lists of floats. Error: {e}")

        # Check if the target vectors are compatible with the ligand geometries.
        for target_vector_list, ligand in zip(target_vectors, ligands):
            ligand_geom_vector = np.asarray(all_geometries[ligand.n_eff_denticities][ligand.geometry][0])
            _, rssd = align_vectors(target_vectors=target_vector_list, donor_vectors=ligand_geom_vector)
            if not np.isclose(rssd, 0.0):
                # Check if any combination of the target vectors matches the geometry of the ligand.
                n = len(target_vectors)
                target_vector_list_permutations = list(itertools.permutations(target_vector_list, n))
                any_other_order_matches = False
                for target_vector_list_test in target_vector_list_permutations:
                    _, rssd = align_vectors(target_vectors=target_vector_list_test, donor_vectors=ligand_geom_vector)
                    if np.isclose(rssd, 0.0):
                        any_other_order_matches = True
                        break
                if any_other_order_matches:
                    warn_string = f'Most likely, this is due to a wrong order of the provided target vectors for each donor (see documentation), because if we change the order, they fit perfectly. '
                else:
                    warn_string = f'Most likely, this is due to erroneously provided ligands/target vectors, or simply because the provided target vectors are intended to be different from the ideal geometry of the ligand. '
                logging.warning(f"WARNING: Provided target vectors `{target_vector_list}` do not perfectly match the ligand geometry `{ligand.geometry}`, which has ideal target vectors of {ligand_geom_vector.tolist()}. The assembler will continue with the input you provided, but the assembled complexes may not have the intended geometry. {warn_string}If this is intended, you can ignore this warning.")

        return metal_centers, ligand_origins, ligands, target_vectors

    # This is the method used in the DART workflow to generate isomers.
    def generate_isomers(self) -> List['Isomer']:
        """
        Generates all possible isomers from the ligands and metal centers provided in the constructor.
        :return: List of Isomer objects.
        """
        # Generate all possible geometric isomers to be generated via exchanging ligands (or, as here implemented, exchanging the target vectors of the ligands).
        target_vector_combs = assign_ligands_to_vectors(ligands=self.target_vectors, swap_groups=self.swap_groups)
        ligand_origin_combs = assign_ligands_to_vectors(ligands=self.ligand_origins, swap_groups=self.swap_groups)

        isomers = []
        ligands = self.ligands
        for target_vectors, ligand_origins in zip(target_vector_combs, ligand_origin_combs):
            rotated_ligands = self._get_rotated_ligands(ligands=ligands, target_vectors=target_vectors, ligand_origins=ligand_origins)

            # Generate all combinations. Each combination is a tuple with one isomer per ligand.
            combinations = list(itertools.product(*rotated_ligands))
            unique_metal_centers = self._get_all_unique_metal_centers()
            ase_isomers = []
            for combo in combinations:
                combined = Atoms()  # Start with an empty Atoms object.

                for atom in unique_metal_centers:
                    combined += atom
                for ligand in combo:  # Iterate over the ligands in the combination.
                    combined += ligand  # combining Atoms objects.
                ase_isomers.append(combined)  # Store all the new isomers
            metal_idc = [idx for idx in range(len(unique_metal_centers))]

            for ase_isomer in ase_isomers:
                # Merge the graphs of the ligands and the metal centers to get the full graph of the complex.
                graph, ligand_indices, donor_indices = self._get_merged_graph_from_ligands_and_metal_centers()
                global_props = {}  # Will be populated during the DART workflow.
                # To save disk space, each complex is saved with only the most important information about its ligands. The most important one, the ligand_idc, will be saved in the isomer object anyway.
                ligand_info = {
                    # Important info for making Ligands() objects in the Isomer().
                    'unique_names': [lig.unique_name for lig in ligands],
                    'geometries': [lig.geometry for lig in ligands],
                    'donor_idcs': [lig.donor_idc for lig in ligands],
                    'charges': [lig.charge for lig in ligands],
                    'stoichiometries': [lig.stoichiometry for lig in ligands],
                    'hapdent_idcs': [lig.hapdent_idc for lig in ligands],
                    'geometric_isomers_hapdent_idcs': [lig.geometric_isomers_hapdent_idc for lig in ligands],
                    'target_vectors': target_vectors,
                    # Convenience information for the output csv.
                    'donors': ['-'.join(sorted(lig.donor_elements)) for lig in ligands]
                }
                isomer = Isomer(
                    atomic_props=ase_isomer,
                    graph=graph,
                    metal_idc=metal_idc,
                    ligand_idc=ligand_indices,
                    donor_idc=donor_indices,
                    global_props=global_props,
                    ligand_info=ligand_info,
                    warning='',  # Initially no warning, will be updated later if needed
                    validity_check=True,
                )
                isomer = AxialOptModifier(isomers=[isomer], opt=self.opt_mono_rot).modify(
                    target_vectors=target_vectors, ligand_origins=ligand_origins)[0]
                isomers.append(isomer)

        # Warnings for each isomer. If an isomer has no issues, the note is ''. If an isomer is excluded because of clashing ligands or because it's equivalent to another one, the note is `clashes' or `duplicate`.
        for idx, isomer in enumerate(isomers):
            if isomer.warning == '' and self.filter_clashing_structures:
                clashfilter = IsomerClashFilter(
                    buffer=self.filter_clashing_structures_cov_radii_buffer,
                    check_metal_clashes=self.check_metal_clashes
                )
                clashing = clashfilter.has_clashing_atoms(
                    atoms=isomer.atoms,
                    ligand_idc=isomer.ligand_idc,
                    metal_idc=isomer.metal_idc
                )
                if clashing:
                    isomer.warning = 'clashes'
                    continue

        # Add inplace warnings for isomers that are duplicates of each other. Very important to do this after the clash filter, because otherwise a clashing isomer might be seen as the "first duplicate" and thus be kept, while the other, not-clashing isomer would be discarded.
        not_clashing_isomers = [isomer for isomer in isomers if isomer.warning == '']
        DuplicateIsomerFilter(
                                isomers=not_clashing_isomers,
                                method=self.filter_duplicate_isomers_method,
                                grid_size=self.filter_duplicate_isomers_grid_size,
                                isomer_comparison_mode=self.isomer_comparison_mode,
                                isomer_comparison_grouping_mode=self.isomer_comparison_grouping_mode,
                                fingerprint_grouping_cutoff=self.isomer_comparison_grouping_cutoff,
                                metal_centres=self.metal_centers
                                ).filter()

        # Sort isomers so that the ones without warnings are first, so that the indices of the names are consecutive.
        successful_isomers = [isomer for isomer in isomers if isomer.warning == '']
        unsuccessful_isomers = [isomer for isomer in isomers if isomer.warning != '']
        isomers = successful_isomers + unsuccessful_isomers  # Put successful isomers first, then unsuccessful ones, while keeping the order of isomers within each group.

        return isomers

    def _get_rotated_ligands(self, ligands, target_vectors, ligand_origins) -> list[list[Atoms]]:
        """
        Rotates the ligands according to the target vectors and returns a list of rotated ligands.
        :return: List of lists of ASE Atoms objects, where the outer list corresponds to each ligand and the inner lists correspond to the isomers of that ligand.
        """
        rotated_ligands = []
        for ligand, target_vector_list, origin in zip(ligands, target_vectors, ligand_origins):
            # Extract the geometry and donor atoms of the ligand
            atoms, donor_atoms = ligand.get_isomers_effective_ligand_atoms_with_effective_donor_indices()
            # Cast the target vectors to numpy arrays
            target_vector_list = [np.array(v) for v in target_vector_list]

            # Align the donor atoms of the ligand to the target vectors. Either make all possible geometrical isomers or just the ones specified in the MetaLig. For most cases these two should be identical, but making all combinations has two consequences: (a) the order of input target vectors does not matter (the one with the lowest error is always assembled) and (b) some geometries have more isomers, e.g. the `trigonal` geometry has in theory three isomers in which simply the ligand is rotated, but for the MetaLig we had decided to filter out these isomers so that only one is kept.
            if self.all_combinations:
                ligand_isomers, _, _ = try_all_geometrical_isomer_possibilities(atoms=atoms, donor_idc=donor_atoms[0], target_vectors=target_vector_list)
            else:
                ligand_isomers = [align_donor_atoms(atoms, donor_idc=idc, target_vectors=target_vector_list, return_rssd=False) for idc in donor_atoms]

            # Remove the dummy atom from the haptic ligands
            ligand_isomers = self._remove_haptic_dummy_atom(atoms_list=ligand_isomers, dummy_atom="Cu", donor_atoms_idc=ligand.hapdent_idc)

            # Translate the ligand to its correct location in the complex
            for ligand_isomer in ligand_isomers:
                # Note: This method assumes that the ligand has been pre-translated to 0,0,0, which is the case for all ligands in the MetaLig database.
                ligand_isomer.set_positions(ligand_isomer.get_positions() + np.array(origin))

            # Append the rotated ligands to the list
            rotated_ligands.append(ligand_isomers)

        return rotated_ligands

    @staticmethod
    def _remove_haptic_dummy_atom(atoms_list: List[ase.Atoms], dummy_atom: str, donor_atoms_idc: Tuple[Tuple[int]]):
        """
        Removes the dummy atom from the generated isomers.
        :param atoms_list: List of ASE Atoms objects representing the isomers.
        :param dummy_atom: The symbol of the dummy atom to remove (e.g., "Cu").
        :param donor_atoms_idc: Indices of the donor atoms in the isomers. If the donor atoms are tuples, it indicates haptic coordination.
        :return: List[Atoms]
        """
        # Check to see if there is haptic coordination
        haptic_coordination = False
        for donor_atoms in donor_atoms_idc:
            if type(donor_atoms) == tuple:
                haptic_coordination = True
                break
            else:
                pass

        # If there is no haptic coordination, return the atoms list as is
        if not haptic_coordination:
            return atoms_list

        # If there is haptic coordination, remove the dummy atom from the donor atoms
        else:
            for atoms in atoms_list:
                dummy_idc = [i for i, atom in enumerate(atoms) if
                             atom.symbol == dummy_atom]
                dummy_idc.sort(
                    reverse=True)  # This is important so that the larger index is removed first so as not to change the index of the other atoms
                for dummy_idx in dummy_idc:
                    atoms.pop(dummy_idx)
            return atoms_list

    def _get_all_unique_metal_centers(self) -> List[ase.Atom]:
        """
        Get a list of all unique metal centers.
        :return: List of ase.Atom objects
        """
        unique_metal_centers = [self.metal_centers[0][0]]  # initialize the list with the first metal center
        for metal_list in self.metal_centers:
            for metal in metal_list:
                metal_in_list = any([are_atoms_equal(metal, m) for m in unique_metal_centers])
                if not metal_in_list:
                    unique_metal_centers.append(metal)

        return unique_metal_centers

    def _get_merged_graph_from_ligands_and_metal_centers(self) -> tuple[nx.Graph, list, list]:
        """
        Merges the graphs from the ligands into one graph. The metal is added as a node with index 0 and connected to the donor atoms of the ligands.
        # todo
            This function does not yet work perfectly for multidentate bridging ligands. The donor atoms of the ligands are not correctly connected to the metal center. E.g. for a bidentate bridging atom with two metal donors, all donor atoms are connected to each metal of the two metal centers. This needs to be either fixed or documented.
        :return: Tuple of the merged graph of the complex, the indices of the ligand atoms and the indices of the ligand donor atoms
        """
        ligand_graphs = [deepcopy(lig.graph) for lig in self.ligands]
        unique_metal_centers = self._get_all_unique_metal_centers()

        # Create the new graph by merging everything
        graph = nx.Graph()
        for i, unique_metal_center in enumerate(unique_metal_centers):
            graph.add_nodes_from([(i, {"node_label": unique_metal_center.symbol})])

        # Relabel the nodes of the old graphs so that they are unique for the next step
        i = len(unique_metal_centers)  # start after the metals
        ligand_indices = []
        for ligand_graph in ligand_graphs:
            node_mapping = {node: i + k for k, node in enumerate(sorted(ligand_graph.nodes))}
            nx.relabel_nodes(ligand_graph, mapping=node_mapping, copy=False)
            ligand_indices.append(list(node_mapping.values()))
            i += len(ligand_graph.nodes)

        # Copy the ligand graphs
        for ligand_graph in ligand_graphs:
            graph.add_nodes_from(ligand_graph.nodes(data=True))  # add ligand nodes
            graph.add_edges_from(ligand_graph.edges())  # add ligand edges

        # Connect the metal centers to the ligands
        ligand_donor_indices = [[] for _ in self.ligands]
        for i, (ligand, ligand_metal_centers, ligand_graph) in enumerate(zip(self.ligands, self.metal_centers, ligand_graphs)):
            for metal_center in ligand_metal_centers:
                unique_metal_center_idx = \
                [i for i, atom in enumerate(unique_metal_centers) if are_atoms_equal(atom, metal_center)][0]
                for atomic_donor_idx in ligand.donor_idc:
                    assert ligand.atomic_props['atoms'][
                               atomic_donor_idx] in ligand.donor_elements, f"Atom {ligand.atomic_props['atoms'][atomic_donor_idx]} is not a donor atom of ligand."
                    graph_donor_idx = sorted(ligand_graph.nodes)[atomic_donor_idx]
                    graph.add_edge(unique_metal_center_idx, graph_donor_idx)
                    if graph_donor_idx not in ligand_donor_indices[i]:
                        ligand_donor_indices[i].append(graph_donor_idx)

        # Check if everything is valid
        assert nx.is_connected(graph), "The graph is not fully connected!"
        assert all([set(ligand_donor_indices[i]).issubset(set(ligand_indices[i])) for i in
                    range(len(ligand_indices))]), "The ligand donor indices are not subset of the ligand indices!"
        assert sorted(graph.nodes) == list(
            range(len(graph.nodes))), f"The graphs indices are not in order: {list(graph.nodes)}"

        all_atomic_elements = [unique_metal_center.symbol for unique_metal_center in unique_metal_centers]
        for ligand in self.ligands:
            all_atomic_elements += ligand.atomic_props['atoms']
        all_graph_elements = [graph.nodes[node]['node_label'] for node in sorted(graph.nodes)]
        assert all_graph_elements == all_atomic_elements, f"The graph elements do not match the atomic elements: {all_graph_elements} vs {all_atomic_elements}!"

        atomic_donor_elements = sorted([el for lig in self.ligands for el in lig.donor_elements])
        graph_donor_elements = sorted(
            [graph.nodes[node]['node_label'] for idc in ligand_donor_indices for node in sorted(graph.nodes) if
             node in idc])
        assert atomic_donor_elements == graph_donor_elements, f"The atomic donor elements do not match the graph donor elements: {atomic_donor_elements} vs {graph_donor_elements}!"

        # For debugging: Plot the graph only for the metals and the coordination atoms
        # plot_graph = deepcopy(graph)
        # keep_idc = list(range(len(unique_metal_centers))) + [idx for idc in ligand_donor_indices for idx in idc]
        # for node in list(plot_graph.nodes):
        #     if node not in keep_idc:
        #         plot_graph.remove_node(node)
        # view_graph(plot_graph)

        # Flatten the ligand donor indices
        donor_idc = [idx for idc in ligand_donor_indices for idx in idc]

        return graph, ligand_indices, donor_idc

class AxialOptModifier:
    def __init__(self, isomers: List['Isomer'], opt: bool = True, distance_cutoff: Optional[float] = 4.0, use_cutoff: bool = False):
        """
        This class will take a list of ASE Atoms objects and optimize mono-coordinating ligands around their coordination axis.
        Optionally, a distance_cutoff (Å) can be supplied so the penalty is only evaluated for pairs closer than this threshold.
        """
        self.input_isomers = isomers
        self.opt_command = opt
        # Optional distance cutoff for the objective function (None = use all pairs)
        self.distance_cutoff: Optional[float] = distance_cutoff if use_cutoff else None
        self.output_isomers = []
        logging.debug(f"AxialOpt initialized with {len(self.input_isomers)} Isomer objects.")

    def modify(self, target_vectors, ligand_origins) -> List['Isomer']:
        """
        Optimize the rotation of each mono-coordinating ligand around their respective coordination axis simultaneously
        """
        # If opt_command is False, return the input isomers without optimization
        if not self.opt_command:
            return self.input_isomers

        # Loop through each of the inputted complexes
        for isomer in self.input_isomers:
            atoms = isomer.atoms.copy()

            # Run the global optimizer.
            bounds = [[0, 360] for _ in target_vectors]
            geometries = [ligand.geometry for ligand in isomer.ligands]
            result = differential_evolution(
                self.objective_function, bounds=bounds,
                args=(target_vectors, ligand_origins, atoms.copy(), isomer.ligand_idc, geometries),
                seed=42, # For reproducibility
                maxiter = 10, popsize=5, polish=False
            )

            # Apply the best rotation angles to the atoms.
            best_ligand_angles = list(result.x)
            for angle, axis, origin, idc, ligand in zip(best_ligand_angles, target_vectors, ligand_origins,
                                                   isomer.ligand_idc, isomer.ligands):
                if ligand.geometry not in ['1_monodentate', '2_trans']:
                    continue
                self.rotate(atoms=atoms, vector=axis[0], origin=origin, idc=idc, angle=angle)

            # Update the Isomer() with the new 3D coordinates
            isomer.atoms = atoms
            isomer.atomic_props = get_atomic_props_from_ase_atoms(atoms)
            for ligand, idc in zip(isomer.ligands, isomer.ligand_idc):
                ligand.atoms = atoms[idc]
                ligand.atomic_props = get_atomic_props_from_ase_atoms(ligand.atoms)

            # Append the new Isomer to the output complexes
            self.output_isomers.append(isomer)

        logging.debug(f"Optimized {len(self.output_isomers)} complexes with mono-coordinating ligand rotations.")
        return self.output_isomers

    def objective_function(self, x: np.ndarray, vectors_in: List[np.array], origins_in: List[np.array],
                           TMC_in: ase.Atoms, ligand_idc: list[list[int]], geometries: List[str]) -> float:
        """
        Objective function to optimize the position of the ligands in the TMC complex.
        :param: x:  Array of angles to rotate each ligand around its respective vector.
        :param: vectors_in: List of vectors for each ligand, where each vector is a dictionary with a single key-value pair.
        :param: origins_in: List of origins for each ligand, where each origin is a list of three floats.
        :param: TMC_in: ASE Atoms object representing the transition metal complex (TMC) to be optimized.
        :return: float: The penalty score based on interatomic distances.
        """
        # Generate a copy of the input complex
        TMC_worker = TMC_in.copy()

        for angle, axis, origin, idc, geometry in zip(list(x), vectors_in, origins_in, ligand_idc, geometries):
            axis = np.asarray(axis).squeeze()   # Ensure axis is a 1D vector for monodentate ligands
            if geometry not in ['1_monodentate', '2_trans']:
                continue
            elif geometry == '2_trans':
                # Reduce the 2D axis of two (hopefully) collinear vectors to a single vector to rotate around.
                if not np.allclose(axis[0], -axis[1]):
                    warnings.warn(f"Ligands with geometry '2_trans' have target vectors that are not collinear. This may lead to unexpected results in the rotation of these ligands around their axis to maximize inter-ligand distance.")
                axis = axis[0]
            assert axis.ndim == 1, 'Axis must be a 1D vector for rotation.'
            self.rotate(atoms=TMC_worker, vector=axis, origin=origin, idc=idc, angle=angle)

        # Vectorised distance computation using condensed matrix
        positions = TMC_worker.positions
        d = pdist(positions)        # condensed distance matrix as flat array (each pair once)

        # Apply optional distance cutoff
        if self.distance_cutoff is not None:
            d = d[d <= self.distance_cutoff]

        # Guard against zero distances
        d = d[d > 0.0]

        if d.size == 0:
            return 0.0

        penalty = np.sum(1.0 / d**2)
        return penalty

    @staticmethod
    def rotate(atoms: Atoms, vector: np.array, origin: np.array, idc: List[int], angle: int):
        """
        Rotate the atoms in the Atoms object (only atoms with indices=idc) around the vector by the specified angle.
        Vectorised implementation: rotation is applied to all selected atoms at once.
        """
        # Normalise rotation vector
        vector = np.asarray(vector, dtype=float)
        vector /= np.linalg.norm(vector)
        origin = np.asarray(origin, dtype=float)

        # Create rotation object
        rotation = R.from_rotvec(np.radians(angle) * vector)

        # Vectorised rotation of the chosen indices
        idc_arr = np.asarray(idc, dtype=int)
        rel = atoms.positions[idc_arr] - origin            # translate
        atoms.positions[idc_arr] = rotation.apply(rel) + origin  # rotate & translate back
        return atoms

    def visualize_structures(self):
        """
        Visualize input and output complexes interleaved using ASE's GUI.
        Each input is followed by its corresponding optimized output.
        """
        if not self.output_isomers:
            print("No output complexes found. Run opt_mono_rotation() first.")
            return

        structures_to_view = []

        for i, (in_isomer, out_isomer) in enumerate(zip(self.input_isomers, self.output_isomers)):
            input_copy = in_isomer.atoms.copy()
            input_copy.info["label"] = f"Input {i}"

            output_copy = out_isomer.atoms.copy()
            output_copy.info["label"] = f"Optimized {i}"

            structures_to_view.extend([input_copy, output_copy])  # interleave input/output

        print(f"Launching viewer for {len(structures_to_view)} structures...")
        view(structures_to_view)

class DuplicateIsomerFilter:
    """
    Class to reduce the number of isomers based on alignment or fingerprint similarity.
    """

    def __init__(self, isomers: List['Isomer'],
                 method: str = "fingerprint",
                 grid_size=9,
                 isomer_comparison_mode: str = "max_diff",
                 isomer_comparison_grouping_mode: str = "cluster",  # 'cluster' or 'cutoff'
                 fingerprint_grouping_cutoff: float = 1.0,
                 metal_centres: List[List[float]] = None,
                 energy_heuristic_mode: str = "max",
                 ):
        """
        Initialize the isomer reduction class.
        :param: isomers: The list of ASE Atoms objects representing isomers.
        :param: method: The method to use for reduction, either 'alignment' or 'fingerprint'.
        :param: grid_size: The number of grid points when scanning from 0 to 360 for exact alignment.
        :param: isomer_comparison_mode: The mode for comparing fingerprints, e.g., 'max_diff', etc.
        :param: isomer_comparison_grouping_mode: The mode for grouping fingerprints, either 'cluster' or 'cutoff'.
        :param: fingerprint_grouping_cutoff: The cutoff value for grouping fingerprints when using 'cutoff' mode.
        """
        self.isomers = isomers
        self.method = method.lower()
        self.grid_size = grid_size  # The number of grid points when scanning from 0 to 360 for exact alignment
        self.isomer_comparison_mode = isomer_comparison_mode
        self.isomer_comparison_grouping_mode = isomer_comparison_grouping_mode
        self.fingerprint_grouping_cutoff = fingerprint_grouping_cutoff
        self.metal_centres = metal_centres
        self.diff_matrix = None  # Placeholder for the fingerprint difference matrix
        self.energy_heuristic_mode = energy_heuristic_mode
        self.output_isomers = []
        self.similarity_cutoff_used = None


    def filter(self) -> List['Isomer']:
        """
        Reduce the number of isomers based on the specified method.
        :return: unique isomers as a list of ASE Atoms objects
        """
        if len(self.isomers) <= 1:
            return self.isomers

        if self.method == "alignment":
            self.output_isomers = self._reduce_by_alignment()
        elif self.method == "fingerprint":
            self.output_isomers = self._reduce_by_fingerprint()
        else:
            raise ValueError(f"Fatal Error: Unsupported reduction method '{self.method}. Supported methods are 'alignment' and 'fingerprint'.")

        logging.debug(f"Reduced isomers from {len(self.isomers)} to {len(self.output_isomers)} using method '{self.method}'.")
        return self.output_isomers

    def _reduce_by_alignment(self) -> List['Isomer']:
        """
        Reduce isomers by aligning them and calculating RMSD or another distance metric.
        :return:
        """
        n = len(self.isomers)
        diff_matrix = np.zeros((n, n))
        counter = 0
        total = n * (n - 1) // 2

        for i in range(n):
            isomer1 = self.isomers[i]
            for j in range(i + 1, n):
                isomer2 = self.isomers[j]
                logging.debug(f"Comparing isomers [{i}] and [{j}] ... [{counter}/{total}]")
                score = self.align_isomers(isomer1.atoms, isomer2.atoms)
                diff_matrix[i, j] = score
                diff_matrix[j, i] = score
                counter += 1

        self.diff_matrix = diff_matrix

        group_labels_matrix = self._analyze_similarity(
            self.diff_matrix,
            quantile=0.2,
            method=self.isomer_comparison_grouping_mode,
            cutoff=self.fingerprint_grouping_cutoff
        )
        self.output_isomers = self._assign_duplicate_warnings(group_labels_matrix, self.isomers)
        return self.output_isomers

    def energy_heuristic(self, stat_atoms: ase.Atoms, rot_atoms: ase.Atoms):
        """
        Heuristic to determine the energy of the isomer to be minimized.

        For each element type, pair each atom in stat_atoms with an atom in rot_atoms
        (both having the same number of atoms for that element) such that:
          - If mode="sum": the sum of distances is minimized.
          - If mode="max": the maximum paired distance is returned.

        :param: stat_atoms: ase.Atoms object (stationary isomer)
        :param: rot_atoms: ase.Atoms object (rotated isomer)
        :return: The energy metric (either sum or max of paired distances).
        """
        # Initialize the total energy and overall max distance.
        total_energy = 0.0
        overall_max_distance = 0.0

        # Get chemical symbols from both isomers
        stat_symbols = stat_atoms.get_chemical_symbols()
        rot_symbols = rot_atoms.get_chemical_symbols()

        # Precompute indices for each element.
        stat_index_dict = {element: [i for i, sym in enumerate(stat_symbols) if sym == element]
                           for element in set(stat_symbols)}
        rot_index_dict = {element: [i for i, sym in enumerate(rot_symbols) if sym == element]
                          for element in set(rot_symbols)}

        # Loop over each unique element in the stationary isomer.
        for element, stat_indices in stat_index_dict.items():
            if element not in rot_index_dict:
                raise ValueError(f"Element {element} missing in rotated isomer")
            rot_indices = rot_index_dict[element]
            if len(stat_indices) != len(rot_indices):
                raise ValueError(f"Mismatch in number of '{element}' atoms between isomers: "
                                 f"{len(stat_indices)} vs. {len(rot_indices)}")
            # Get positions for the current element.
            stat_positions = stat_atoms.positions[stat_indices]
            rot_positions = rot_atoms.positions[rot_indices]

            # Compute the distance matrix using C-optimized function.
            distance_matrix = cdist(stat_positions, rot_positions)

            # Solve the assignment problem to pair atoms optimally
            row_ind, col_ind = linear_sum_assignment(distance_matrix)

            if self.energy_heuristic_mode == "sum":
                total_energy += distance_matrix[row_ind, col_ind].sum()
            elif self.energy_heuristic_mode == "max":
                element_max = distance_matrix[row_ind, col_ind].max()
                overall_max_distance = max(overall_max_distance, element_max)
            else:
                raise ValueError("Invalid mode. Use 'sum' or 'max'.")

        return total_energy if self.energy_heuristic_mode == "sum" else overall_max_distance

    def align_isomers(self, stationary_atoms: ase.Atoms, rotated_atoms: ase.Atoms) -> float:
        """
        Using Scipy's global optimizer to find the optimal rotation that minimizes the distance
        between like atoms in the two isomers. When the two isomers are optimally aligned,
        the RMSD or another measure of inter-atomic distance can be calculated.

        This method, with a large enough grid density, is considered to be the ground truth
        for the optimal rotation of the two isomers. Other Scipy global optimizers can be
        easily integrated here for testing purposes.

        :param stationary_atoms: ASE Atoms object to remain fixed.
        :param rotated_atoms: ASE Atoms object to be rotated.
        :return: float — alignment score (as defined by energy_heuristic or objective_function)
        """
        assert hasattr(self, "metal_centres"), ValueError("Fatal Error: metal_centres must be defined before calling align_isomers. ")
        logging.debug(f"Aligning isomers based on {len(self.metal_centres)} metal centre(s).")

        # Here the number of metal centres and the fact that metal centres of different isomers
        # must always be aligned is taken into account to determine if the isomers are similar or not
        if len(self.metal_centres) == 1:
            # There are 3 axes which an isomer can be rotated around to align it with another isomer
            bounds = [[0, 360] for _ in range(3)]
            # The 3 cardinal axes, properly centered on the metal centre
            axes = np.eye(3)  # Standard Cartesian axes (x, y, z)
            logging.debug("Performing 3D brute-force alignment over x, y, z axes.")

        elif len(self.metal_centres) == 2:
            # The isomer can only be rotated around the metal-metal axis to determine if the isomers are similar or not
            bounds = [[0, 360] for _ in range(1)]
            axis_vector = np.array(self.metal_centres[1]) - np.array(self.metal_centres[0])
            axis_vector /= np.linalg.norm(axis_vector)  # Normalize
            axes = [axis_vector]
            logging.debug(f"Performing 1D brute-force alignment around axis: {axis_vector.tolist()}")

        elif len(self.metal_centres) >= 3:
            # 3 metal centres means each isomer is fixed in space and their geometries can be directly compared
            # We simply return the energy heuristic
            logging.debug("Three or more metal centres detected — skipping alignment and using direct heuristic.")
            return self.energy_heuristic(stat_atoms=stationary_atoms, rot_atoms=rotated_atoms)

        else:
            raise ValueError(f"Fatal Error: Unsupported number of metal centres ({len(self.metal_centres)}). ")

        # Perform brute-force global optimization using the configured grid density
        result_angles = brute(
            func=self.objective_function,
            ranges=bounds,
            args=(stationary_atoms, rotated_atoms, axes),
            Ns=self.grid_size
        )

        min_score = self.objective_function(np.array(result_angles), stationary_atoms, rotated_atoms, axes)
        logging.debug(f"Best alignment score after brute-force search: {min_score:.4f}")
        return min_score

    def objective_function(self, x: np.ndarray, atoms1: ase.Atoms, atoms2: ase.Atoms, axes: np.array):
        """
        Objective function to optimize the position of the ligands in the TMC complex.
        :param x:
        :param atoms1:
        :param atoms2:
        :param axes:
        :return:
        """
        # Copy the input atoms to avoid modifying the original
        stationary_isomer = atoms1.copy()
        rotated_isomer = atoms2.copy()

        # Precompute the combined rotation matrix.
        R_total = self.combined_rotation_matrix(x, axes)

        # Apply the rotation to the rotated isomer
        self.apply_combined_rotation(atoms=rotated_isomer,
                                     R_total=R_total,
                                     center=np.array(self.metal_centres[0]))

        # calculate the energy heuristic that will be minimized
        val = self.energy_heuristic(stat_atoms=stationary_isomer,
                                    rot_atoms=rotated_isomer)
        return val

    @staticmethod
    def combined_rotation_matrix(angles, axes):
        """
        Compute the combined rotation matrix from a list of angles and corresponding axes.
        :param: angles: List of angles in degrees.
        :param: axes: List of rotation axes.
        """
        # Start with identity matrix.
        R_total = np.eye(3)
        for angle, axis in zip(angles, axes):
            # Create a rotation for this angle and axis.
            r = R.from_rotvec(np.deg2rad(angle) * np.array(axis))
            # Combine rotations. Note that the order matters!
            R_total = r.as_matrix() @ R_total
        return R_total

    @staticmethod
    def apply_combined_rotation(atoms, R_total, center):
        """
        Apply the combined rotation to the positions of an ASE Atoms object.
        :param: atoms: ASE Atoms object to rotate.
        :param: R_total: Combined rotation matrix.
        :param: center: Center of rotation
        """
        # Shift positions relative to center, apply rotation, then shift back.
        shifted = atoms.positions - center
        atoms.positions = center + (shifted @ R_total.T)

    def _reduce_by_fingerprint(self):
        """
        Reduce isomers using fingerprint-based similarity clustering.
        Select one representative per cluster labeled 'Close'.
        :return: List of unique isomers after fingerprint-based reduction.
        """
        self.diff_matrix = self._compute_fingerprint_matrix(self.isomers)

        # Cluster labels: "Close" or "Far" for each pair
        group_labels_matrix = self._analyze_similarity(self.diff_matrix, quantile=0.2,
                                                       method=self.isomer_comparison_grouping_mode,
                                                       cutoff=self.fingerprint_grouping_cutoff)
        self.output_isomers = self._assign_duplicate_warnings(group_labels_matrix, self.isomers)
        return self.output_isomers

    @staticmethod
    def _assign_duplicate_warnings(group_labels_matrix: np.ndarray, isomers: List['Isomer']):
        """
        Assign 'duplicate' warning to all but the first occurrence in each group
        :param group_labels_matrix: 2D numpy array with labels 'Close' or 'Far'.
        """
        n = len(isomers)
        for i in range(n):
            for j in range(i + 1, n):
                if group_labels_matrix[i, j] == "Close":
                    # The jth isomer is close to the ith one and j>i, so we mark the jth isomer as a duplicate
                    isomers[j].warning = "duplicate"

        return isomers


    def _analyze_similarity(self, matrix: np.ndarray, quantile: float = 0.2, method: str = "cluster", cutoff: Optional[float] = None) -> np.ndarray:
        """
        Analyze pairwise similarity between isomers using either MeanShift clustering or a hard cutoff.
        :param matrix: 2D numpy array representing the fingerprint difference matrix.
        :param quantile: Quantile for estimating bandwidth in MeanShift mode.
        :param method: 'meanshift' or 'cutoff' to define grouping method.
        :param cutoff: Optional float threshold to use when method='cutoff'.
        :return: 2D numpy array of group labels ('Close' or 'Far').
        """
        if method == "cluster":
            from sklearn.cluster import MeanShift, estimate_bandwidth
            # Flatten the matrix values for clustering.
            # Upper triangle indices are used to avoid redundancy.
            triu_indices = np.triu_indices_from(matrix, k=1)
            values = matrix[triu_indices].reshape(-1, 1)
            # bandwidth = float(estimate_bandwidth(values, quantile=quantile, n_samples=len(values)))
            ms = MeanShift(bandwidth=0.5)
            ms.fit(values)
            cluster_labels = ms.labels_
            cluster_centers = ms.cluster_centers_

            # Identify the cluster whose center is closest to zero.
            close_cluster_idx = np.argmin(np.abs(cluster_centers))

            # Extract the true cutoff as the largest value still labeled "Close"
            close_values = values[np.where(cluster_labels == close_cluster_idx)]
            self.similarity_cutoff_used = float(np.max(close_values))  # ensure it's a float, not a 0-d array

            # Reconstruct group label matrix
            full_labels = np.full(matrix.shape, "Far", dtype=object)
            full_labels[triu_indices] = np.where(cluster_labels == close_cluster_idx, "Close", "Far")
            full_labels[(triu_indices[1], triu_indices[0])] = full_labels[triu_indices]  # mirror
            np.fill_diagonal(full_labels, "Far")  # enforce diagonal is "Far"

            return full_labels

        elif method == "cutoff":
            if cutoff is None:
                raise ValueError("Fatal Error: Cutoff value must be provided when using 'cutoff' method.")
            group_labels = np.where(matrix <= cutoff, "Close", "Far")
            self.similarity_cutoff_used = cutoff
            return group_labels

        else:
            raise ValueError(
                f"Fatal Error: Unsupported similarity analysis method '{method}'. "
                f"Supported options: 'cluster', 'cutoff'.")

    def _compute_fingerprint_matrix(self, isomers: list) -> np.ndarray:
        """
        Efficiently compute the upper-triangle of the fingerprint difference matrix.
        This result is then reflected to the lower triangle to create a full symmetric matrix.
        """
        # Generate fingerprints for each isomer
        fingerprints = [self._compute_sorted_distance_fingerprint(isomer) for isomer in isomers]
        # Ensure all fingerprints have the same length
        assert all(len(fp) == len(fingerprints[0]) for fp in fingerprints)
        # Initialize a square matrix to hold the differences
        n = len(fingerprints)
        diff_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):  # Only compute upper triangle
                diff = self._fingerprint_comparison(fingerprints[i], fingerprints[j], mode=self.isomer_comparison_mode)
                diff_matrix[i, j] = diff
                diff_matrix[j, i] = diff  # Reflect to lower triangle

        return diff_matrix

    @staticmethod
    def _compute_sorted_distance_fingerprint(isomer) -> np.ndarray:
        """
        Compute an inter-atomic distance matrix for an isomer and sort the distances under to two conditions:
        1. entries in the matrix (atom-atom distances) are sorted in ascending order.
        2. the rows of the matrix are sorted in lexicographical order of the element symbols.
            i.e. C-C, C-H, H-H, etc. This ensures like inter-atomic distances are compared between isomers.
        The resulting 1D array acts as a fingerprint for the isomer which can be compared to other isomers.

        :param isomer: ASE Atoms object.
        :return: sorted inter-atomic distance fingerprint as a 1D numpy array.
        """
        positions = isomer.atoms.get_positions()  # (N, 3)
        elements = np.asarray(isomer.atoms.get_chemical_symbols())

        # Get pairwise distances and elements in a fully vectorized manner
        dists = pdist(positions)  # length N*(N-1)//2
        i_idx, j_idx = np.triu_indices(len(elements), k=1)  # matching i<j indices
        e_i, e_j = elements[i_idx], elements[j_idx]

        # Sort elements to ensure that (C-H) and (H-C) are treated the same
        swap = e_i > e_j
        first_elem = np.where(swap, e_j, e_i)
        second_elem = np.where(swap, e_i, e_j)

        # Global sort: primary key = elements, secondary key = distance
        order = np.lexsort((dists, second_elem, first_elem))  # last key is primary
        sorted_dists = dists[order]  # sort distances according to the order

        return sorted_dists

    @staticmethod
    def _fingerprint_comparison(fp1: np.ndarray, fp2: np.ndarray, mode: str = "max_diff"):
        """
        Compute the maximum absolute difference between two sorted fingerprint vectors.
        :param fp1: 1D numpy array fingerprint for isomer 1.
        :param fp2: 1D numpy array fingerprint for isomer 2.
        :param mode: Comparison mode: 'max_diff', 'sum_diff', 'mean_diff', or 'rmsd'.
        :return: Maximum absolute difference.
        """
        logging.debug(f"Comparing fingerprints with mode '{mode}'.")
        if mode == "max_diff":
            return np.max(np.abs(fp1 - fp2))
        elif mode == "sum_diff":
            return np.sum(np.abs(fp1 - fp2))
        elif mode == "mean_diff":
            return np.mean(np.abs(fp1 - fp2))
        elif mode == "rmsd":
            return np.sqrt(np.mean((fp1 - fp2) ** 2))
        else:
            raise ValueError(
                f"Fatal Error: Unsupported fingerprint comparison mode '{mode}'. "
                f"Supported modes are 'max_diff', 'sum_diff', 'mean_diff', and 'rmsd'.")

    def plot_fingerprint_difference_matrix(self,
                                           write_svg: bool = True,
                                           plot_plotly: bool = True,
                                           colorscale_min: float = 0.0,
                                           colorscale_mid: float = 0.5,
                                           color_scale_max: float = 1.0,
                                           min_color: Optional[dict] = None,
                                           mid_color: Optional[dict] = None,
                                           max_color: Optional[dict] = None,
                                           cell_label_mode: str = "value",  # "value", "group", or "none"
                                           ) -> None:
        """
        Visualizes the fingerprint difference matrix as a heatmap using Matplotlib and/or Plotly.

        :param: write_svg: If True, saves the Matplotlib figure as an SVG file.
        :param: plot_plotly: If True, creates an interactive Plotly heatmap.
        :param: colorscale_min: Minimum value for the color scale.
        :param: colorscale_mid: Midpoint value for the color scale.
        :param: color_scale_max: Maximum value for the color scale.
        :param: min_color: Dictionary with RGB values for the minimum color.
        :param: mid_color: Dictionary with RGB values for the midpoint color.
        :param: max_color: Dictionary with RGB values for the maximum color.
        :param: cell_label_mode: 'value' (float), 'group' (Close/Far), or 'none'.
        """
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.pyplot as plt

        assert cell_label_mode in {"value", "group", "none"}, ValueError(
            f"Fatal Error: Unsupported cell label mode '{cell_label_mode}'. "
            "Supported modes are 'value', 'group', and 'none'.")

        # Fallback defaults for colors
        min_color = {"r": 238, "g": 100, "b": 97} if min_color is None else min_color
        mid_color = {"r": 255, "g": 255, "b": 255} if mid_color is None else mid_color
        max_color = {"r": 12, "g": 171, "b": 185} if max_color is None else max_color

        if self.diff_matrix is None:
            self.filter()

        labels_list = [f"{i}" for i in range(len(self.isomers))]
        df = pd.DataFrame(self.diff_matrix, index=labels_list, columns=labels_list)

        group_labels_matrix = self._analyze_similarity(
            matrix=self.diff_matrix,
            method=self.isomer_comparison_grouping_mode,
            cutoff=self.fingerprint_grouping_cutoff
        )

        # --- Matplotlib ---
        if write_svg:
            cmap = LinearSegmentedColormap.from_list("custom_rwb", [
                (colorscale_min, tuple(np.array(list(min_color.values())) / 255)),
                (colorscale_mid, tuple(np.array(list(mid_color.values())) / 255)),
                (color_scale_max, tuple(np.array(list(max_color.values())) / 255))
            ])
            fig, ax = plt.subplots(figsize=(8, 8))
            cax = ax.imshow(df.values, cmap=cmap, vmin=0, vmax=1.0)

            def get_contrasting_text_color(value, vmin=0.0, vmax=0.5):
                norm_val = (value - vmin) / (vmax - vmin)
                r, g, b = cmap(norm_val)[:3]
                luminance = 0.299 * r + 0.587 * g + 0.114 * b
                return 'black' if luminance > 0.5 else 'white'

            for i in range(len(df)):
                for j in range(len(df)):
                    label = ""
                    if cell_label_mode == "value":
                        label = f"{df.values[i, j]:.2f}"
                    elif cell_label_mode == "group":
                        label = group_labels_matrix[i, j]

                    if label:
                        color = get_contrasting_text_color(df.values[i, j])
                        ax.text(j, i, label, ha='center', va='center', color=color, fontsize=8)

            ax.set_xticks(np.arange(len(df)))
            ax.set_yticks(np.arange(len(df)))
            ax.set_xticklabels(labels_list, rotation=0, ha='center', fontsize=12)
            ax.set_yticklabels(labels_list, rotation=0, ha='right', fontsize=12)
            ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

            fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04, label="Difference")
            ax.set_title("Comparison Matrix (Matplotlib)")
            plt.tight_layout()
            plt.savefig("Isomer_comparison_matrix.svg", format="svg")
            plt.close()

        # --- Plotly ---
        if plot_plotly:
            color_scale = [
                [colorscale_min, f"rgb({min_color['r']},{min_color['g']},{min_color['b']})"],
                [colorscale_mid, f"rgb({mid_color['r']},{mid_color['g']},{mid_color['b']})"],
                [color_scale_max, f"rgb({max_color['r']},{max_color['g']},{max_color['b']})"]
            ]
            self.launch_interactive_heatmap(df, group_labels_matrix, color_scale, cell_label_mode)

    def launch_interactive_heatmap(self, df: pd.DataFrame, group_labels_matrix: np.ndarray,
                                   color_scale: list, cell_label_mode: str = "value"):
        """
        Launch an interactive Dash app to render the heatmap and enable isomer alignment viewing.
        """
        # Plotly/Dash imports
        import plotly.express as px
        import dash
        app = dash.Dash(__name__)
        fig = px.imshow(
            df,
            color_continuous_scale=color_scale,
            zmin=0, zmax=1.0,
            title="Comparison Matrix (Interactive)",
            labels={"x": "Isomer Index", "y": "Isomer Index", "color": "Difference"}
        )

        # Get the upper triangle indices (excluding the diagonal)
        triu_indices = np.triu_indices_from(df.values, k=1)

        # Extract the values from the upper triangle (excluding diagonal)
        upper_triangle_values = df.values[triu_indices]

        # Create the histogram
        fig_hist = px.histogram(
            x=upper_triangle_values,
            nbins=400,
            title="Histogram of Comparison Values (Upper Triangle Only)",
            labels={'x': 'Difference Value', 'y': 'Count'}
        )

        if self.similarity_cutoff_used is not None:
            fig_hist.add_vline(
                x=self.similarity_cutoff_used,
                line_width=3,
                line_dash="dash",
                line_color="red",
                annotation_text=f"Cutoff = {self.similarity_cutoff_used:.2f}",
                annotation_position="top left"
            )

        annotations = []
        for i in range(len(df)):
            for j in range(len(df)):
                value = df.values[i, j]
                label = ""
                if cell_label_mode == "value":
                    label = f"{value:.2f}"
                elif cell_label_mode == "group":
                    label = group_labels_matrix[i, j]

                if label:
                    # Compute contrast-aware color
                    r, g, b = [int(c) for c in color_scale[min(2, int(3 * value))][1][4:-1].split(',')]
                    luminance = 0.299 * r + 0.587 * g + 0.114 * b
                    font_color = "black" if luminance > 186 else "white"
                    annotations.append(dict(
                        x=j, y=i, text=label,
                        showarrow=False,
                        font=dict(color=font_color, size=10),
                        xanchor="center", yanchor="middle"
                    ))

        fig.update_layout(
            annotations=annotations,
            dragmode="zoom",
            hovermode="closest",
            height=800,
            width=800
        )

        app.layout = dash.html.Div([
            dash.dcc.Graph(id="heatmap", figure=fig),
            dash.dcc.Graph(id="histogram", figure=fig_hist),

        ])

        @app.callback(
            dash.Output("heatmap", "figure"),
            dash.Input("heatmap", "clickData")
        )
        def display_click(clickData):
            if clickData:
                point = clickData["points"][0]
                i_idx = int(point["y"])
                j_idx = int(point["x"])
                print(f"Clicked cell: ({i_idx}, {j_idx}) — triggering isomer viewer.")
                self.view_isomer_alignment(i_idx, j_idx)
            return fig  # Return same figure to avoid layout refresh

        app.run(debug=False, use_reloader=False)

    def view_isomer_alignment(self, index1: int, index2: int, grid_size: int = None) -> None:
        """
        Visualize two isomers (by index) and their optimal alignment using ASE's viewer.
        Frame 0: isomer1
        Frame 1: isomer2
        Frame 2: overlaid view after rotating isomer2 to align with isomer1.

        :param: grid_size: The number of grid points when scanning from 0 to 360 for exact alignment.
        :param: index1: Index of the reference isomer in self.isomers
        :param: index2: Index of the isomer to be aligned and overlaid
        """
        assert 0 <= index1 < len(self.isomers), f"Index1 out of range: {index1}"
        assert 0 <= index2 < len(self.isomers), f"Index2 out of range: {index2}"

        isomer1 = self.isomers[index1].DART_atoms.copy()
        isomer2 = self.isomers[index2].DART_atoms.copy()
        isomer2_aligned = self.isomers[index2].DART_atoms.copy()

        # Determine alignment axes
        if len(self.metal_centres) == 1:
            bounds = [[0, 360] for _ in range(3)]
            axes = np.eye(3)
        elif len(self.metal_centres) == 2:
            bounds = [[0, 360]]
            axis_vector = np.array(self.metal_centres[1]) - np.array(self.metal_centres[0])
            axis_vector /= np.linalg.norm(axis_vector)
            axes = [axis_vector]
        elif len(self.metal_centres) >= 3:
            bounds = None
            axes = None
            print("Three or more metal centres — skipping rotation; showing structures unaligned.")
        else:
            raise ValueError("Fatal Error: Invalid number of metal centres for alignment.")

        # Perform rotation if applicable
        if len(self.metal_centres) < 3:
            result_angles = brute(
                func=self.objective_function,
                ranges=bounds,
                args=(isomer1, isomer2_aligned, axes),
                Ns=grid_size if grid_size else self.grid_size
            )
            R_total = self.combined_rotation_matrix(result_angles, axes)
            self.apply_combined_rotation(atoms=isomer2_aligned, R_total=R_total, center=np.array(self.metal_centres[0]))

        # Create overlaid image: isomer1 + rotated isomer2
        overlaid = isomer1.copy() + isomer2_aligned.copy()

        # Assign tags so colors are distinguishable
        for atom in isomer1:
            atom.tag = 1
        for atom in isomer2:
            atom.tag = 2
        for atom in overlaid:
            atom.tag = 3

        view([isomer1, isomer2, overlaid], viewer="ase")


class IsomerClashFilter:
    """
    Filters out isomers that have clashes between atoms. Optionally considers ligand-metal and metal-metal clashes.
    """
    def __init__(
                    self,
                    buffer: float = -0.3,
                    check_metal_clashes: bool = False
    ):
        """
        Check if there are any clashing atoms in the isomer.
        :param buffer: Distance buffer below sum of covalent radii to allow (negative means atoms can be closer).
        :param check_metal_clashes: Whether to check for clashes with metal atoms, i.e. metal-ligand or metal-metal clashes.
        """
        self.buffer = buffer
        self.check_metal_clashes = check_metal_clashes

    def has_clashing_atoms(
                            self,
                            atoms: ase.Atoms,
                            ligand_idc: list[list[int]],
                            metal_idc: list[int],
                            ) -> bool:
        """
        Check if there are any clashing atoms in the isomer.
        :param atoms: ASE Atoms object of the isomer
        :param ligand_idc: List of lists containing atom indices for each ligand
        :param metal_idc: List of metal atom indices
        :return: True if there are clashing atoms, False otherwise
        """
        n = len(atoms)
        if n <= 1:                                   # nothing to clash
            return False

        pos = atoms.positions                       # (N, 3)
        symbols = atoms.get_chemical_symbols()      # list[str]

        # Get vectorized covalent radii
        radii_map = {s: Element(s).covalent_radius_angstrom for s in set(symbols)}
        radii = np.fromiter((radii_map[s] for s in symbols), dtype=float)

        # Get arrays of pairwise distances and minimum allowed distances
        dists = pdist(pos)
        i_idx, j_idx = np.triu_indices(n, k=1)   # matches pdist order
        min_allowed = radii[i_idx] + radii[j_idx] + self.buffer

        # Build masks to exclude intra-ligand distances and optionally metal distances
        mask = np.ones_like(dists, dtype=bool)
        # (a) skip intra-ligand pairs
        ligand_id_mask = np.full(n, -1, dtype=np.int32)
        for ligand_id, ligand_indices in enumerate(ligand_idc):
            ligand_id_mask[ligand_indices] = ligand_id
        same_ligand = (ligand_id_mask[i_idx] != -1) & (ligand_id_mask[i_idx] == ligand_id_mask[j_idx])
        mask &= ~same_ligand
        # (b) optionally skip any pair that involves a metal atom, i.e. ligand-metal or metal-metal clashes
        metal_idc = np.asarray(metal_idc, dtype=np.int32)
        if not self.check_metal_clashes and metal_idc.size:
            metal_pair = np.isin(i_idx, metal_idc) | np.isin(j_idx, metal_idc)
            mask &= ~metal_pair

        if not np.any(mask):                        # nothing left to check
            return False

        # Check if any distances are below the allowed minimum which means atoms are clashing
        has_clashes = np.any(dists[mask] < min_allowed[mask])

        return has_clashes

"""
This file contains the classes and methods that are used to process the input data and generate the assembled transition metal complex isomers.
"""
import sys
import warnings
from typing import Dict, Any, List, Optional, Tuple, Union, Iterable
from ase.visualize import view
from ase import Atoms, Atom
import numpy as np
import itertools
import ase
import networkx as nx
from copy import deepcopy
from typing import List
import logging
from scipy.spatial.distance import cdist, pdist
from scipy.spatial.transform import Rotation as R
from scipy.optimize import linear_sum_assignment, differential_evolution
import pandas as pd
from scipy.optimize import brute
from DARTassembler.src.assembler.utils import are_atoms_equal, get_list_with_all_possible_swappings, \
    remove_haptic_dummy_atom, get_complex_name, join_duplicate_groups_by_union
from DARTassembler.src.metalig.archetype import try_all_geometrical_isomer_possibilities, all_archetypes, align_vectors, align_donor_atoms
from DARTassembler.src.constants.chem import Element
from DARTassembler.src.metalig.mol import BaseMolecule, Ligand
from DARTassembler.src.metalig.utils_molecule import get_atomic_props_from_ase_atoms
from DARTassembler.src.metalig.utils_graph import  graph_to_dict_with_node_labels, get_graph_hash

try:    # optional imports for visualization. Don't print if not available (since not essential for the main functionality).
    import plotly.graph_objects as go
    import plotly.express as px
    import plotly.io as pio
    import dash
    pio.renderers.default = 'browser'
except ImportError:
    pass


class AssembledIsomer(BaseMolecule):

    def __init__(self,
                 atomic_props: Union[ase.Atoms, Dict[str, Any]],
                 graph: nx.Graph,
                 metal_idc: List[int],
                 donor_idc: List[List[int]],
                 ligand_idc: List[List[int]],
                 ligand_info: Dict[str, Any] = None,
                 global_props: Dict[str, Any] = None,
                 validity_check: bool = False,
                 target_vectors=None,
                 ligand_origins=None,
                 warning: str = '',
                 isomer_name: str = None
                 ):
        if global_props is None:
            global_props = {}
        if ligand_info is None:
            ligand_info = {}

        super().__init__(
            atomic_props=atomic_props,
            global_props=global_props,
            graph=graph,
            validity_check=False  # will be performed later if required
        )

        self.metal_idc = metal_idc
        self.donor_idc = donor_idc
        self.ligand_idc = ligand_idc
        self.ligand_info = ligand_info
        self.warning = warning
        self.target_vectors = target_vectors
        self.ligand_origins = ligand_origins
        self.isomer_name = isomer_name
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
                global_props={'archetype': self.ligand_info['archetypes'][idx]},
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
        self._check_if_molecule_valid()  # Checks basic molecular properties like atomic_props and graph
        # Doublecheck if all the metals are really metals. Don't raise an error in case it's intentional.
        all_metals = all(Element(metal).is_metal for metal in self.metals)
        if not all_metals:
            warnings.warn(f"Any of the metal centers {self.metals} in AssembledIsomer() is not a metal. Providing a chemical element as `metal center` that is not a metal is not a problem for DART, this warning is just to make you aware.")

        return

    def to_dict(self):
        """
        Converts the AssembledIsomer object to a dictionary.
        :return: Dictionary representation of the AssembledIsomer object
        """
        d = super().to_dict()  # Base class takes care of atomic_props, global_props, and graph
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

    # todo: update these methods using the new output format of the AssembledComplex.to_dict() method.
    # @classmethod
    # def from_json(cls, filepath) -> 'AssembledIsomer':
    #     """
    #     Loads an AssembledIsomer object from a .json file.
    #     :param filepath: Path to the .json file
    #     :return: AssembledIsomer object
    #     """
    #     data = load_json(filepath)
    #     return cls.from_dict(data)
    #
    # @classmethod
    # def from_dict(cls, data: Dict[str, Any]) -> 'AssembledIsomer':
    #     """
    #     Creates an AssembledIsomer object from a dictionary in the correct format.
    #     :param data: Dictionary containing the AssembledIsomer data
    #     :return: AssembledIsomer object
    #     """
    #     data['graph'] = graph_from_graph_dict(data['graph'])
    #     data['charge'] = 0  # todo
    #     return cls(**data)


class AssembledComplex(object):

    def __init__(
            self,
            ligands: List[Ligand],
            target_vectors: list[list[list[float]]],
            metal_centers: Union[List[List[Union[str, List[float]]]], str],
            ligand_origins: List[List[float]] = None,
    ):
        """
        Generates isomers from a list of ligands, target vectors and metal centers. The orientation of the ligands relative to its metal center is determined by the target vectors. The ligand_origins can be used to shift the ligand with respect to the metal center. If `metal_centers` is a chemical element such as 'Ru', it is assumed to be a mono-metallic complex at the origin.
        :param ligands: List of Ligand objects from the MetaLig database.
        :param target_vectors: List of target vectors for each ligand.
        :param metal_centers: List of tuple with element and position for each metal center. If a string is provided, it is assumed to be the chemical element of a mono-metallic complex at the origin.
        :param ligand_origins: List of the origin for each ligand.

        Example usage for assembling a mono-metallic square-planar Pd complex with 2 cis bidentate ligands in the xy-plane:
            factory = AssembledComplex(
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
            factory = AssembledComplex(
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

        # Check if the target vectors are compatible with the ligand archetypes.
        for target_vector_list, ligand in zip(target_vectors, ligands):
            ligand_arch_vector = np.asarray(all_archetypes[ligand.n_eff_denticities][ligand.archetype][0])
            _, rssd = align_vectors(target_vectors=target_vector_list, donor_vectors=ligand_arch_vector)
            if not np.isclose(rssd, 0.0):
                # Check if any combination of the target vectors matches the archetype of the ligand.
                n = len(target_vectors)
                target_vector_list_permutations = list(itertools.permutations(target_vector_list, n))
                any_other_order_matches = False
                for target_vector_list_test in target_vector_list_permutations:
                    _, rssd = align_vectors(target_vectors=target_vector_list_test, donor_vectors=ligand_arch_vector)
                    if np.isclose(rssd, 0.0):
                        any_other_order_matches = True
                        break
                if any_other_order_matches:
                    warn_string = f'Most likely, this is due to a wrong order of the provided target vectors for each donor (see documentation), because if we change the order, they fit perfectly. '
                else:
                    warn_string = f'Most likely, this is due to erroneously provided ligands/target vectors, or simply because the provided target vectors are intended to be different from the ideal archetype of the ligand. '
                logging.warning(
                    f"WARNING: Provided target vectors `{target_vector_list}` do not perfectly match the ligand archetype `{ligand.archetype}`, which has ideal target vectors of {ligand_arch_vector.tolist()}. The assembler will continue with the input you provided, but the assembled complexes may not have the intended archetype. {warn_string}If this is intended, you can ignore this warning.")

        return metal_centers, ligand_origins, ligands, target_vectors

    # This is the method used in the DART workflow to generate isomers.
    def generate_isomers(self,
                         check_duplicate: bool = True,
                         check_clashing: bool = True,
                         swap_groups: Optional[List[int]] = None,
                         optimize_monoaxial: Optional[bool] = True,
                         force_all_isomers: bool = False,
                         clashing_buffer: float = -0.3,
                         clashing_metal: bool = False,
                         duplicate_cutoff: float = 0.5,
                         complex_name_length: int = 8,
                         complex_name_suffix: str = '',
                         avoid_names: Optional[Iterable[str]] = None,
                         ):
        """
        Generates all possible isomers from the ligands and metal centers provided.
        :return: List of AssembledIsomer objects.
        """
        if avoid_names is None:
            avoid_names = set()
        self.check_duplicate = check_duplicate
        self.check_clashing = check_clashing
        self.clashing_buffer = clashing_buffer
        self.clashing_metal = clashing_metal
        self.duplicate_cutoff = duplicate_cutoff
        self.swap_groups = swap_groups
        self.optimize_monoaxial = optimize_monoaxial
        self.force_all_isomers = force_all_isomers
        self.complex_name_length = complex_name_length
        self.complex_name_suffix = complex_name_suffix
        self.avoid_names = avoid_names

        unique_metal_centers = self._get_all_unique_metal_centers()
        self.metal_idc = [idx for idx in range(len(unique_metal_centers))]
        self.graph, self.ligand_indices, self.donor_indices = self._get_merged_graph_from_ligands_and_metal_centers()
        self.graph_hash = get_graph_hash(self.graph)
        self.complex_name = self._get_complex_name(avoid_names=avoid_names)

        # Generate all possible geometric isomers to be generated via exchanging ligands (or, as here implemented, exchanging the target vectors of the ligands).
        target_vector_combs = get_list_with_all_possible_swappings(objects=self.target_vectors, swap_groups=self.swap_groups)
        ligand_origin_combs = get_list_with_all_possible_swappings(objects=self.ligand_origins, swap_groups=self.swap_groups)

        isomers = []
        same_length_target_vectors = []
        same_length_ligand_origins = []
        isomer_idx = 1
        for target_vectors, ligand_origins in zip(target_vector_combs, ligand_origin_combs):
            rotated_ligands = self._get_rotated_ligands(target_vectors=target_vectors, ligand_origins=ligand_origins)

            # Generate all combinations. Each combination is a tuple with one isomer per ligand.
            combinations = list(itertools.product(*rotated_ligands))
            ase_isomers = []
            for combo in combinations:
                combined = Atoms()  # Start with an empty Atoms object.
                for atom in unique_metal_centers:
                    combined += atom
                for ligand in combo:  # Iterate over the ligands in the combination.
                    combined += ligand  # combining Atoms objects.
                ase_isomers.append(combined)  # Store all the new isomers

            for ase_isomer in ase_isomers:
                isomer = AssembledIsomer(
                    atomic_props=ase_isomer,
                    graph=self.graph,
                    metal_idc=self.metal_idc,
                    ligand_idc=self.ligand_indices,
                    donor_idc=self.donor_indices,
                    global_props={},
                    ligand_info=self._get_ligandinfo(),
                    target_vectors=target_vectors,
                    ligand_origins=ligand_origins,
                    warning='',  # Initially no warning, will be updated later if needed
                    validity_check=True,
                    isomer_name=self.complex_name + str(isomer_idx)  # Assign a name based on the complex name and index
                )
                isomers.append(isomer)
                same_length_ligand_origins.append(ligand_origins)
                same_length_target_vectors.append(target_vectors)
                isomer_idx += 1

        pre_isomers_duplicate_groups = DuplicateIsomerFilter(isomers=isomers, fingerprint_grouping_cutoff=self.duplicate_cutoff, metal_centers=self.metal_centers).get_duplicate_groups()
        pre_isomers_duplicate_group_names = [set([isomer.isomer_name for isomer in isomer_group]) for isomer_group in pre_isomers_duplicate_groups]

        # Do a mono-axial optimization of the isomers and afterward check for clashing ligands.
        for idx, isomer, target_vectors, ligand_origins in zip(range(len(isomers)), isomers, same_length_target_vectors, same_length_ligand_origins):
            isomer = AxialOptModifier(isomers=[isomer], opt=self.optimize_monoaxial).modify(target_vectors_list=[target_vectors], ligand_origins_list=[ligand_origins])[0]
            isomers[idx] = isomer # important: copy changed object over to list

            if isomer.warning == '' and self.check_clashing:
                clashfilter = IsomerClashFilter(buffer=self.clashing_buffer, check_metal_clashes=self.clashing_metal)
                clashing = clashfilter.has_clashing_atoms(atoms=isomer.atoms, ligand_idc=isomer.ligand_idc, metal_idc=isomer.metal_idc)
                if clashing:
                    isomer.warning = 'clashing'  # this now updates the optimized object

        # Check for duplicates again after the mono-axial optimization.
        joined_isomers_duplicate_group_names = self.get_duplicate_isomers_group_names(isomers=isomers, pre_isomers_duplicate_group_names=pre_isomers_duplicate_group_names, duplicate_cutoff=self.duplicate_cutoff,  metal_centers=self.metal_centers)

        self.successful_isomers, self.unsuccessful_isomers = self.divide_into_successful_and_unsuccessful_isomers(isomers, joined_isomers_duplicate_group_names)
        self.success = len(self.successful_isomers) > 0
        self.isomers = isomers

        return

    @staticmethod
    def get_duplicate_isomers_group_names(isomers: list[Any], pre_isomers_duplicate_group_names: list[set[str | None]], duplicate_cutoff: float, metal_centers: list[list[Atom]]) -> list[list[str]]:
        """
        Get the joined duplicate isomer group names from pre- and post-monoaxial optimization duplicate groups. Sort the joined groups by the order of `isomers` to preserve the output order.
        :param isomers: List of isomer objects, each having an 'isomer_name' attribute.
        :param pre_isomers_duplicate_group_names: List of sets of isomer names that are considered duplicates before monoaxial optimization.
        :param duplicate_cutoff: Cutoff for duplicate isomer filtering.
        :param metal_centers: List of metal center atoms.
        :return: List of lists of isomer names that are considered duplicates.
        """
        post_isomers_duplicate_groups = DuplicateIsomerFilter(isomers=isomers, fingerprint_grouping_cutoff=duplicate_cutoff, metal_centers=metal_centers).get_duplicate_groups()
        post_isomers_duplicate_group_names = [set([isomer.isomer_name for isomer in isomer_group]) for isomer_group in post_isomers_duplicate_groups]

        # Join the pre- and post-isomers duplicate groups. If an isomer is a duplicate in either the pre- or post-isomers duplicate groups, it is considered a duplicate.
        joined_isomers_duplicate_group_names = join_duplicate_groups_by_union(pre_isomers_duplicate_group_names, post_isomers_duplicate_group_names, mode='post')
        assert sorted(name for group in joined_isomers_duplicate_group_names for name in group) == sorted(name for isomer in isomers for name in [isomer.isomer_name]), "Joined isomer groups do not contain all isomers."

        # Sort the joint isomer names by the order of `isomers` and convert to lists, so that the output order of isomers is preserved. That is particularly important so that the duplicate filter always keeps the same, "first" isomer in the group.
        isomer_names_order = [isomer.isomer_name for isomer in isomers]
        joined_isomers_duplicate_group_names = [sorted(list(group), key=lambda x: isomer_names_order.index(x)) for group
                                                in joined_isomers_duplicate_group_names]
        return joined_isomers_duplicate_group_names

    def divide_into_successful_and_unsuccessful_isomers(self, isomers: list[Any], joined_isomers_duplicate_group_names: list[list[str]]) -> tuple[list[Any], list[Any]]:
        """
        Divide the isomers into successful and unsuccessful isomers. If an isomer has a 'clashing' warning, or if it is a duplicate, it is considered unsuccessful. For duplicates, we use the joint duplicate groups and let through only the first isomer in each group, marking the others as duplicates. It is important that this is done after the clash filter, so that clashing isomers are not considered duplicates.
        :param isomers: List of isomer objects, each having an 'isomer_name' attribute and a 'warning' attribute.
        :param joined_isomers_duplicate_group_names: List of lists of isomer names that are considered duplicates.
        :return: A tuple containing two lists: (successful_isomers, unsuccessful_isomers).
        """
        successful_isomers = []
        unsuccessful_isomers = []
        isomer_dict = {isomer.isomer_name: isomer for isomer in isomers}
        for isomer_group in joined_isomers_duplicate_group_names:
            assert len(isomer_group) == len(set(isomer_group)) and len(
                isomer_group) > 0, f"Duplicate isomer group contains duplicates or is empty: {isomer_group}"
            first_isomer_in_group_added = False
            for isomer_name in isomer_group:
                if isomer_dict[isomer_name].warning == 'clashing':
                    # If the isomer is clashing, add it to the unsuccessful isomers.
                    unsuccessful_isomers.append(isomer_dict[isomer_name])
                elif isomer_dict[isomer_name].warning == '':
                    if not first_isomer_in_group_added:
                        # If the isomer is not clashing and is the first in the group, add it to the successful isomers.
                        successful_isomers.append(isomer_dict[isomer_name])
                        first_isomer_in_group_added = True
                    else:
                        # If the isomer is not clashing and is not the first in the group, add it to the unsuccessful isomers and mark it as a duplicate.
                        duplicate_indices = ','.join(name.removeprefix(self.complex_name) for name in isomer_group)
                        isomer_dict[isomer_name].warning = f'duplicate({duplicate_indices})'
                        unsuccessful_isomers.append(isomer_dict[isomer_name])
                else:
                    raise ValueError(
                        f"Isomer {isomer_name} has an unexpected warning: {isomer_dict[isomer_name].warning}. Expected 'clashing' or ''.")
        return successful_isomers, unsuccessful_isomers


    def to_dict(self):
        """
        Converts the AssembledComplex object to a dictionary.
        :return: Dictionary representation of the AssembledComplex object
        """
        isomer_data = {}
        for isomer in self.isomers:
            isomer_data[isomer.isomer_name] = {
                'atomic_props': isomer.atomic_props,
                'warning': isomer.warning,
                'target_vectors': isomer.target_vectors,
                'ligand_origins': isomer.ligand_origins,
            }
        return {
            "complex_name": self.complex_name,
            "isomers": isomer_data,
            "graph": graph_to_dict_with_node_labels(self.graph),
            "graph_hash": self.graph_hash,
            "metal_idc": self.metal_idc,
            "donor_idc": self.donor_indices,
            "ligand_idc": self.ligand_indices,
            "ligand_info": self._get_ligandinfo(),
            "input": {
                "check_duplicate": self.check_duplicate,
                "check_clashing": self.check_clashing,
                "clashing_buffer": self.clashing_buffer,
                "clashing_metal": self.clashing_metal,
                "duplicate_cutoff": self.duplicate_cutoff,
                "swap_groups": self.swap_groups,
                "optimize_monoaxial": self.optimize_monoaxial,
                "force_all_isomers": self.force_all_isomers,
                "complex_name_length": self.complex_name_length,
                "complex_name_suffix": self.complex_name_suffix,
            }

        }

    def _get_complex_name(self, avoid_names: Optional[Iterable[str]]) -> str:
        return get_complex_name(seed=self.graph_hash, length=self.complex_name_length, suffix=self.complex_name_suffix, avoid_names=avoid_names)

    def _get_ligandinfo(self) -> Dict[str, Any]:
        return {
            # Important info for making Ligands() objects in the AssembledIsomer().
            'unique_names': [lig.unique_name for lig in self.ligands],
            'archetypes': [lig.archetype for lig in self.ligands],
            'donor_idcs': [lig.donor_idc for lig in self.ligands],
            'charges': [lig.charge for lig in self.ligands],
            'stoichiometries': [lig.stoichiometry for lig in self.ligands],
            'hapdent_idcs': [lig.hapdent_idc for lig in self.ligands],
            'geometric_isomers_hapdent_idcs': [lig.geometric_isomers_hapdent_idc for lig in self.ligands],
            # Convenience information for the output csv.
            'donors': ['-'.join(sorted(lig.donor_elements)) for lig in self.ligands]
        }

    def _get_rotated_ligands(self, target_vectors, ligand_origins) -> list[list[Atoms]]:
        """
        Rotates the ligands according to the target vectors and returns a list of rotated ligands.
        :return: List of lists of ASE Atoms objects, where the outer list corresponds to each ligand and the inner lists correspond to the isomers of that ligand.
        """
        rotated_ligands = []
        for ligand, target_vector_list, origin in zip(self.ligands, target_vectors, ligand_origins):
            # Extract the archetype and donor atoms of the effective ligand, potentially with 'Cu' dummy atoms for haptic ligands.
            atoms, donor_atoms = ligand.get_isomers_effective_ligand_atoms_with_effective_donor_indices(dummy='Cu')
            # Cast the target vectors to numpy arrays
            target_vector_list = [np.array(v) for v in target_vector_list]

            # Align the donor atoms of the ligand to the target vectors. Either make all possible geometrical isomers or just the ones specified in the MetaLig. For most cases these two should be identical, but making all combinations has two consequences: (a) the order of input target vectors does not matter (the one with the lowest error is always assembled) and (b) some archetypes have more isomers, e.g. the `trigonal` archetype has in theory three isomers in which simply the ligand is rotated, but for the MetaLig we had decided to filter out these isomers so that only one is kept.
            if self.force_all_isomers:
                ligand_isomers, _, _ = try_all_geometrical_isomer_possibilities(atoms=atoms, donor_idc=donor_atoms[0], target_vectors=target_vector_list)
            else:
                ligand_isomers = [align_donor_atoms(atoms, donor_idc=idc, target_vectors=target_vector_list, return_rssd=False) for idc in donor_atoms]

            # Remove the dummy atom from the haptic ligands
            if ligand.n_haptic_atoms > 0:
                ligand_isomers = [remove_haptic_dummy_atom(atoms=atoms, dummy_atom='Cu') for atoms in ligand_isomers]

            # Translate the ligand to its correct location in the complex
            for ligand_isomer in ligand_isomers:
                # Note: This method assumes that the ligand has been pre-translated to 0,0,0, which is the case for all ligands in the MetaLig database.
                ligand_isomer.set_positions(ligand_isomer.get_positions() + np.array(origin))

            # Append the rotated ligands to the list
            rotated_ligands.append(ligand_isomers)

        return rotated_ligands

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
    def __init__(self, isomers: List['AssembledIsomer'], opt: bool = True, distance_cutoff: Optional[float] = 4.0, use_cutoff: bool = False):
        """
        This class will take a list of ASE Atoms objects and optimize mono-coordinating ligands around their coordination axis.
        Optionally, a distance_cutoff (Å) can be supplied so the penalty is only evaluated for pairs closer than this threshold.
        """
        self.input_isomers = isomers
        self.opt_command = opt
        # Optional distance cutoff for the objective function (None = use all pairs)
        self.distance_cutoff: Optional[float] = distance_cutoff if use_cutoff else None
        self.output_isomers = []
        logging.debug(f"AxialOpt initialized with {len(self.input_isomers)} AssembledIsomer objects.")

    def modify(self, target_vectors_list, ligand_origins_list, maxiter=1000, popsize=15) -> List['AssembledIsomer']:
        """
        Optimize each isomer independently, with its own target_vectors and ligand_origins.
        """
        if not self.opt_command:
            return self.input_isomers

        # Clear output isomers each run
        self.output_isomers = []

        # Sanity check lengths
        if len(self.input_isomers) != len(target_vectors_list) or len(self.input_isomers) != len(ligand_origins_list):
            raise ValueError("Each isomer must have its own set of target_vectors and ligand_origins.")

        # Run the optimizer.
        # Optimize each isomer independently
        for isomer, target_vectors, ligand_origins in zip(self.input_isomers, target_vectors_list, ligand_origins_list):
            atoms = isomer.atoms.copy()

            # Each ligand rotation angle gets its own bound
            bounds = [[0, 360] for _ in target_vectors]

            archetypes = [ligand.archetype for ligand in isomer.ligands]

            # Run the optimizer
            result = differential_evolution(
                self.objective_function,
                bounds=bounds,
                args=(target_vectors, ligand_origins, atoms.copy(), isomer.ligand_idc, archetypes),
                seed=42,
                maxiter=maxiter,
                popsize=popsize,
                polish=True
            )

            best_ligand_angles = list(result.x)

            # Correctly apply rotations to this isomer's ligands
            for angle, axis, origin, idc, ligand in zip(best_ligand_angles, target_vectors, ligand_origins, isomer.ligand_idc, isomer.ligands):
                if ligand.archetype not in ['1-mono', '2-trans']:
                    continue
                self.rotate(atoms=atoms, vector=axis, origin=origin, idc=idc, angle=angle)

            # Copy isomer before modification to avoid unintended side effects
            new_isomer = deepcopy(isomer)
            new_isomer.atoms = atoms
            new_isomer.atomic_props = get_atomic_props_from_ase_atoms(atoms)

            for ligand, idc in zip(new_isomer.ligands, new_isomer.ligand_idc):
                ligand.atoms = deepcopy(atoms[idc])
                ligand.atomic_props = get_atomic_props_from_ase_atoms(ligand.atoms)

            self.output_isomers.append(new_isomer)

        logging.debug(f"Optimized {len(self.output_isomers)} complexes correctly.")
        return self.output_isomers

    def objective_function(self, x: np.ndarray, vectors_in: List[np.array], origins_in: List[np.array],
                           TMC_in: ase.Atoms, ligand_idc: list[list[int]], archetypes: List[str]) -> float:
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

        for angle, axis, origin, idc, archetype in zip(list(x), vectors_in, origins_in, ligand_idc, archetypes):
            if archetype not in ['1-mono', '2-trans']:
                continue
            self.rotate(atoms=TMC_worker, vector=axis, origin=origin, idc=idc, angle=angle)

        # Vectorised distance computation using condensed matrix
        positions = TMC_worker.positions
        d = pdist(positions)  # condensed distance matrix as flat array (each pair once)

        # Apply optional distance cutoff
        if self.distance_cutoff is not None:
            d = d[d <= self.distance_cutoff]

        # Guard against zero distances
        d = d[d > 0.0]

        if d.size == 0:
            return 0.0

        penalty = np.sum(1.0 / d ** 2)
        return penalty

    @staticmethod
    def rotate(atoms: Atoms, vector: np.array, origin: np.array, idc: List[int], angle: int):
        """
        Rotate selected atoms around the given vector by the specified angle.
        Robust to common axis shapes for ligand archetypes:
          - (3,)   : standard single axis vector
          - (1,3)  : single axis wrapped in an extra list
          - (2,3)  : 2-trans case (two anti-parallel donor vectors) -> collapse to the first one
        """
        # Normalize vector shape and cast the list to a numpy array
        vector = np.asarray(vector, dtype=float)

        # If vector is 2D, reduce to 1D
        if vector.ndim == 2:
            if vector.shape == (2, 3): # If 2-trans archetype we should have two near opposite donor atom vectors
                # Check if donor atom vectors are opposite
                # if not np.allclose(vector[0], -vector[1], atol=1e-3):
                #     logging.warning("rotate(): 2-trans target vectors not perfectly opposite; using first vector.")
                vector = vector[0]  # warn if not opposite but still use the first vector
            elif vector.shape == (1, 3):
                vector = vector[0]
            else:
                raise ValueError(f"rotate(): Unexpected axis shape {vector.shape}. Expected (3,), (1,3), or (2,3).")

        elif vector.ndim != 1 or vector.shape[0] != 3:
            raise ValueError(f"rotate(): Axis must be a single 3D vector, got shape {vector.shape}")

        # --- Normalize direction (unit vector) ---
        norm = np.linalg.norm(vector)
        if norm == 0:
            raise ValueError("rotate(): Rotation axis has zero length.")
        vector = vector / norm

        origin = np.asarray(origin, dtype=float)

        # --- Perform rotation ---
        rotation = R.from_rotvec(np.radians(angle) * vector)
        idc_arr = np.asarray(idc, dtype=int)
        rel = atoms.positions[idc_arr] - origin  # translate
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

    def __init__(self, isomers: List['AssembledIsomer'],
                 method: str = "distances",
                 grid_size=9,
                 isomer_comparison_mode: str = "max_diff",
                 isomer_comparison_grouping_mode: str = "cutoff",  # 'cluster' or 'cutoff'
                 fingerprint_grouping_cutoff: float = 0.5,
                 metal_centers: List[List[ase.Atom]] = None,
                 energy_heuristic_mode: str = "max",
                 ):
        """
        Initialize the isomer reduction class.
        :param: isomers: The list of ASE Atoms objects representing isomers.
        :param: method: The method to use for reduction, either 'alignment' or 'distances'.
        :param: grid_size: The number of grid points when scanning from 0 to 360 for exact alignment.
        :param: duplicate_distances_metric: The mode for comparing fingerprints, e.g., 'max_diff', etc.
        :param: duplicate_distances_classifier: The mode for grouping fingerprints, either 'cluster' or 'cutoff'.
        :param: fingerprint_grouping_cutoff: The cutoff value for grouping fingerprints when using 'cutoff' mode.
        """
        if metal_centers is None:
            metal_centers = []
        self.isomers = isomers
        self.method = method.lower()
        self.grid_size = grid_size  # The number of grid points when scanning from 0 to 360 for exact alignment
        self.isomer_comparison_mode = isomer_comparison_mode
        self.isomer_comparison_grouping_mode = isomer_comparison_grouping_mode
        self.fingerprint_grouping_cutoff = fingerprint_grouping_cutoff
        self.metal_centres = metal_centers
        self.unique_metal_centers = list({(atom.symbol, tuple(atom.position)): atom for sublist in self.metal_centres for atom in sublist}.values())
        self.diff_matrix = None  # Placeholder for the fingerprint difference matrix
        self.energy_heuristic_mode = energy_heuristic_mode
        self.isomer_group = []
        self.similarity_cutoff_used = None

    def get_duplicate_groups(self) -> List[List['AssembledIsomer']]:
        """
        Returns all input isomers as a list of lists, in which each sublist contains isomers that are considered duplicates.
        """
        if len(self.isomers) <= 1:
            return [self.isomers]

        if self.method == "alignment":
            self.isomer_group = self._reduce_by_alignment()
        elif self.method == "distances":
            self.isomer_group = self._reduce_by_fingerprint()
        else:
            raise ValueError(f"Fatal Error: Unsupported reduction method '{self.method}. Supported methods are 'alignment' and 'distances'.")

        return self.isomer_group

    def _reduce_by_alignment(self) -> List[List['AssembledIsomer']]:
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
        self.isomer_group = self._group_isomers(group_labels_matrix, self.isomers)
        return self.isomer_group

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
        assert hasattr(self, "metal_centers"), ValueError("Fatal Error: metal_centers must be defined before calling align_isomers. ")
        logging.debug(f"Aligning isomers based on {len(self.unique_metal_centers)} metal centre(s).")

        # Here the number of metal centres and the fact that metal centres of different isomers
        # must always be aligned is taken into account to determine if the isomers are similar or not
        if len(self.unique_metal_centers) == 1:
            # There are 3 axes which an isomer can be rotated around to align it with another isomer
            bounds = [[0, 360] for _ in range(3)]
            # The 3 cardinal axes, properly centered on the metal centre
            axes = np.eye(3)  # Standard Cartesian axes (x, y, z)
            logging.debug("Performing 3D brute-force alignment over x, y, z axes.")

        elif len(self.unique_metal_centers) == 2:
            # The isomer can only be rotated around the metal-metal axis to determine if the isomers are similar or not
            bounds = [[0, 360] for _ in range(1)]
            axis_vector = self.unique_metal_centers[1].position - self.unique_metal_centers[0].position
            axis_vector /= np.linalg.norm(axis_vector)  # Normalize
            axes = [axis_vector]
            logging.debug(f"Performing 1D brute-force alignment around axis: {axis_vector.tolist()}")

        elif len(self.unique_metal_centers) >= 3:
            # 3 metal centres means each isomer is fixed in space and their archetypes can be directly compared
            # We simply return the energy heuristic
            logging.debug("Three or more metal centres detected — skipping alignment and using direct heuristic.")
            return self.energy_heuristic(stat_atoms=stationary_atoms, rot_atoms=rotated_atoms)

        else:
            raise ValueError(f"Fatal Error: Unsupported number of metal centres ({len(self.unique_metal_centers)}). ")

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
                                     center=np.array(self.unique_metal_centers[0].position))

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
        self.isomer_group = self._group_isomers(group_labels_matrix, self.isomers)
        return self.isomer_group

    @staticmethod
    def _group_isomers(group_labels_matrix: np.ndarray, isomers: List['AssembledIsomer']) -> List[List['AssembledIsomer']]:
        """
        Group isomers based on which ones are labeled 'Close' in the group_labels_matrix.
        :param group_labels_matrix: 2D numpy array with labels 'Close' or 'Far'. It's possible that this list is not consistent, i.e. if isomers 1 and 2 are 'Close', and isomers 2 and 3 are 'Close', it is still possible that isomers 1 and 3 are 'Far', even though this is illogical. This will be handled here such that isomers 1 and 3 are grouped together as well in this case.
        :param isomers: List of AssembledIsomer objects to be grouped.
        :return: List of lists, where each sublist contains isomers that are considered duplicates.
        """
        n = len(isomers)
        if group_labels_matrix.shape != (n, n):
            raise ValueError(
                "group_labels_matrix must be a square matrix with the same "
                "dimension as the number of isomers."
            )

        # Build an undirected adjacency list
        adjacency = [[] for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if (
                        group_labels_matrix[i, j] == "Close"
                        or group_labels_matrix[j, i] == "Close"
                ):
                    adjacency[i].append(j)
                    adjacency[j].append(i)

        # Depth-first search for connected components
        visited = [False] * n
        groups: List[List["AssembledIsomer"]] = []

        for start in range(n):
            if visited[start]:
                continue
            stack = [start]
            component_indices: List[int] = []

            while stack:
                node = stack.pop()
                if visited[node]:
                    continue
                visited[node] = True
                component_indices.append(node)
                stack.extend(adjacency[node])

            # Maintain original order inside the component
            component_indices.sort()
            groups.append([isomers[idx] for idx in component_indices])

        assert len([_ for sublist in groups for _ in sublist]) == n, "Grouped indices do not cover all isomers or have duplicates."

        return groups

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
        fingerprints = [self._compute_sorted_distance_fingerprint(isomer)[0] for isomer in isomers]
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
    def _compute_sorted_distance_fingerprint(isomer) -> tuple[np.ndarray, list[tuple[str, str]]]:
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

        #  also get sorted element pairs
        sorted_pairs = list(zip(first_elem[order], second_elem[order]))

        return sorted_dists, sorted_pairs

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
            self.get_duplicate_groups()
            for isomer_group in self.isomer_group:
                for i in range(1, len(isomer_group)):
                    isomer_group[i].warning = 'duplicate'  # Mark all but the first isomer in the group as duplicate.

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
        if "plotly" not in sys.modules:
            print("Plotly is not installed. Please install it to use the interactive heatmap feature.")
            return None
        app = dash.Dash(__name__)
        fig = px.imshow(
            df,
            color_continuous_scale=color_scale,
            zmin=0, zmax=1.0,
            title="Comparison Matrix (Interactive)",
            labels={"x": "AssembledIsomer Index", "y": "AssembledIsomer Index", "color": "Difference"}
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
                print(f"Clicked cell: ({i_idx}, {j_idx}) — opening viewer.")
                # Direct call to viewer
                self.view_isomer_alignment(i_idx, j_idx, grid_size=self.grid_size)
            return fig

        app.run(
            debug=True,
            use_reloader=False,
            **{
                "threaded": False,
                "processes": 1,
                "use_debugger": True,
                "dev_tools_silence_routes_logging": False,
                "dev_tools_prune_errors": False
            }
        )

    def view_isomer_alignment(self, index1: int, index2: int, grid_size=None) -> None:
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

        print("1")
        isomer1 = self.isomers[index1].atoms.copy()
        print("2")
        isomer2 = self.isomers[index2].atoms.copy()
        print("3")
        isomer2_aligned = self.isomers[index2].atoms.copy()
        print("4")

        # Determine alignment axes
        if len(self.unique_metal_centers) == 1:
            bounds = [[0, 360] for _ in range(3)]
            axes = np.eye(3)
        elif len(self.unique_metal_centers) == 2:
            bounds = [[0, 360]]
            axis_vector = self.unique_metal_centers[1].position - self.unique_metal_centers[0].position
            axis_vector /= np.linalg.norm(axis_vector)
            axes = [axis_vector]
        elif len(self.unique_metal_centers) >= 3:
            bounds = None
            axes = None
            print("Three or more metal centres — skipping rotation; showing structures unaligned.")
        else:
            raise ValueError("Fatal Error: Invalid number of metal centres for alignment.")

        # Perform rotation if applicable
        if len(self.unique_metal_centers) < 3:
            result_angles = brute(
                func=self.objective_function,
                ranges=bounds,
                args=(isomer1, isomer2_aligned, axes),
                Ns=grid_size if grid_size else self.grid_size
            )
            R_total = self.combined_rotation_matrix(result_angles, axes)
            self.apply_combined_rotation(atoms=isomer2_aligned, R_total=R_total, center=np.array(self.unique_metal_centers[0].position))

        # Create overlaid image: isomer1 + rotated isomer2
        overlaid = isomer1.copy() + isomer2_aligned.copy()

        # Assign tags so colors are distinguishable
        for atom in isomer1:
            atom.tag = 1
        for atom in isomer2:
            atom.tag = 2
        for atom in overlaid:
            atom.tag = 3

        print("Launching ASE viewer with aligned isomers...")
        view([isomer1, isomer2, overlaid], viewer="ase")

    def debug_fingerprints(self, idx1, idx2):
        """
        A tool that displays the fingerprints of two isomers side by side and highlights
        which difference is used to generate the score
        :param idx1: index of first isomer
        :param idx2: index of second isomer
        :return: None
        """

        # Get the fingerprints of the two isomers
        fp1, pairs1 = self._compute_sorted_distance_fingerprint(self.isomers[idx1])
        fp2, pairs2 = self._compute_sorted_distance_fingerprint(self.isomers[idx2])

        # Ensure both fingerprints are of the same length
        if len(fp1) != len(fp2):
            raise ValueError(f"Fingerprints of isomers {idx1} and {idx2} have different lengths: {len(fp1)} vs {len(fp2)}")

        # Create a DataFrame for visualization
        df = pd.DataFrame({
            "Distance": np.arange(len(fp1)),
            f"Isomer {idx1}": fp1,
            f"Isomer {idx2}": fp2,
            "Pair": [f"{pair[0]}-{pair[1]}" for pair in pairs1],
            "Difference": np.abs(fp1 - fp2),
        })
        # display the results
        print(df)

        return df


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
        if n <= 1:  # nothing to clash
            return False

        pos = atoms.positions  # (N, 3)
        symbols = atoms.get_chemical_symbols()  # list[str]

        # Get vectorized covalent radii
        radii_map = {s: Element(s).covalent_radius_angstrom for s in set(symbols)}
        radii = np.fromiter((radii_map[s] for s in symbols), dtype=float)

        # Get arrays of pairwise distances and minimum allowed distances
        dists = pdist(pos)
        i_idx, j_idx = np.triu_indices(n, k=1)  # matches pdist order
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

        if not np.any(mask):  # nothing left to check
            return False

        # Check if any distances are below the allowed minimum which means atoms are clashing
        has_clashes = np.any(dists[mask] < min_allowed[mask])

        return has_clashes

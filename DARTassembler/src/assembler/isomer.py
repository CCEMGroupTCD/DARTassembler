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
    """
    Represent an assembled transition-metal complex isomer.

    This class wraps the assembled ASE Atoms object together with connectivity graph,
    ligand/metal indices and metadata for downstream filtering and output.

    :param atomic_props: ASE Atoms object or dict of atomic properties describing the complex.
    :type atomic_props: Union[ase.Atoms, dict[str, Any]]
    :param graph: NetworkX graph representing atomic connectivity and node labels.
    :type graph: nx.Graph
    :param metal_idc: Indices of metal center atoms in the ASE Atoms object.
    :type metal_idc: list[int]
    :param donor_idc: Flattened indices of donor atoms in the merged graph.
    :type donor_idc: list[list[int]]
    :param ligand_idc: Indices of atoms belonging to each ligand in the merged graph.
    :type ligand_idc: list[list[int]]
    :param ligand_info: Dictionary containing ligand metadata required to reconstruct Ligand objects.
    :type ligand_info: dict[str, Any] | None
    :param global_props: Global properties for the molecule.
    :type global_props: dict[str, Any] | None
    :param validity_check: If True, run transition-metal-complex specific validity checks.
    :type validity_check: bool
    :param target_vectors: Target donor vectors used to assemble the isomer (one list per ligand).
    :type target_vectors: list[list[list[float]]] | None
    :param ligand_origins: Origin coordinates (shape (3,)) used to translate each ligand to its metal center.
    :type ligand_origins: list[list[float]] | None
    :param warning: Warning string (e.g., 'clashing' or 'duplicate(...)') to annotate the isomer.
    :type warning: str
    :param isomer_name: Unique name assigned to the isomer for identification.
    :type isomer_name: str | None
    """

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
        """
        Construct Ligand objects from stored ligand metadata and ASE Atoms slices.

        The method builds Ligand() wrappers for each ligand using the ligand indices
        and ligand_info provided during assembly. These Ligand objects are convenient
        for downstream geometry operations and will not be persisted by default.

        :param validity_check: Whether to run validity checks when instantiating Ligand objects.
        :type validity_check: bool
        :return: List of Ligand objects corresponding to each ligand in the isomer.
        :rtype: list[ Ligand ]
        """
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
        Return the chemical symbols of the metal center atoms.

        :return: Chemical symbols of metals in the same order as metal_idc.
        :rtype: list[str]
        """
        return self.atoms[self.metal_idc].get_chemical_symbols()

    def _tmc_validity_checks(self) -> None:
        """
        Run transition-metal-complex specific validity checks.

        Checks basic molecular validity via the BaseMolecule helper and warns if any
        declared metal center is not classified as a metal element by the Element utility.

        :raises AssertionError: If base molecule validity checks fail.
        :return: None
        :rtype: None
        """
        self._check_if_molecule_valid()  # Checks basic molecular properties like atomic_props and graph
        # Doublecheck if all the metals are really metals. Don't raise an error in case it's intentional.
        all_metals = all(Element(metal).is_metal for metal in self.metals)
        if not all_metals:
            warnings.warn(f"Any of the metal centers {self.metals} in AssembledIsomer() is not a metal. Providing a chemical element as `metal center` that is not a metal is not a problem for DART, this warning is just to make you aware.")

        return

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the AssembledIsomer to a serializable dictionary.

        The returned dictionary includes atomic properties, connectivity graph (handled
        by the base class), and assembly-specific metadata such as metal/donor/ligand indices.

        :return: Dictionary representation suitable for JSON or other serialization.
        :rtype: dict[str, Any]
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
        Return an ASE Atoms object containing only the metal center atoms.

        :return: ASE Atoms with the metal center atoms (order corresponds to metal_idc).
        :rtype: ase.Atoms
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
    """
    Assemble transition-metal complex isomers from Ligand objects, target vectors and metal centers.

    The class stores the input ligands, their target donor vectors and metal center definitions and
    provides methods to generate assembled isomers, perform mono-axial optimization and filter duplicates.
    """

    def __init__(
            self,
            ligands: List[Ligand],
            target_vectors: list[list[list[float]]],
            metal_centers: Union[List[List[Union[str, List[float]]]], str],
            ligand_origins: List[List[float]] = None,
    ):
        """
        The assembler will align and translate each ligand such that its donor atoms point along the
        provided target_vectors and place metal center atoms at the specified coordinates.
        If metal_centers is a string (element symbol), a single metal center at the origin is assumed
        for each ligand (mono-metallic complex).

        :param ligands: Sequence of Ligand objects to assemble.
        :type ligands: list[ Ligand ]
        :param target_vectors: For each ligand, a list of donor target vectors (each vector length 3).
        :type target_vectors: list[list[list[float]]]  # shape: (n_ligands, n_donors_per_ligand, 3)
        :param metal_centers: Either a symbol string for a single metal at origin, or a list where each entry
                              corresponds to the metal center(s) the ligand binds to. For the latter, each
                              metal is a tuple [element_symbol, [x,y,z]].
        :type metal_centers: Union[str, list[list[ Union[str, list[float]] ]]]
        :param ligand_origins: Optional list of origin coordinates (3 floats) for each ligand; if None,
                               defaults to the corresponding metal center positions (or the mean for bridging ligands).
        :type ligand_origins: list[list[float]] | None


        :param ligand_origins: Optional list of origin coordinates (3 floats) for each ligand; if None,
                               defaults to the corresponding metal center positions (or the mean for bridging ligands).
        :type ligand_origins: list[list[float]] | None

        Example usage for assembling a mono-metallic square-planar Pd complex with 2 cis bidentate ligands in the xy-plane:
            factory = AssembledComplex(
                                    ligands=..., (list of two bidentate Ligand() objects, usually from the MetaLig database)
                                    metal_centers='Pd',
                                    target_vectors=[
                                                    [[1, 0, 0], [0, 1, 0]],     # the donor atoms of the first bidentate ligand are oriented in (+x,+y) direction
                                                    [[-1, 0, 0], [0, -1, 0]],   # the donor atoms of the second bidentate ligand are oriented in (-x,-y) direction
                                                    ],
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
        """
        Validate inputs and set sensible defaults for metal_centers and ligand_origins.

        This method:
          - Expands a single-element metal_centers string into ASE Atoms at the origin for each ligand.
          - Converts provided metal center specifications into ASE Atom objects.
          - Sets default ligand_origins to the metal center(s) position (or mean position for bridging ligands).
          - Validates dimensions and numeric types of target_vectors and matches them against ligand archetypes,
            issuing warnings if the provided vectors do not match the ligand archetype.

        :param metal_centers: Metal center specification (see __init__).
        :type metal_centers: Union[str, list[list[ Union[str, list[float]] ]]]
        :param ligand_origins: Optional ligand origins to use for translation.
        :type ligand_origins: list[list[float]] | None  # shape: (n_ligands, 3)
        :param ligands: List of Ligand objects.
        :type ligands: list[ Ligand ]
        :param target_vectors: Target donor vectors per ligand.
        :type target_vectors: list[list[list[float]]]  # shape: (n_ligands, n_donors, 3)
        :raises ValueError: If list lengths mismatch or target vector shapes are invalid.
        :raises TypeError: If target_vectors has an invalid type structure.
        :return: Normalized (metal_centers, ligand_origins, ligands, target_vectors)
        :rtype: tuple[ list[list[ase.Atom]], list[list[float]], list[ Ligand ], list[list[list[float]]] ]
        """
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
                         permutable_ligands: Optional[List[int]] = None,
                         monoaxial_optimization: Optional[bool] = True,
                         force_all_isomers: bool = False,
                         clashing_tolerance: float = -0.3,
                         clashing_metal: bool = False,
                         duplicate_tolerance: float = 0.5,
                         complex_name_length: int = 8,
                         complex_name_suffix: str = '',
                         avoid_names: Optional[Iterable[str]] = None,
                         ):
        """
        Generate assembled isomers from the provided ligands and metal centers.

        The method constructs all geometric combinations by assigning rotated ligand instances
        to metal centers, applies optional mono-axial optimization, filters clashing isomers,
        and groups duplicates via fingerprint or alignment-based comparison.

        :param permutable_ligands: Indices of ligands allowed to be permuted when generating isomers.
        :type permutable_ligands: list[int] | None
        :param monoaxial_optimization: Whether to perform mono-axial rotation optimization for mono-coordinating ligands.
        :type monoaxial_optimization: bool
        :param force_all_isomers: If True, generate all geometrical isomers for each ligand archetype.
        :type force_all_isomers: bool
        :param clashing_tolerance: Distance buffer below covalent radii sum allowed before flagging a clash (Å).
        :type clashing_tolerance: float | None
        :param clashing_metal: If True, also check ligand-metal and metal-metal clashes.
        :type clashing_metal: bool
        :param duplicate_tolerance: Cutoff used by duplicate grouping (fingerprint threshold).
        :type duplicate_tolerance: float
        :param complex_name_length: Length of the generated complex name seed.
        :type complex_name_length: int
        :param complex_name_suffix: Optional suffix appended to generated complex names.
        :type complex_name_suffix: str
        :param avoid_names: Iterable of names to avoid when generating unique complex names.
        :type avoid_names: Iterable[str] | None
        :return: None. Results are stored on the AssembledComplex instance as `isomers`,
                 `successful_isomers`, `unsuccessful_isomers` and `success` boolean.
        :rtype: None
        """
        if avoid_names is None:
            avoid_names = set()
        self.clashing_tolerance = clashing_tolerance
        self.clashing_metal = clashing_metal
        self.duplicate_tolerance = duplicate_tolerance
        self.permutable_ligands = permutable_ligands
        self.monoaxial_optimization = monoaxial_optimization
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
        target_vector_combs = get_list_with_all_possible_swappings(objects=self.target_vectors, permutable_ligands=self.permutable_ligands)
        ligand_origin_combs = get_list_with_all_possible_swappings(objects=self.ligand_origins, permutable_ligands=self.permutable_ligands)

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

        pre_isomers_duplicate_groups = DuplicateIsomerFilter(isomers=isomers, fingerprint_grouping_cutoff=self.duplicate_tolerance, metal_centers=self.metal_centers).get_duplicate_groups()
        pre_isomers_duplicate_group_names = [set([isomer.isomer_name for isomer in isomer_group]) for isomer_group in pre_isomers_duplicate_groups]

        # Do a mono-axial optimization of the isomers and afterward check for clashing ligands.
        for idx, isomer, target_vectors, ligand_origins in zip(range(len(isomers)), isomers, same_length_target_vectors, same_length_ligand_origins):
            isomer = AxialOptModifier(isomers=[isomer], opt=self.monoaxial_optimization).modify(target_vectors_list=[target_vectors], ligand_origins_list=[ligand_origins])[0]
            isomers[idx] = isomer # important: copy changed object over to list

            if isomer.warning == '' and pd.notna(self.clashing_tolerance):
                clashfilter = IsomerClashFilter(buffer=self.clashing_tolerance, check_metal_clashes=self.clashing_metal)
                clashing = clashfilter.has_clashing_atoms(atoms=isomer.atoms, ligand_idc=isomer.ligand_idc, metal_idc=isomer.metal_idc)
                if clashing:
                    isomer.warning = 'clashing'  # this now updates the optimized object

        # Check for duplicates again after the mono-axial optimization.
        joined_isomers_duplicate_group_names = self.get_duplicate_isomers_group_names(isomers=isomers, pre_isomers_duplicate_group_names=pre_isomers_duplicate_group_names, duplicate_tolerance=self.duplicate_tolerance,  metal_centers=self.metal_centers)

        self.successful_isomers, self.unsuccessful_isomers = self.divide_into_successful_and_unsuccessful_isomers(isomers, joined_isomers_duplicate_group_names)
        self.success = len(self.successful_isomers) > 0
        self.isomers = isomers

        return

    @staticmethod
    def get_duplicate_isomers_group_names(isomers: list[Any], pre_isomers_duplicate_group_names: list[set[str | None]], duplicate_tolerance: float, metal_centers: list[list[Atom]]) -> list[list[str]]:
        """
        Join pre- and post-optimization duplicate groups and preserve isomer ordering.

        Uses union of pre- and post-monoaxial-optimization duplicate sets and sorts each
        joined group according to the original `isomers` order so the first element can be retained.

        :param isomers: List of AssembledIsomer objects.
        :type isomers: list[Any]
        :param pre_isomers_duplicate_group_names: List of sets containing isomer names flagged as duplicates before optimization.
        :type pre_isomers_duplicate_group_names: list[ set[str] ]
        :param duplicate_tolerance: Cutoff used by duplicate grouping (fingerprint threshold).
        :type duplicate_tolerance: float
        :param metal_centers: Metal center atom definitions used for duplicate comparison.
        :type metal_centers: list[list[ase.Atom]]
        :return: Ordered list of duplicate groups where each group is a list of isomer names.
        :rtype: list[list[str]]
        """
        post_isomers_duplicate_groups = DuplicateIsomerFilter(isomers=isomers, fingerprint_grouping_cutoff=duplicate_tolerance, metal_centers=metal_centers).get_duplicate_groups()
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
        Split isomers into successful and unsuccessful sets based on warnings and duplicate groups.

        Successful isomers are the first non-clashing member of each duplicate group. Others are
        marked unsuccessful and labelled with a suitable warning ('clashing' or 'duplicate(...)').

        :param isomers: List of AssembledIsomer objects (each must have 'isomer_name' and 'warning' attributes).
        :type isomers: list[Any]
        :param joined_isomers_duplicate_group_names: Joined duplicate groups with isomer names (post optimization).
        :type joined_isomers_duplicate_group_names: list[list[str]]
        :return: Tuple (successful_isomers, unsuccessful_isomers).
        :rtype: tuple[list[Any], list[Any]]
        :raises AssertionError: If duplicate group contains duplicate names or is empty.
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
        Convert the AssembledComplex factory state to a dictionary.

        The resulting dictionary contains the generated isomers (by name), the merged
        connectivity graph, graph hash, index mappings and the input parameters used to assemble.

        :return: Dictionary containing assembly metadata and serialized isomer entries.
        :rtype: dict[str, Any]
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
                "clashing_tolerance": self.clashing_tolerance,
                "clashing_metal": self.clashing_metal,
                "duplicate_tolerance": self.duplicate_tolerance,
                "permutable_ligands": self.permutable_ligands,
                "monoaxial_optimization": self.monoaxial_optimization,
                "force_all_isomers": self.force_all_isomers,
                "complex_name_length": self.complex_name_length,
                "complex_name_suffix": self.complex_name_suffix,
            }

        }

    def _get_complex_name(self, avoid_names: Optional[Iterable[str]]) -> str:
        """
        Generate a (pseudo-)random complex name based on the graph hash.

        :param avoid_names: Iterable of names that should be avoided when generating a new name.
        :type avoid_names: Iterable[str] | None
        :return: Generated complex name string.
        :rtype: str
        """
        return get_complex_name(seed=self.graph_hash, length=self.complex_name_length, suffix=self.complex_name_suffix, avoid_names=avoid_names)

    def _get_ligandinfo(self) -> Dict[str, Any]:
        """
        Return ligand metadata required to reconstruct Ligand objects in AssembledIsomer.

        :return: Dictionary with keys such as 'unique_names', 'archetypes', 'donor_idcs', etc.
        :rtype: dict[str, Any]
        """
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
        Produce rotated and translated ASE Atoms for each ligand according to target vectors.

        Each ligand yields a list of ASE Atoms objects representing the possible geometric isomers
        compatible with the supplied target_vectors. Returned structure is a list of lists:
        outer list over ligands, inner list over isomer instances (ASE Atoms).

        :param target_vectors: List of target vectors for each ligand.
        :type target_vectors: list[list[list[float]]]  # shape: (n_ligands, n_donors, 3)
        :param ligand_origins: List of origin coordinates used to translate each ligand instance.
        :type ligand_origins: list[list[float]]  # shape: (n_ligands, 3)
        :return: Rotated and translated ligand isomer instances.
        :rtype: list[list[ ase.Atoms ]]
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
        Return a list of unique ASE Atom metal centers used across all ligand definitions.

        The method deduplicates identical metal atoms (same element and position).

        :return: List of unique ASE Atom objects for all metal centers.
        :rtype: list[ ase.Atom ]
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
        Merge per-ligand connectivity graphs and connect donor atoms to metal nodes.

        Returns a merged NetworkX Graph in which metal centers occupy the first node indices
        followed by ligand atom nodes. Also returns ligand atom index lists (per ligand) and
        flattened donor index list.

        :return: Tuple (graph, ligand_indices, donor_idc)
        :rtype: tuple[ nx.Graph, list[list[int]], list[int] ]
        :raises AssertionError: If merged graph is not connected or node labels mismatch expected atom list.
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
    """
    Optimize mono-coordinating ligand rotations around their coordination axis.

    Uses a global differential evolution optimizer to find per-ligand rotation angles
    minimizing a distance-based penalty (short interatomic distances penalized).
    """

    def __init__(self, isomers: List['AssembledIsomer'], opt: bool = True, distance_cutoff: Optional[float] = 4.0, use_cutoff: bool = False):
        """
        Initialize the mono-axial optimization modifier.

        :param isomers: List of AssembledIsomer objects to optimize.
        :type isomers: list[ AssembledIsomer ]
        :param opt: Whether to perform optimization or return inputs unchanged.
        :type opt: bool
        :param distance_cutoff: Optional Å cutoff to restrict penalty to interatomic distances below this value.
        :type distance_cutoff: float | None
        :param use_cutoff: Whether to apply the provided distance_cutoff.
        :type use_cutoff: bool
        """
        self.input_isomers = isomers
        self.opt_command = opt
        # Optional distance cutoff for the objective function (None = use all pairs)
        self.distance_cutoff: Optional[float] = distance_cutoff if use_cutoff else None
        self.output_isomers = []
        logging.debug(f"AxialOpt initialized with {len(self.input_isomers)} AssembledIsomer objects.")

    def modify(self, target_vectors_list, ligand_origins_list, maxiter=1000, popsize=15) -> List['AssembledIsomer']:
        """
        Optimize rotation angles for each provided isomer independently.

        For each isomer, a separate differential evolution run is executed to find the set
        of ligand rotation angles that minimize the interatomic distance penalty. Only ligands
        with archetypes '1-mono' or '2-trans' are rotated.

        :param target_vectors_list: List of target_vectors corresponding to each isomer in input_isomers.
        :type target_vectors_list: list[list[list[list[float]]]]  # shape: (n_isomers, n_ligands, n_donors, 3)
        :param ligand_origins_list: Corresponding ligand origins for each isomer.
        :type ligand_origins_list: list[list[list[float]]]  # shape: (n_isomers, n_ligands, 3)
        :param maxiter: Maximum iterations passed to differential_evolution.
        :type maxiter: int
        :param popsize: Population size passed to differential_evolution.
        :type popsize: int
        :return: List of AssembledIsomer objects with updated atom positions and atomic_props.
        :rtype: list[ AssembledIsomer ]
        :raises ValueError: If lengths of input_isomers, target_vectors_list, and ligand_origins_list differ.
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
        Compute the penalty for a given set of rotation angles.

        The penalty is computed as the sum over 1.0 / d^2 for all atomic pairs (or pairs within the distance_cutoff),
        where d is the interatomic distance after applying ligand rotations defined by x.
        Rotations are only applied to ligands whose archetype is in ['1-mono', '2-trans'].

        :param x: Array of rotation angles in degrees for each ligand.
        :type x: np.ndarray  # shape: (n_rotations,)
        :param vectors_in: List of rotation axes (unit or non-unit vectors) per ligand.
        :type vectors_in: list[ np.ndarray ]  # each of shape (3,) or (2,3) for 2-trans
        :param origins_in: Origins (3 floats) for each ligand rotation.
        :type origins_in: list[list[float]]  # shape: (n_rotations, 3)
        :param TMC_in: ASE Atoms object representing the assembled complex to act upon.
        :type TMC_in: ase.Atoms
        :param ligand_idc: List of lists containing atom indices for each ligand in the complex.
        :type ligand_idc: list[list[int]]
        :param archetypes: List of ligand archetypes in the same order as ligand_idc.
        :type archetypes: list[str]
        :return: Scalar penalty value to minimize.
        :rtype: float
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
        Rotate selected atoms around a given axis by the specified angle (degrees).

        Accepts axis definitions of shape (3,), (1,3) or (2,3). For (2,3) the first row is used (2-trans).
        Rotation is performed in-place on the provided ASE Atoms object.

        :param atoms: ASE Atoms object whose subset of atoms will be rotated in-place.
        :type atoms: ase.Atoms
        :param vector: Rotation axis vector (or matrix of axis vectors).
        :type vector: np.ndarray  # shape: (3,) or (1,3) or (2,3)
        :param origin: Rotation center coordinates.
        :type origin: list[float]  # length 3
        :param idc: Indices of atoms to rotate.
        :type idc: list[int]
        :param angle: Rotation angle in degrees.
        :type angle: int | float
        :return: The modified ASE Atoms object (same object passed in).
        :rtype: ase.Atoms
        :raises ValueError: If the axis has zero length or unexpected shape.
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
        Visualize input and optimized output complexes interleaved using ASE's viewer.

        The ASE viewer will show alternating frames: Input_0, Optimized_0, Input_1, Optimized_1, ...
        If no optimization has been performed, the method prints a message and returns.

        :return: None
        :rtype: None
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
    Reduce the number of assembled isomers by detecting duplicates via fingerprint or alignment.

    The class supports two main deduplication strategies:
      - 'distances' (fingerprint-based, fast)
      - 'alignment' (rotational alignment and heuristic comparison)
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
        Initialize the duplicate isomer filter parameters.

        :param isomers: List of AssembledIsomer objects to analyze.
        :type isomers: list[ AssembledIsomer ]
        :param method: Deduplication strategy: 'alignment' or 'distances'.
        :type method: str
        :param grid_size: Grid density used by brute-force alignment search (points per angle).
        :type grid_size: int
        :param isomer_comparison_mode: Mode for fingerprint comparison ('max_diff', 'sum_diff', 'mean_diff', 'rmsd').
        :type isomer_comparison_mode: str
        :param isomer_comparison_grouping_mode: Grouping approach: 'cluster' or 'cutoff'.
        :type isomer_comparison_grouping_mode: str
        :param fingerprint_grouping_cutoff: Numeric cutoff used in 'cutoff' grouping mode. Use NaN to disable.
        :type fingerprint_grouping_cutoff: float
        :param metal_centers: Metal center definitions used to determine rotation degrees of freedom.
        :type metal_centers: list[list[ ase.Atom ]] | None
        :param energy_heuristic_mode: Mode for energy heuristic in alignment ('max' or 'sum').
        :type energy_heuristic_mode: str
        """
        if metal_centers is None:
            metal_centers = []
        self.isomers = isomers
        self.method = method.lower()
        self.grid_size = grid_size  # The number of grid points when scanning from 0 to 360 for exact alignment
        self.isomer_comparison_mode = isomer_comparison_mode
        self.isomer_comparison_grouping_mode = isomer_comparison_grouping_mode
        self.fingerprint_grouping_cutoff = fingerprint_grouping_cutoff
        self.check_duplicate = pd.notna(fingerprint_grouping_cutoff)
        self.metal_centres = metal_centers
        self.unique_metal_centers = list({(atom.symbol, tuple(atom.position)): atom for sublist in self.metal_centres for atom in sublist}.values())
        self.diff_matrix = None  # Placeholder for the fingerprint difference matrix
        self.energy_heuristic_mode = energy_heuristic_mode
        self.isomer_group = []
        self.similarity_cutoff_used = None

    def get_duplicate_groups(self) -> List[List['AssembledIsomer']]:
        """
        Group input isomers into duplicate clusters.

        If duplicate checking is disabled or only one isomer exists, returns trivial groups.
        Otherwise, delegates to the configured reduction method.

        :return: List of groups; each group is a list of AssembledIsomer objects considered duplicates.
        :rtype: list[list[ AssembledIsomer ]]
        """
        if len(self.isomers) <= 1:
            return [self.isomers]
        if not self.check_duplicate:
            # If duplicate checking is disabled, return each isomer as its own group.
            return [[isomer] for isomer in self.isomers]

        if self.method == "alignment":
            self.isomer_group = self._reduce_by_alignment()
        elif self.method == "distances":
            self.isomer_group = self._reduce_by_fingerprint()
        else:
            raise ValueError(f"Fatal Error: Unsupported reduction method '{self.method}. Supported methods are 'alignment' and 'distances'.")

        return self.isomer_group

    def _reduce_by_alignment(self) -> List[List['AssembledIsomer']]:
        """
        Reduce isomers by searching for an alignment that minimizes a distance heuristic.

        Rotational alignment between pairs is performed using a brute-force grid over
        allowed rotation axes and angles. The resulting pairwise scores are clustered
        to form duplicate groups.

        :return: Duplicate groups as lists of AssembledIsomer objects.
        :rtype: list[list[ AssembledIsomer ]]
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
        Compute a heuristic 'energy' that quantifies similarity between two isomers.

        For each element type the method solves an optimal assignment between like atoms
        (Hungarian algorithm) and either returns the sum of assigned distances ('sum' mode)
        or the maximum assigned distance among element-types ('max' mode).

        :param stat_atoms: Stationary ASE Atoms object (reference).
        :type stat_atoms: ase.Atoms
        :param rot_atoms: Rotated ASE Atoms object (to be compared).
        :type rot_atoms: ase.Atoms
        :return: Energy heuristic: scalar representing sum or maximum of assigned distances.
        :rtype: float
        :raises ValueError: If elements mismatch or counts per element differ between isomers.
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
        Align two isomers by searching for an alignment that minimizes a distance heuristic.

        Depending on the number of unique metal centres, the function reduces the rotational
        degrees of freedom (3D rotation for one metal, 1D rotation around metal-metal axis for two metals,
        and direct heuristic for three or more metals).

        :param stationary_atoms: ASE Atoms kept fixed during alignment.
        :type stationary_atoms: ase.Atoms
        :param rotated_atoms: ASE Atoms that will be rotated to align with stationary_atoms.
        :type rotated_atoms: ase.Atoms
        :return: Alignment score computed by objective_function and energy_heuristic.
        :rtype: float
        :raises AssertionError: If metal_centers was not set prior to calling this method.
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
        Objective used by the brute-force alignment routine.

        Applies combined rotations described by x around axes and returns the energy heuristic
        comparing atoms1 (reference) and the rotated version of atoms2.

        :param x: Rotation angles in degrees.
        :type x: np.ndarray
        :param atoms1: Stationary reference ASE Atoms.
        :type atoms1: ase.Atoms
        :param atoms2: ASE Atoms to be rotated.
        :type atoms2: ase.Atoms
        :param axes: List/array of rotation axes.
        :type axes: np.ndarray | list[np.ndarray]
        :return: Scalar value from energy_heuristic to minimize.
        :rtype: float
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
        Compute a combined 3x3 rotation matrix from sequences of angles and axes.

        :param angles: Iterable of angles in degrees.
        :type angles: Iterable[ float ]
        :param axes: Corresponding iterable of 3D axis vectors.
        :type axes: Iterable[ list[float] ] | np.ndarray
        :return: Combined rotation matrix (3x3).
        :rtype: np.ndarray
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
        Apply a combined rotation matrix to an ASE Atoms object about a center.

        Positions are rotated in-place.

        :param atoms: ASE Atoms to rotate.
        :type atoms: ase.Atoms
        :param R_total: 3x3 rotation matrix.
        :type R_total: np.ndarray
        :param center: Center of rotation (3 floats).
        :type center: list[float] | np.ndarray
        :return: None
        :rtype: None
        """
        # Shift positions relative to center, apply rotation, then shift back.
        shifted = atoms.positions - center
        atoms.positions = center + (shifted @ R_total.T)

    def _reduce_by_fingerprint(self):
        """
        Reduce isomers using a sorted interatomic-distance fingerprint and clustering.

        Computes a symmetric difference matrix and clusters pairs labelled 'Close' to form groups.

        :return: List of duplicate groups (lists of AssembledIsomer objects).
        :rtype: list[list[ AssembledIsomer ]]
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
        Group isomers into connected components based on a pairwise label matrix.

        The input matrix should contain labels 'Close'/'Far' for each pair; connected 'Close' entries
        are grouped into duplicate clusters. The method ensures that groups preserve original ordering.

        :param group_labels_matrix: Square 2D numpy array with values 'Close' or 'Far'.
        :type group_labels_matrix: np.ndarray
        :param isomers: List of AssembledIsomer objects to group.
        :type isomers: list[ AssembledIsomer ]
        :return: Grouped isomers as a list of lists preserving original order within groups.
        :rtype: list[list[ AssembledIsomer ]]
        :raises ValueError: If group_labels_matrix is not square or does not match number of isomers.
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
        Convert a pairwise difference matrix to a matrix of labels ('Close'/'Far').

        Two modes are supported:
          - 'cluster': MeanShift clustering on upper-triangle values to identify the 'Close' cluster.
          - 'cutoff': Hard cutoff thresholding of matrix values.

        :param matrix: Square pairwise difference matrix.
        :type matrix: np.ndarray
        :param quantile: Quantile used to estimate bandwidth in clustering mode (unused default kept).
        :type quantile: float
        :param method: 'cluster' or 'cutoff' defining grouping strategy.
        :type method: str
        :param cutoff: Numeric cutoff required when method='cutoff'.
        :type cutoff: float | None
        :return: Square matrix with entries 'Close' or 'Far'.
        :rtype: np.ndarray
        :raises ValueError: If cutoff is None when method='cutoff'.
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
        Compute symmetric fingerprint difference matrix for a list of isomers.

        :param isomers: List of AssembledIsomer objects.
        :type isomers: list[ AssembledIsomer ]
        :return: Symmetric numpy array (n,n) of pairwise fingerprint differences.
        :rtype: np.ndarray
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
        Compute a sorted inter-atomic distance fingerprint for an isomer.

        The fingerprint vector is the ordered list of pairwise distances between atoms,
        sorted first by element-pair (lexicographically) and then by ascending distance.
        Also returns the corresponding ordered element-pair labels.

        :param isomer: AssembledIsomer object with an `.atoms` ASE Atoms attribute.
        :type isomer: Any
        :return: Tuple (sorted_distances, sorted_element_pairs).
        :rtype: tuple[ np.ndarray, list[ tuple[str, str] ] ]
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
        Compare two fingerprint vectors with various metrics.

        Supported modes:
          - 'max_diff': maximum absolute element-wise difference
          - 'sum_diff': sum of absolute differences
          - 'mean_diff': mean absolute difference
          - 'rmsd': root mean square difference

        :param fp1: Fingerprint array for isomer 1.
        :type fp1: np.ndarray
        :param fp2: Fingerprint array for isomer 2.
        :type fp2: np.ndarray
        :param mode: Comparison mode string.
        :type mode: str
        :return: Scalar score representing fingerprint difference according to mode.
        :rtype: float
        :raises ValueError: If an unsupported mode is supplied.
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
        Visualize the fingerprint difference matrix as a heatmap (Matplotlib and optional Plotly).

        Produces an SVG file and optionally launches an interactive Plotly Dash app. The
        function will compute similarity grouping if not already available.

        :param write_svg: Whether to save a Matplotlib SVG of the heatmap.
        :type write_svg: bool
        :param plot_plotly: Whether to launch an interactive Plotly/Dash visualization.
        :type plot_plotly: bool
        :param colorscale_min: Lower bound for custom colorscale interpolation.
        :type colorscale_min: float
        :param colorscale_mid: Midpoint for colorscale interpolation.
        :type colorscale_mid: float
        :param color_scale_max: Upper bound for colorscale interpolation.
        :type color_scale_max: float
        :param min_color: Dictionary with 'r','g','b' ints for minimum color.
        :type min_color: dict | None
        :param mid_color: Dictionary with 'r','g','b' ints for mid color.
        :type mid_color: dict | None
        :param max_color: Dictionary with 'r','g','b' ints for maximum color.
        :type max_color: dict | None
        :param cell_label_mode: 'value' to show numeric differences, 'group' to show Close/Far, or 'none'.
        :type cell_label_mode: str
        :return: None
        :rtype: None
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
        Launch an interactive Dash application showing the difference matrix and histogram.

        Clicking on a heatmap cell opens the ASE viewer with the corresponding isomer alignment.

        :param df: Pandas DataFrame representing the difference matrix (square).
        :type df: pd.DataFrame
        :param group_labels_matrix: Square matrix of 'Close'/'Far' labels.
        :type group_labels_matrix: np.ndarray
        :param color_scale: Plotly-compatible color scale definition.
        :type color_scale: list
        :param cell_label_mode: Mode for cell labels ('value', 'group', 'none').
        :type cell_label_mode: str
        :return: None
        :rtype: None
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
        Visualize two isomers and their optimal alignment using ASE viewer.

        Displays three frames: reference isomer (index1), unaligned isomer (index2),
        and overlaid result after aligning isomer2 to isomer1.

        :param index1: Index of the reference isomer in self.isomers.
        :type index1: int
        :param index2: Index of the isomer to align and overlay.
        :type index2: int
        :param grid_size: Optional grid density for brute alignment; defaults to object's grid_size.
        :type grid_size: int | None
        :return: None
        :rtype: None
        :raises AssertionError: If indices are out of range.
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
        Display and return a DataFrame comparing two isomer fingerprints for debugging.

        The DataFrame contains distances, element pairs and absolute differences to help
        interpret which entries contribute most to the fingerprint score.

        :param idx1: Index of the first isomer to compare.
        :type idx1: int
        :param idx2: Index of the second isomer to compare.
        :type idx2: int
        :return: Pandas DataFrame summarizing distances, pairs and differences.
        :rtype: pd.DataFrame
        :raises ValueError: If fingerprints differ in length.
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
    Filter assembled isomers for interatomic clashes based on covalent radii.

    Clashes are determined by comparing pairwise distances to the sum of covalent radii plus a buffer.
    Pairs within the same ligand are ignored. Metal-involving pairs can optionally be excluded from checks.
    """

    def __init__(
            self,
            buffer: float = -0.3,
            check_metal_clashes: bool = False
    ):
        """
        Initialize clash filter configuration.

        :param buffer: Distance buffer (Å) added to sum of covalent radii. Negative allows closer approach.
        :type buffer: float
        :param check_metal_clashes: If False, pairs involving any metal atom are ignored.
        :type check_metal_clashes: bool
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
        Determine whether the provided ASE Atoms contains any clashing atom pairs.

        Ignores intra-ligand pairs and optionally ignores any pair involving metal atoms.
        Clashing is True if any remaining interatomic distance is smaller than the corresponding
        sum of covalent radii plus the configured buffer.

        :param atoms: ASE Atoms object representing the assembled isomer.
        :type atoms: ase.Atoms
        :param ligand_idc: List of lists of atom indices corresponding to each ligand.
        :type ligand_idc: list[list[int]]
        :param metal_idc: List of indices of metal center atoms within `atoms`.
        :type metal_idc: list[int]
        :return: True if any non-ignored pair is closer than allowed, otherwise False.
        :rtype: bool
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

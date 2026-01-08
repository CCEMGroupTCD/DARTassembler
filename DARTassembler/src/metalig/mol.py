"""
The MetaLig contains ligands, which are represented by the :class:`Ligand` class.
"""

import warnings
from functools import cached_property

import ase
import networkx as nx
import numpy as np
from copy import deepcopy
from typing import Union
from ase.visualize import view
from ase import Atoms
from DARTassembler.src.misc.refactor_v1_0_0 import refactor_metalig_entry_from_v1_0_0_to_v1_1_0
from DARTassembler.src.metalig.utils_molecule import get_planarity, get_denticities_and_hapticities_idc, \
    get_isomers_effective_ligand_atoms_with_effective_donor_indices, \
    get_all_effective_ligand_atoms_with_effective_donor_indices, format_hapdent_idc, has_smarts_pattern, \
    check_metal_center_format, get_atomic_props_from_ase_atoms, get_ase_atoms_from_atomic_props, stoichiometry2atomslist
from DARTassembler.src.constants.chem import Element
from DARTassembler.src.metalig.utils_graph import graph_from_graph_dict, graph_to_dict_with_node_labels, view_graph, \
    get_sorted_atoms_and_indices_from_graph, get_reindexed_graph, get_graph_fragments, count_atoms_with_n_bonds, get_graph_hash, get_heavy_atoms_graph, \
    get_adjacency_matrix, assert_graph_and_coordinates_are_consistent, \
    remove_node_features_from_graph, make_multigraph_to_graph, graph_to_smiles
from DARTassembler.src.metalig.utils import get_stable_sorted_value_counts, check_equal, sort_dict_as
from DARTassembler.src.metalig.utils_molecule import get_standardized_stoichiometry_from_atoms_list, \
    unknown_rdkit_bond_orders, get_max_deviation_from_coplanarity, if_same_stoichiometries
from DARTassembler.src.metalig.archetype import assign_archetype

pseudo_metal = 'Cu'     # Pseudo metal for display in ligand xyz files and for use in the SMARTS filter.
# The set of properties that are stored in the ligand dictionary of the MetaLig database.
ligand_dict_props = {'atomic_props', 'global_props', 'graph', 'donor_idc', 'ligand_instances', 'hapdent_idc', 'geometric_isomers_hapdent_idc'}
# The set of properties in the ligand_instances dictionary of each ligand in the MetaLig database.
ligand_instance_props = {'ligand_name', 'parent_complex_id', 'parent_complex_charge', 'parent_metal', 'parent_metal_os'}
# The properties expected in the ligand `global_props`. Must be a list instead of a set, because the order matters for the csv output.
ligand_global_props = [
    # General properties
    'unique_name',
    'stoichiometry',
    'archetype',
    'charge',
    'smiles',
    'smiles_with_metal',
    # Integer numerical properties
    'n_donors',
    'n_eff_denticities',
    'n_denticities',
    'n_haptic_atoms',
    'n_haptic_groups',
    'n_atoms',
    'n_elements',
    'n_bonds',
    'n_electrons',
    'n_protons',
    'n_beta_hydrogens',
    # Float numerical properties
    'molecular_weight',
    'planarity',
    'donor_planarity',
    'donor_metal_planarity',
    'min_interatomic_distance',
    'max_ligand_extension',
    'archetype_rssd',
    'archetype_confidence',
    # Boolean properties
    'is_2D_symmetrical',
    'has_all_bond_orders_valid',
    'has_confident_charge',
    # Graph hashes
    'graph_hash',
    'graph_hash_with_metal',
    'heavy_atoms_graph_hash',
    'heavy_atoms_graph_hash_with_metal',
    'bond_order_graph_hash',
    # Parent complex properties from CSD
    'n_ligand_instances',
]

class _cached_global_props(cached_property):
    """
    Cache a computed property and mirror it into the instance global_props when applicable.

    This descriptor behaves like functools.cached_property but additionally stores the computed
    value in instance.global_props if the property name is present in the canonical ligand_global_props
    list. This ensures derived global properties are persisted for serialization.

    The descriptor returns the cached value on subsequent accesses.
    """
    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        name = self.func.__name__
        # 1. Check if the value is already in `global_props`
        if name in instance.global_props:
            return instance.global_props[name]
        # 2. Otherwise, delegate to cached_property’s __get__, which will compute (if needed) and stash in __dict__
        value = super().__get__(instance, owner)
        # 3. Store the computed value in `global_props` as well
        if name in ligand_global_props:
            instance.global_props[name] = value
        return value

class BaseMolecule(object):

    def __init__(self,
                 graph: nx.Graph,
                 atomic_props: Union[dict, ase.Atoms],
                 global_props: dict = None,
                 validity_check = True,
                 node_label: str = 'node_label',
                 bond_label: str = 'bond_type',
                 ):
        """
        Initialize a :class:`BaseMolecule` instance, which represents a molecule with molecular graph (atomic connectivity) and atomic properties such as 3D coordinates.

        :param nx.Graph graph: NetworkX Graph describing connectivity. Nodes should correspond to atomic_props ordering.
        :param atomic_props: ASE Atoms object or atomic properties dictionary. Example format:

            .. code-block:: python

                {
                    'atoms': ['C','H','O',...],
                    'x': [x1, x2, ...],
                    'y': [y1, y2, ...],
                    'z': [z1, z2, ...],
                    ... (optional per-atom arrays)
                }

        :type atomic_props: dict | ase.Atoms
        :param global_props: Optional flat dictionary of global properties (charge, stoichiometry, etc.). Example format:

            .. code-block:: python

                {
                    'charge': -1,
                    'stoichiometry': 'C6H6O',
                    ...
                }

        :type global_props: dict | None
        :param validity_check: If True, perform internal consistency checks comparing graph and coordinates.
        :type validity_check: bool
        :param node_label: Node attribute key used for element symbols in graphs.
        :type node_label: str
        :param bond_label: Edge attribute key used for bond type information in graphs.
        :type bond_label: str
        :raises AssertionError: If graph and atomic_props are inconsistent when validity_check is True.
        """
        if global_props is None:
            global_props = {}

        # Convert ase.Atoms to `atomic_props` dictionary if needed
        try:
            atomic_props = get_atomic_props_from_ase_atoms(atomic_props)
        except AttributeError:
            pass

        self.atomic_props = atomic_props
        self.global_props = global_props
        self.node_label = node_label  # node label in the graph specifying the atom type
        self.bond_label = bond_label  # edge label in the graph specifying the bond type
        self.graph = get_reindexed_graph(graph)
        self.atoms = self.get_ase_atoms()

        if validity_check:
            self._check_if_molecule_valid()

    @_cached_global_props
    def n_atoms(self) -> int:
        """
        Return the number of atoms in the molecule.

        :return: Atom count.
        :rtype: int
        """
        return len(self.atomic_props['atoms'])

    @_cached_global_props
    def n_hydrogens(self) -> int:
        """
        Return the number of hydrogen atoms in the molecule.

        :return: Count of hydrogen atoms.
        :rtype: int
        """
        return sum(1 for atom in self.atomic_props['atoms'] if atom == 'H')

    @_cached_global_props
    def n_protons(self) -> int:
        """
        Return the total number of protons (sum of atomic numbers) present in the molecule.

        :return: Sum of atomic numbers across all atoms.
        :rtype: int
        """
        return sum([Element(el).atomic_number for el in self.atomic_props['atoms']])

    @_cached_global_props
    def n_bonds(self) -> int:
        """
        Return the number of bonds as given by the molecular graph.

        :return: Number of edges in the molecular graph.
        :rtype: int
        """
        return len(self.graph.edges)

    @_cached_global_props
    def has_bond_order_attribute(self) -> bool:
        """
        Indicate whether graph edges expose the configured bond order attribute.

        :return: True when every edge contains the configured bond_label.
        :rtype: bool
        """
        return all(self.bond_label in self.graph.edges[edge] for edge in self.graph.edges)

    @_cached_global_props
    def has_unknown_bond_orders(self) -> bool:
        """
        Indicate whether the graph contains unknown or placeholder RDKit bond orders.

        If the bond order attribute is missing, the property returns True (unknown).
        :return: True when any edge has an unknown bond order or bond orders are missing/invalid.
        :rtype: bool
        """
        if not self.has_bond_order_attribute:
            return True
        return self.count_bond_types(unknown_rdkit_bond_orders) > 0

    @_cached_global_props
    def has_all_bond_orders_valid(self) -> bool:
        """
        Indicate whether all bond orders in the graph are present and valid.

        :return: True if bond orders exist and none are unknown.
        :rtype: bool
        """
        return self.has_bond_order_attribute and not self.has_unknown_bond_orders

    @_cached_global_props
    def graph_hash(self) -> str:
        """
        Compute the Weisfeiler–Lehman graph hash for the molecular graph (element labels only).

        :return: String fingerprint of the graph topology.
        :rtype: str
        """
        return get_graph_hash(self.graph)

    @_cached_global_props
    def heavy_atoms_graph_hash(self) -> str:
        """
        Compute the WL graph hash after removing hydrogen atoms.

        :return: Hash string for heavy-atom-only graph.
        :rtype: str
        """
        heavy_graph = get_heavy_atoms_graph(self.graph, element_label=self.node_label)
        return get_graph_hash(heavy_graph)

    @_cached_global_props
    def bond_order_graph_hash(self) -> Union[str, None]:
        """
        Compute a graph hash that includes bond order information when available.

        :return: Hash string including bond orders, or None if bond orders are missing/invalid.
        :rtype: str | None
        """
        if not self.has_all_bond_orders_valid:
            return None
        return get_graph_hash(self.graph, node_attr=self.node_label, edge_attr=self.bond_label)

    @_cached_global_props
    def stoichiometry(self) -> str:
        """
        Return the canonical stoichiometry string for the molecule (e.g. 'C6H6').

        :return: Standardized stoichiometry string.
        :rtype: str
        """
        return get_standardized_stoichiometry_from_atoms_list(self.atomic_props['atoms'])

    def get_reindexed_graph(self) -> nx.Graph:
        """
        Return a reindexed copy of the internal molecular graph (0..n-1 node indices).

        :return: NetworkX Graph with contiguous node indices.
        :rtype: nx.Graph
        """
        return get_reindexed_graph(self.graph)

    def get_ase_atoms(self, remove_elements: list=None, add_atoms: list=None) -> Atoms:
        """
        Produce an ASE Atoms object from internal atomic_props with optional filtering and additions.

        :param remove_elements: List of element symbols to remove from the returned ASE Atoms (e.g. ['H']).
        :type remove_elements: list[str] | None
        :param add_atoms: List of tuples (element_symbol, (x, y, z)) specifying atoms to append.
        :type add_atoms: list[ tuple[str, tuple[float, float, float]] ] | None
        :return: ASE Atoms object representing the molecule with requested modifications.
        :rtype: ase.Atoms
        """
        return get_ase_atoms_from_atomic_props(atomic_props=self.atomic_props, remove_elements=remove_elements,add_atoms=add_atoms)

    def get_smiles(self) -> Union[str,None]:
        """
        Return a SMILES string of the molecular graph if bond orders are available.

        :return: SMILES string or None when bond orders are unknown/invalid.
        :rtype: str | None
        """
        if not self.has_all_bond_orders_valid: # if the molecule has unknown bond orders, we cannot calculate the SMILES
            return None

        smiles = graph_to_smiles(self.graph, element_label=self.node_label, bond_label=self.bond_label)

        return smiles

    def count_C_H_bonds(self) -> int:
        """
        Count the number of C–H bonds in the molecular graph.

        Iterates over edges and counts edges between 'C' and 'H'.

        :return: Number of C–H bonds.
        :rtype: int
        """
        n_bonds = 0
        atoms = self.graph.nodes(data=True)
        for idx1, idx2 in self.graph.edges():
            el1 = atoms[idx1][self.node_label]
            el2 = atoms[idx2][self.node_label]
            if sorted([el1, el2]) == ['C', 'H']:
                n_bonds += 1

        return n_bonds

    def count_bond_types(self, bond_types: list) -> Union[int, float]:
        """
        Count occurrences of specified RDKit bond-type codes in the molecular graph.

        :param bond_types: Iterable of RDKit bond-type integers to count (e.g. aromatic, single, double codes).
        :type bond_types: list[int]
        :return: Number of edges whose bond_type attribute matches any entry in bond_types.
        :rtype: int
        :raises ValueError: If the graph lacks the configured bond order attribute.
        """
        if not self.has_bond_order_attribute:
            raise ValueError("The molecular graph does not have a bond order attribute.")

        n = 0
        for idx1, idx2, bond_dict in self.graph.edges(data=True):
            bond_type = bond_dict[self.bond_label]
            if bond_type in bond_types:
                n += 1

        return n

    def get_graph_fragments(self) -> tuple[list, list]:
        """
        Return connected component indices and corresponding element-lists for each fragment.

        :return: Tuple (fragment_indices, fragment_elements) where each entry corresponds to a disconnected component.
        :rtype: tuple[list, list]
        """
        return get_graph_fragments(graph=self.graph, node_label=self.node_label)

    def view_3D(self) -> None:
        """
        Open an ASE-based 3D viewer for the current ASE Atoms representation.

        :return: None
        :rtype: None
        """
        view(self.atoms)

    def view_graph(self, node_size=150) -> None:
        """
        Visualize the molecular graph using the project helper.

        :param node_size: Display size for graph nodes in the viewer.
        :type node_size: int
        :return: None
        :rtype: None
        """
        view_graph(self.graph, node_label=self.node_label, node_size=node_size)

    def get_coordinates(self) -> list:
        """
        Return the atomic coordinates as a list of 3‑tuples.

        :return: Coordinates in format [[x1, y1, z1], [x2, y2, z2], ...].
        :rtype: list[list[float]]
        """
        return [[self.atomic_props['x'][i], self.atomic_props['y'][i], self.atomic_props['z'][i]] for i in range(len(self.atomic_props['x']))]

    def get_interatomic_distances(self, skip_elements=None) -> tuple[float, float, np.ndarray]:
        """
        Compute pairwise interatomic distances and return minimum, maximum and full distance matrix.

        :param skip_elements: Element symbol or list of symbols to exclude from distance calculations.
        :type skip_elements: str | list[str] | None
        :return: Tuple (min_distance, max_distance, distances_matrix).
        :rtype: tuple[ float, float, np.ndarray ]
        """
        if skip_elements is None:
            skip_elements = []
        elif isinstance(skip_elements, str):
            skip_elements = [skip_elements]

        atoms = self.atoms.get_chemical_symbols()
        valid = np.array([not (el in skip_elements) for el in atoms])

        distances = self.atoms.get_all_distances()
        distances = distances[valid,:][:,valid]

        if len(distances) > 1:
            min_dist = np.where(distances > 0, distances, np.inf).min()
            max_dist = np.where(distances > 0, distances, -np.inf).max()
        else:
            min_dist = 0
            max_dist = 0

        return min_dist, max_dist, distances

    def get_all_interatomic_distances_flat(self) -> list:
        """
        Return unique pairwise interatomic distances as a flattened list (upper triangle).

        :return: List of pairwise distances (each pair reported once).
        :rtype: list[float]
        """
        distances = self.get_interatomic_distances()[2]
        unique_distances = []
        for i in range(len(distances)):
            for j in range(i+1, len(distances)):
                unique_distances.append(distances[i,j])

        return unique_distances

    def get_xyz_string(self, comment: str = '') -> str:
        """
        Produce a string in XYZ file format for the molecule.

        :param str comment: Optional single-line comment to include in the XYZ header.
        :return: Multiline string conforming to the XYZ format.
        :rtype: str
        """
        from io import StringIO
        from ase.io import write

        atoms = get_ase_atoms_from_atomic_props(self.atomic_props)
        buf = StringIO()
        write(buf, atoms, format='xyz', comment=comment)
        xyz = buf.getvalue().rstrip('\n')

        return xyz

    def to_dict(self, include_graph: bool=True) -> dict:
        """
        Serialize the molecule to a dictionary suitable for JSON output.

        :param include_graph: If True, include graph representation with node labels.
        :type include_graph: bool
        :return: Dictionary with keys 'atomic_props', 'global_props' and optionally 'graph'.
        :rtype: dict
        """
        d = {
                'atomic_props': self.atomic_props,
                'global_props': self.global_props,
        }
        if include_graph:
            d['graph'] = graph_to_dict_with_node_labels(self.graph)

        return d

    def _append_to_file(self, key: str, writer) -> None:
        """
        Append the serialized molecule to a writer (file-like) under a key.

        :param key: Identifier under which the serialized object will be written.
        :type key: str
        :param writer: Writer object implementing a write(dict) method.
        :type writer: object
        :return: None
        :rtype: None
        """
        data = {'key': key, 'value': self.to_dict()}
        writer.write(data)

        return

    def _count_atoms_with_n_bonds(self, element: Union[str, None], n_bonds: int) -> int:
        """
        Count atoms of a given element that have exactly n_bonds connections.

        :param element: Element symbol to filter by (e.g. 'C'). Use None to count all elements.
        :type element: str | None
        :param n_bonds: Number of bonds that qualify an atom to be counted.
        :type n_bonds: int
        :return: Integer count of atoms matching the criteria.
        :rtype: int
        """
        return count_atoms_with_n_bonds(graph=self.graph, element=element, n_bonds=n_bonds, graph_element_label=self.node_label)

    def _contains_only(self, atoms: Union[str, list], except_elements=None) -> bool:
        """
        Return whether the molecule contains only elements from a given set (optionally excluding some).

        :param atoms: Element symbol or list of symbols allowed.
        :type atoms: str | list[str]
        :param except_elements: Elements to ignore when testing containment.
        :type except_elements: list[str] | None
        :return: True if every non-excluded atom is in `atoms`, False otherwise.
        :rtype: bool
        """
        if except_elements is None:
            except_elements = []
        if isinstance(atoms, str):
            atoms = [atoms]

        own_atoms = [atom for atom in self.atomic_props['atoms'] if not atom in except_elements]
        contains_only_atoms = all(np.isin(own_atoms, atoms))

        return contains_only_atoms

    def _check_if_molecule_valid(self) -> None:
        """
        Validate internal consistency between graph and coordinate/atomic lists.

        Raises an AssertionError if atom lists or ordering mismatch, or if graph/coordinates are inconsistent.

        :return: None
        :rtype: None
        :raises AssertionError: If atoms in graph and atomic_props mismatch or adjacency assertions fail.
        """
        atoms, _ = get_sorted_atoms_and_indices_from_graph(self.graph)

        both_atom_lists_printed = f'\ngraph atoms: {atoms}\natomic_props["atoms"]: {self.atomic_props["atoms"]}'
        same_atoms_contained = sorted(atoms) == sorted(self.atomic_props['atoms'])
        assert same_atoms_contained, f'The atoms from the graph and the atoms from the atomic_props don\'t match:{both_atom_lists_printed}'
        same_order_of_atoms = atoms == self.atomic_props['atoms']
        assert same_order_of_atoms, f'Order of atoms in graph and in atomic_props doesn\'t match:{both_atom_lists_printed}'

        try:    # Check if the molecule is a ligand
            donor_idc = self.donor_idc
        except AttributeError:  # Molecule is not a ligand
            donor_idc = None
        assert_graph_and_coordinates_are_consistent(
                                                        graph=self.graph,
                                                        graph_hash=self.graph_hash,
                                                        atoms=self.atomic_props['atoms'],
                                                        donor_idc=donor_idc,
                                                        node_label=self.node_label
                                                    )

        return

    def _if_same_stoichiometry(self, other: str) -> bool:
        """
        Compare this molecule's stoichiometry against another stoichiometry string.

        :param other: Stoichiometry string to compare with (e.g. 'C6H6').
        :type other: str
        :return: True when stoichiometries are identical, False otherwise.
        :rtype: bool
        """
        return if_same_stoichiometries(self.stoichiometry, other)

    def _get_rdkit_mol_from_smiles(self, sanitize: bool=False):
        """
        Construct an RDKit Mol object from the internal SMILES string.

        :param sanitize: If True, run RDKit sanitization on the molecule.
        :type sanitize: bool
        :return: RDKit Chem.Mol object or None if SMILES cannot be computed.
        :rtype: rdkit.Chem.Mol | None
        :raises ImportError: If RDKit is not installed.
        """
        from rdkit import Chem
        smiles = self.get_smiles()
        if smiles is None:
            return None

        return Chem.MolFromSmiles(smiles, sanitize=sanitize)

    def _get_planarity(self) -> float:
        """
        Compute a unit-scaled planarity score based on maximum deviation from a best-fit plane.

        The returned score is in (0, 1], where 1.0 indicates perfect coplanarity.

        :return: Planarity score between 0 and 1 (rounded to 10 decimals).
        :rtype: float
        """
        coordinates = self.get_coordinates()
        deviation = get_max_deviation_from_coplanarity(points=coordinates)  # deviation is a float that is 0 if the molecule is perfectly planar and > 0 if it is not. The higher the value, the less planar the molecule is.
        planarity = 1/ (1+ deviation)   # planarity is a float between 0 and 1. 0 means not planar at all (a sphere), 1 means perfectly planar.
        planarity = round(planarity, 10)    # round to 10 decimal places to avoid floating point errors which happen with np.linalg.svd() in different versions of numpy

        return planarity

    def _get_xyz_array(self) -> np.ndarray:
        """
        Return atomic coordinates as a numpy array with shape (n_atoms, 3).

        :return: Coordinates array with rows corresponding to atoms.
        :rtype: np.ndarray
        """
        return np.array([self.atomic_props['x'], self.atomic_props['y'], self.atomic_props['z']]).T

    def _ligand_naming(self, n_donors: int, ligand_list: list) -> (str, str):
        """
        Generate a unique ligand name based on donor count and existing ligand list.

        :param n_donors: Number of donor atoms for the ligand.
        :type n_donors: int
        :param ligand_list: Sequence of existing ligand objects to ensure uniqueness.
        :type ligand_list: list[ object ]
        :return: Tuple (ligand_name, csd) where ligand_name is unique and csd is parent CSD code if present.
        :rtype: tuple[str, str]
        """
        try:
            lig_key = f"CSD-{self.global_props['parent_complex_id']}"
            csd = self.global_props["CSD_code"]
        except KeyError:
            lig_key = f"CSD-{self.global_props['parent_complex_id']}"
            csd = self.global_props["CSD_code"]

        from DARTassembler.src._extraction.constants import mini_alphabet
        j = 0
        while True:
            ligand_name = f'{lig_key}-0{n_donors}-{mini_alphabet[j]}'
            if ligand_name not in [lig.name for lig in ligand_list]:
                break
            else:
                j += 1

        return ligand_name, csd

    def _check_input_inherit_global_properties(self, inherit_global_properties: list) -> list:
        """
        Validate and normalize an 'inherit_global_properties' specification.

        If None is provided, returns all keys from self.global_props. Otherwise ensure all requested
        properties exist in self.global_props.

        :param inherit_global_properties: List of property names to inherit or None.
        :type inherit_global_properties: list[str] | None
        :return: Validated list of property names to inherit.
        :rtype: list[str]
        :raises ValueError: If any requested property is missing from self.global_props.
        """
        if inherit_global_properties is None:
            inherit_global_properties = list(self.global_props.keys())
        else:
            unknown_global_property = [prop for prop in inherit_global_properties if not prop in self.global_props]
            if unknown_global_property:
                raise ValueError(
                    f'Unknown values {unknown_global_property}. All properties in inherit_global_properties must be found in `self.global_properties`.')

        return inherit_global_properties

    def _remove_node_features_from_molecular_graphs_inplace(self, keep: list = None) -> None:
        """
        Remove node feature attributes from the internal graph except those explicitly kept.

        :param keep: List of node attribute names to retain. Defaults to [self.node_label].
        :type keep: list[str] | None
        :return: None
        :rtype: None
        """
        if keep is None:
            keep = [self.node_label]

        remove_node_features_from_graph(graph=self.graph, keep=keep, inplace=True)

        return

    def _normalize_multigraph_into_graph_inplace(self) -> None:
        """
        Convert a MultiGraph stored in self.graph into a standard Graph in-place.

        This collapses parallel edges into single edges while preserving node attributes.
        :return: None
        :rtype: None
        """
        self.graph = make_multigraph_to_graph(self.graph)

        return

class Ligand(BaseMolecule):

    def __init__(self,
                 atomic_props: Union[dict, ase.Atoms],
                 donor_idc: list[int],
                 graph: nx.Graph,
                 unique_name: str,
                 charge: Union[int, float],
                 global_props: dict = None,
                 ligand_instances: dict = None,
                 hapdent_idc: tuple = None,
                 geometric_isomers_hapdent_idc: list = None,
                 validity_check=False,
                 ):
        """
        Initialize a :class:`Ligand` instance representing a metal-coordinating ligand from the :ref:`MetaLig <metalig>` database.

        The constructor accepts atomic and graph data plus metadata describing donor indices,
        known parent complex instances and optional precomputed hapdent (denticity/hapticity) tuples.

        :param atomic_props: ASE Atoms object or atomic properties dictionary with keys 'atoms', 'x', 'y', 'z', etc.
        :type atomic_props: dict | ase.Atoms
        :param donor_idc: List of atomic indices that act as donor atoms coordinating to a metal.
        :type donor_idc: list[int]
        :param graph: NetworkX Graph describing connectivity (node labels must match atomic_props['atoms'] ordering).
        :type graph: nx.Graph
        :param unique_name: Unique identifier string for the ligand (database key).
        :type unique_name: str
        :param charge: Formal charge of the ligand. May be np.nan when unknown.
        :type charge: int | float
        :param global_props: Optional dictionary of global properties.
        :type global_props: dict | None
        :param ligand_instances: Optional dictionary describing parent complex instances. Example format:

            .. code-block:: python

                {
                    'ligand_name': ['lig1','lig2',...],
                    'parent_complex_id': ['CSD1','CSD2',...],
                    'parent_complex_charge': [0, +1, ...],
                    'parent_metal': ['Cu','Fe',...],
                    'parent_metal_os': [2, 2, ...],
                }

        :type ligand_instances: dict | None
        :param hapdent_idc: Optional denticity/hapticity tuple structure (haptic groups as sub-tuples).
        :type hapdent_idc: tuple | None
        :param geometric_isomers_hapdent_idc: Optional list of hapdent tuples for geometric isomers.
        :type geometric_isomers_hapdent_idc: list[ tuple[int | tuple[int]] ] | None
        :param validity_check: If True, run consistency checks after initialization.
        :type validity_check: bool
        :raises ValueError: If ligand_instances lacks required keys.
        """
        if global_props is None:
            global_props = {}
        if ligand_instances is None:
            ligand_instances = {
                'ligand_name': [],
                'parent_complex_id': [],
                'parent_complex_charge': [],
                'parent_metal': [],
                'parent_metal_os': [],
            }
        elif not ligand_instance_props.issubset(set(ligand_instances.keys())):  # Check if all necessary keys are present
            raise ValueError(f'The dictionary `ligand_instances` must contain the keys {ligand_instance_props}.')

        super().__init__(
                         atomic_props=atomic_props,
                         global_props=global_props,
                         graph=graph,
                         validity_check=False,  # Will be done all together at the end of the __init__ method.
                         )

        self.unique_name = unique_name
        self.other_ligand_instances = ligand_instances
        self.donor_idc = donor_idc
        self.charge = charge
        # Saving the hapdent_idc as json converts the tuples to lists, so we need to convert them back to tuples. If None, let them undefined for later calculation.
        if hapdent_idc is not None:
            self.hapdent_idc = format_hapdent_idc(hapdent_idc)
        if geometric_isomers_hapdent_idc is not None:
            self.geometric_isomers_hapdent_idc = [format_hapdent_idc(isomer_idc) for isomer_idc in geometric_isomers_hapdent_idc]

        if validity_check:
            self._check_if_molecule_valid()

    @_cached_global_props
    def has_confident_charge(self) -> bool:
        """
        Return whether the stored ligand charge is flagged as confident.

        If the flag is present in global_props, that value is returned; otherwise True is assumed.

        :return: Boolean indicator whether the charge is considered confident.
        :rtype: bool
        """
        if 'has_confident_charge' in self.global_props:
            return self.global_props['has_confident_charge']
        else:
            return True

    @_cached_global_props
    def n_ligand_instances(self) -> int:
        """
        Return the number of recorded parent complex instances for this ligand.

        :return: Number of parent occurrences recorded in ligand_instances.
        :rtype: int
        """
        return len(self.other_ligand_instances['ligand_name'])

    @_cached_global_props
    def donor_metal_planarity(self) -> float:
        """
        Return the planarity score of donor atoms including the original metal center.

        :return: Planarity (0..1) of donor atoms with the metal included.
        :rtype: float
        """
        return self._get_donors_planarity(with_metal=True)

    @_cached_global_props
    def donor_planarity(self) -> float:
        """
        Return the planarity score of donor atoms alone, excluding the metal.

        :return: Planarity (0..1) of donor atoms.
        :rtype: float
        """
        return self._get_donors_planarity(with_metal=False)

    @_cached_global_props
    def planarity(self) -> float:
        """
        Return the overall planarity score of the ligand.

        :return: Planarity (0..1) of the ligand.
        :rtype: float
        """
        return self._get_planarity()

    @_cached_global_props
    def n_donors(self):
        """
        Return the number of donor atoms defined for this ligand.

        :return: Integer number of donor atoms.
        :rtype: int
        """
        return len(self.donor_idc)

    @_cached_global_props
    def donor_elements(self) -> list[str]:
        """
        Return the element symbols of the donor atoms in ligand order.

        :return: List of element symbols for donor atoms (length n_donors).
        :rtype: list[str]
        """
        return [self.atomic_props['atoms'][i] for i in self.donor_idc]

    @_cached_global_props
    def donor_positions(self) -> np.ndarray:
        """
        Return the 3D coordinates of donor atoms as an (n_donors, 3) array.

        :return: numpy array of donor coordinates.
        :rtype: np.ndarray
        """
        return np.array([[self.atomic_props["x"][i], self.atomic_props["y"][i], self.atomic_props["z"][i]] for i in self.donor_idc])

    @_cached_global_props
    def parent_complex_id(self) -> str:
        """
        Return the identifier of the first parent complex instance for this ligand.

        :return: Parent complex ID string.
        :rtype: str
        """
        return self.other_ligand_instances['parent_complex_id'][0]

    @_cached_global_props
    def parent_metal_position(self) -> np.ndarray:
        """
        Return the 3D coordinates of the parent metal center observed in the source complex.

        :return: 3‑vector coordinates of the parent metal.
        :rtype: np.ndarray
        """
        return self.other_ligand_instances['parent_metal_position'][0]

    @_cached_global_props
    def parent_metal(self) -> str:
        """
        Return the chemical symbol of the parent metal observed for this ligand.

        :return: Element symbol (e.g. 'Cu', 'Fe').
        :rtype: str
        """
        return self.other_ligand_instances['parent_metal'][0]

    @_cached_global_props
    def hapdent_idc(self) -> tuple:
        """
        Return the donor index tuple structure with haptic groups as sub-tuples.

        :return: Tuple combining integer donor indices and sub-tuples for haptic groups.
        :rtype: tuple
        """
        return self._get_denticities_and_hapticities_idc()

    @_cached_global_props
    def n_eff_denticities(self) -> int:
        """
        Return the effective denticity counting each haptic group as a single donor.

        :return: Effective denticity integer.
        :rtype: int
        """
        return len(self.hapdent_idc)

    @_cached_global_props
    def n_denticities(self) -> int:
        """
        Return the classical denticity (count of integer donor entries, ignoring haptic groups).

        :return: Integer denticity count.
        :rtype: int
        """
        return sum([1 for el in self.hapdent_idc if isinstance(el, int)])

    @_cached_global_props
    def n_haptic_atoms(self) -> int:
        """
        Return the total number of atoms that participate in haptic coordination.

        :return: Integer count of atoms that are part of haptic groups.
        :rtype: int
        """
        return sum([len(sublist) for sublist in self.hapdent_idc if isinstance(sublist, tuple)])

    @_cached_global_props
    def n_haptic_groups(self) -> int:
        """
        Return the number of haptic groups (hapticity units) present in the ligand.

        :return: Integer number of haptic groups.
        :rtype: int
        """
        return len([sublist for sublist in self.hapdent_idc if isinstance(sublist, tuple)])

    @_cached_global_props
    def _archetype_and_geometrical_isomers(self):
        """Cache all the archetype and isomer information when required."""
        archetype, _, isomer_hapdent_idc, rssd, _, archetype_confidence = self.get_ligand_archetype_and_isomers()
        d = {
            'archetype': archetype,                                                       # str
            'geometric_isomers_hapdent_idc': isomer_hapdent_idc,                        # list of hapdent_idc
            'archetype_rssd': rssd,                                                      # float >= 0.0
            'archetype_confidence': archetype_confidence                                  # float > 1.0
        }

        # Add to global_props so that the information is saved when the ligand is written to a file
        self.global_props['archetype'] = archetype
        self.global_props['archetype_rssd'] = rssd
        self.global_props['archetype_confidence'] = archetype_confidence

        return d

    @_cached_global_props
    def archetype(self) -> str:
        """
        Return the best-matching ligand archetype string (e.g. '2-cis').

        :return: Archetype identifier.
        :rtype: str
        """
        return self._archetype_and_geometrical_isomers['archetype']

    @_cached_global_props
    def geometric_isomers_hapdent_idc(self) -> list:
        """
        Return hapdent tuples for each recognised geometric isomer of the ligand.

        :return: List of hapdent tuples (each tuple contains int or sub-tuple entries).
        :rtype: list
        """
        return self._archetype_and_geometrical_isomers['geometric_isomers_hapdent_idc']

    @_cached_global_props
    def archetype_rssd(self) -> float:
        """
        Return the root-sum-of-squared-differences (RSSD) metric of archetype assignment.

        Lower values indicate closer match to the archetype ideal geometry.
        :return: RSSD value (float).
        :rtype: float
        """
        return self._archetype_and_geometrical_isomers['archetype_rssd']

    @_cached_global_props
    def archetype_confidence(self) -> float:
        """
        Return a confidence metric for the archetype assignment.

        Higher values indicate greater confidence in the assigned archetype.
        :return: Archetype confidence float.
        :rtype: float
        """
        return self._archetype_and_geometrical_isomers['archetype_confidence']

    @_cached_global_props
    def min_interatomic_distance(self) -> float:
        """
        Return the minimum interatomic distance observed in the ligand (Å).

        :return: Minimum positive pairwise distance.
        :rtype: float
        """
        return self.get_interatomic_distances()[0]

    @_cached_global_props
    def max_ligand_extension(self) -> float:
        """
        Return the maximal interatomic distance (molecular extension) in the ligand (Å).

        :return: Maximum pairwise distance.
        :rtype: float
        """
        return self.get_interatomic_distances()[1]

    @_cached_global_props
    def smiles(self) -> Union[str,None]:
        """
        Return the ligand SMILES string if bond orders are present.

        :return: SMILES string or None.
        :rtype: str | None
        """
        return self.get_smiles()

    @_cached_global_props
    def metal_counts(self) -> dict:
        """
        Return a dict mapping observed parent metal symbols to their counts.

        :return: Ordered dict-like mapping 'Element' -> count observed in parent complexes.
        :rtype: dict
        """
        return get_stable_sorted_value_counts(self.other_ligand_instances['parent_metal'])

    @_cached_global_props
    def metal_os_counts(self) -> dict:
        """
        Return parent metal oxidation-state counts as formatted strings.

        :return: Mapping like 'Fe+2' -> count.
        :rtype: dict
        """
        mos_counts = [f'{el}+{mos:.0f}' if mos > 0 else f'{el}{mos:.0f}' for el, mos in zip(self.other_ligand_instances['parent_metal'], self.other_ligand_instances['parent_metal_os']) if not np.isnan(mos)]
        return get_stable_sorted_value_counts(mos_counts)

    @_cached_global_props
    def smiles_with_metal(self) -> Union[str,None]:
        """
        Return a SMILES string of the ligand with a pseudo-metal appended (Hg) for visualization/analysis.

        :return: SMILES string including a pseudo-metal or None.
        :rtype: str | None
        """
        return self.get_smiles(with_metal='Hg')

    @_cached_global_props
    def is_2D_symmetrical(self) -> bool:
        """
        Return whether the ligand graph is topologically symmetrical between any pair of donor atoms.

        This is a 2D graph symmetry check and does not guarantee 3D symmetry.

        :return: True if any donor-pair symmetry is detected, False otherwise.
        :rtype: bool
        """
        return self._check_if_2D_symmetrical()

    @_cached_global_props
    def graph_hash_with_metal(self) -> str:
        """
        Return the WL graph hash of the ligand graph after connecting a pseudo-metal.

        :return: Hash string of ligand+metal graph.
        :rtype: str
        """
        return self.get_graph_hash_with_metal(metal_symbol='Hg')

    @_cached_global_props
    def heavy_atoms_graph_hash_with_metal(self) -> str:
        """
        Return the heavy-atom-only graph hash for the ligand with a pseudo-metal attached.

        :return: Hash string for heavy-atom-only ligand+metal graph.
        :rtype: str
        """
        return self.get_heavy_atoms_graph_hash_with_metal(metal_symbol='Hg')

    @_cached_global_props
    def n_beta_hydrogens(self) -> int:
        """
        Return the count of beta-hydrogen atoms (two bonds from coordinating atom, excluding alpha-H).

        :return: Integer number of beta hydrogens.
        :rtype: int
        """
        return self._get_n_beta_hydrogens()

    @_cached_global_props
    def is_haptic(self) -> bool:
        """
        Return whether the ligand contains any haptic coordination atoms.

        :return: True when haptic atoms/groups exist, False otherwise.
        :rtype: bool
        """
        return self.n_haptic_atoms > 0

    def get_smiles(self, with_metal: str=None) -> Union[str,None]:
        """
        Return the SMILES string for the ligand, optionally including a specified metal node.

        :param with_metal: Element symbol of the metal to attach to the ligand graph. If None, do not attach a metal.
        :type with_metal: str | None
        :return: SMILES string or None if bond orders are unknown.
        :rtype: str | None
        :raises ValueError: If with_metal is provided and is not a metal element.
        """
        # If the molecule has unknown bond orders, we cannot calculate the SMILES.
        if not self.has_all_bond_orders_valid:
            return None

        if with_metal is None:
            graph = self.graph
        else:
            if not Element(with_metal).is_metal:
                raise ValueError(f'Invalid input for with_metal: {with_metal}. Must be a metal symbol such as Fe.')

            graph, _ = self.get_graph_with_metal(metal_symbol=with_metal)

        smiles = graph_to_smiles(graph, element_label=self.node_label, bond_label=self.bond_label)

        return smiles
    
    def get_ase_atoms_with_metal(self, metal: str=None) -> Atoms:
        """
        Return an ASE Atoms object of the ligand with the parent metal placed at the recorded metal position.

        :param metal: Element symbol of the metal to add; if None the default parent metal from metadata is used.
        :type metal: str | None
        :return: ASE Atoms with the metal atom appended at parent_metal_position.
        :rtype: ase.Atoms
        """
        return self.get_ase_atoms(add_atoms=[(metal, self.parent_metal_position)])

    def property_filter(self, name, range: list = None, values: list = None) -> bool:
        """
        Test whether a certain simple property falls within a numeric range or matches one of the provided values.

        :param str name: The property name of one of the :ref:`MetaLig properties <metalig_properties>`.
        :param list[ tuple[float, float] ] | None range: Single range tuple or list of (min, max) tuples specifying allowed numerical intervals.
        :param list | None values: List of allowed exact values or None to skip value matching.
        :return: True if property value satisfies any provided test, otherwise False.
        :rtype: bool
        :raises ValueError: If the requested property is not present or a non-numeric property is tested against numeric ranges.

        **Examples :**

        To test whether a ligand has a '2-cis' archetype:

        .. code-block:: python

            ligand.property_filter(name='archetype', values=['2-cis'])

        To test whether a ligand has between 1 and 5 haptic atoms:

        .. code-block:: python

            ligand.property_filter(name='n_haptic_atoms', range=[1,5])

        To test whether a ligand has between 10 and 20 atoms or between 30 and 40 atoms:

        .. code-block:: python

            ligand.property_filter(name='n_atoms', range=[(10,20),(30,40)])
        """
        try:
            value = self.global_props[name]
        except KeyError:
            raise ValueError(f'Property {name} is not in `global_props`. Available properties are: {list(self.global_props.keys())}.')

        # If the value is specified in the list return always True. Also respect None in the list.
        value_in_list = values is not None and any([check_equal(value, val) for val in values])
        if value_in_list:
            return True

        # If the value is not specified in the list, check if it is in the range. Several ranges can be specified.
        if range is not None:
            if not isinstance(value, (int, float)):
                raise ValueError(f'Property {name} is not numerical, but {value}. Please do not specify a range for non-numerical properties.')
            if not isinstance(range[0], (tuple, list)):
                range = [range]

            # Check if all ranges are valid, i.e. numerical and of length 2.
            for r in range:
                correct_length = len(r) == 2
                if not correct_length:
                    raise ValueError(f'Ranges must be specified as tuples of two values (min & max), but got {len(r)} values: {r}.')
                numerical = all([isinstance(val, (int, float)) for val in r])
                if not numerical:
                    raise ValueError(f'Ranges must be specified as numerical values, but got {r}.')

            for min_, max_ in range:
                if min_ <= value <= max_:
                    return True

        return False

    def composition_filter(self, elements: Union[str,list[str]], instruction: str, only_donors: bool=False) -> bool:
        """
        Test whether a ligand satisfies a specified chemical composition criterion.

        Depending on the value of :arg:`instruction`, this method applies one of several
        composition-based filters using the provided :arg:`elements`. If :arg:`only_donors`
        is ``True``, only the donor atoms are considered; otherwise all atoms in the
        ligand are included.

        :param str | list[str] elements: Stoichiometry (e.g., ``'H2O'``) or list of element symbols (e.g., ``['H', 'H', 'O']``) used to define the expected atomic composition.
        :param str instruction: Specifies the composition rule to apply. Supported options are:

            - ``must_contain_and_only_contain`` – The ligand must consist of exactly these atoms in exactly this count.
              Use this to match an exact stoichiometry (e.g., ``C6H6`` for benzene).
            - ``must_at_least_contain`` – The ligand must contain all specified elements, but may also include others.
            - ``must_exclude`` – The ligand must not contain any of the specified elements.
            - ``must_only_contain_in_any_amount`` – The ligand may contain any count, including zero, of the specified elements, but must not contain any other elements.

            .. tip::

                On first glance, these instructions might seem too general to make a useful filter, but by combining the same filter multiple times with different instructions, users can achieve very specific filters.

        :param bool only_donors: If ``True``, evaluate the composition using only donor atoms; otherwise all atoms in the ligand are considered.
        :return: ``True`` if the ligand satisfies the specified composition criterion.
        :rtype: bool
        :raises ValueError: If an unrecognized instruction string is provided.

        **Examples :**

        To select tridentate ligands whose donors are exactly two N and one C:

        .. code-block:: python

            ligand.composition_filter(
                    elements='CN2',
                    instruction='must_contain_and_only_contain',
                    only_donors=True
                )

        To select ligands whose donors may contain only C and N atoms (or zero of either):

        .. code-block:: python

            ligand.composition_filter(
                    elements='CN',
                    instruction='must_only_contain_in_any_amount',
                    only_donors=True
                )

        To select ligands that contain zero or more C, N, and H atoms, but no other elements:

        .. code-block:: python

            ligand.composition_filter(
                elements='CNH',
                instruction='must_only_contain_in_any_amount',
                only_donors=False
            )

        To exclude ligands containing sulfur atoms:

        .. code-block:: python

            ligand.composition_filter(
                elements='S',
                instruction='must_exclude',
                only_donors=False
            )

        To select ligands that contain at least one O atom:

        .. code-block:: python

            ligand.composition_filter(
                elements='O',
                instruction='must_at_least_contain',
                only_donors=False
            )

        """
        if isinstance(elements, str):
            elements = stoichiometry2atomslist(elements)
        atoms_of_interest = [Element(el).symbol for el in elements]
        if only_donors:
            atoms = self.donor_elements
        else:
            atoms = self.atomic_props['atoms']

        if instruction == 'must_contain_and_only_contain':
            atoms_present = sorted(list(atoms)) == sorted(atoms_of_interest)
        elif instruction == 'must_at_least_contain':
            atoms_present = all(elem in list(atoms) for elem in atoms_of_interest)
        elif instruction == 'must_exclude':
            atoms_present = any(elem in list(atoms) for elem in atoms_of_interest) == False
        elif instruction == 'must_only_contain_in_any_amount':
            atoms_present = all(elem in atoms_of_interest for elem in list(atoms))
        else:
            raise ValueError(f'Invalid instruction: {instruction}. Must be one of ["must_contain_and_only_contain", "must_at_least_contain", "must_exclude", "must_only_contain_in_any_amount"].')

        return atoms_present

    def parents_filter(self, metal_centers: list[str]) -> bool:
        """
        Test whether the ligand has been observed in the CSD source complexes with any of the specified parent metal centers.

        While this filter does not directly check for chemistry, it is useful for maximizing compatibility of the ligand with new metal centers by selecting only those ligands that have previously coordinated to similar metals.

        :param list[str] metal_centers: List of metal center strings to check. Each string must be either a valid element symbol or an element symbol followed by an oxidation state.
        :return: True if any provided metal_center string matches recorded parent metals or oxidation states.
        :rtype: bool

        **Examples :**

        To filter for ligands that have been observed coordinating to either Fe or Cu(II):

        .. code-block:: python

            ligand.parents_filter(metal_centers=['Fe', 'Cu+2'])

        """
        for metal_center in metal_centers:
            check_metal_center_format(metal_center)

        parent_metals = list(self.metal_counts.keys())
        parent_mos = list(self.metal_os_counts.keys())
        all_parent_metal_centers = parent_metals + parent_mos
        has_metal_centers = any(metal in all_parent_metal_centers for metal in metal_centers)

        return has_metal_centers

    def smarts_filter(self, smarts: str, should_contain: bool, include_metal=None) -> bool:
        """
        Test whether the ligand matches (or not) a SMARTS substructure pattern.

        If :arg:`include_metal` is True, a pseudo metal center (Hg) is added to the ligand SMILES string before testing. The Hg atom is connected to all donor atoms via single bonds. This allows SMARTS patterns to target only donor atoms.

        :param str smarts: SMARTS pattern string to search for.
        :param bool should_contain: If True, ligand must contain the pattern; if False, ligand must not contain it.
        :param bool | None include_metal: If True include a Hg metal center in SMILES generation (default).
        :return: True if the ligand satisfies the SMARTS pattern condition, otherwise False.
        :rtype: bool

        .. tip::

            While SMARTS patterns are very useful to define chemical substructures, they can be difficult to create. However, AI assistants can help you, and you can verify the correctness of your SMARTS patterns using tools such as the `online SMARTS tester <https://smarts.plus/>`_.

        **Examples :**

        To include only amide ligands where the N donor is coordinated to the metal and two C or Si atoms:

        .. code-block:: python

            ligand.smarts_filter(smarts='[N&D3X3!a](-[Hg])(-[C,Si])(-[C,Si])', should_contain=True, include_metal=True)

        To exclude ligands that contain at least one azo (N=N) group:

        .. code-block:: python

            ligand.smarts_filter(smarts='[N]=[N]', should_contain=False)

        """
        include_metal = True if include_metal is None else include_metal
        smiles = self.smiles_with_metal if include_metal else self.smiles
        if smiles is None:  # If the ligand has no valid SMILES string, always fail the ligand.
            return False

        has_pattern = has_smarts_pattern(smarts=smarts, smiles=smiles)
        return has_pattern == should_contain

    def get_graph_with_metal(self, metal_symbol: Union[str, None]) -> tuple[nx.Graph, int]:
        """
        Return a copy of the ligand graph with an added metal node connected to donor atoms.

        The metal node receives a node attribute self.node_label equal to metal_symbol and
        metal‑donor edges are assigned bond_type = 1.

        :param metal_symbol: Element symbol of the metal to attach (e.g. 'Fe').
        :type metal_symbol: str | None
        :return: Tuple (graph_with_metal, metal_node_index).
        :rtype: tuple[nx.Graph, int]
        """
        graph_with_metal = nx.Graph(self.graph)   # unfreeze graph
        node_kwargs = {self.node_label: metal_symbol}
        bond_kwargs = {self.bond_label: 1}

        # Add metal node and bonds of coordinating atoms
        metal_idx = max(self.graph.nodes) + 1
        graph_with_metal.add_node(metal_idx, **node_kwargs)
        for atom_idx in self.donor_idc:
            # Indices of graph and atomic indices don't match
            graph_idx = atom_idx
            graph_with_metal.add_edge(metal_idx, graph_idx, **bond_kwargs)

        return graph_with_metal, metal_idx

    def get_graph_hash_with_metal(self, metal_symbol) -> str:
        """
        Return a graph hash for the ligand after attaching a specified metal node.

        :param metal_symbol: Metal element symbol used for hashing.
        :type metal_symbol: str
        :return: Graph hash string.
        :rtype: str
        """
        graph_with_metal, _ = self.get_graph_with_metal(metal_symbol=metal_symbol)
        return get_graph_hash(graph_with_metal, node_attr=self.node_label)

    def get_heavy_atoms_graph_hash_with_metal(self, metal_symbol) -> str:
        """
        Return the heavy-atom-only graph hash for the ligand when a metal is attached.

        :param metal_symbol: Metal element symbol used for hashing.
        :type metal_symbol: str
        :return: Hash for heavy-atom-only ligand+metal graph.
        :rtype: str
        """
        graph_with_metal, _ = self.get_graph_with_metal(metal_symbol=metal_symbol)
        heavy_graph_with_metal = get_heavy_atoms_graph(graph_with_metal)
        return get_graph_hash(heavy_graph_with_metal, node_attr=self.node_label)

    def get_xyz_string(self, comment: str='', with_metal: bool=True) -> str:
        """
        Return an XYZ-format string for the ligand, optionally including a pseudo-metal.

        :param comment: Optional comment line; if None a default informative comment is generated.
        :type comment: str
        :param with_metal: If True include a pseudo metal at recorded parent_metal_position.
        :type with_metal: bool
        :return: XYZ-format string (header + coordinates).
        :rtype: str
        """
        if comment is None: # default comment specifying important properties of the ligand
            donors = '-'.join(self.donor_elements)
            comment = f'Ligand ID: {self.unique_name}  ===  Stoichiometry: {self.stoichiometry}  ===  Charge: {self._get_charge_as_int()}  ===  Donors: {donors}'

        n_ligand_atoms = len(self.atomic_props['x'])
        if with_metal:
            str_ = f"{n_ligand_atoms+1}\n" # +1 for the metal
            str_ += comment + '\n'
            str_ += f"{pseudo_metal}  {self.parent_metal_position[0]}  {self.parent_metal_position[1]}  {self.parent_metal_position[2]} \n"     # metal atom
        else:
            str_ = f"{n_ligand_atoms}\n"
            str_ += comment + '\n'

        # Add ligand atoms
        for i, _ in enumerate(self.atomic_props['x']):
            str_ += f"{self.atomic_props['atoms'][i]}  {self.atomic_props['x'][i]}  {self.atomic_props['y'][i]}  {self.atomic_props['z'][i]} \n"

        # Remove the last newline character to match the .xyz format
        str_ = str_.rstrip('\n')

        return str_

    def get_all_effective_ligand_atoms_with_effective_donor_indices(self, dummy='Cu') -> tuple[ase.Atoms, list[int]]:
        """
        Build an effective ligand ASE Atoms object that replaces haptic groups by dummy atoms.

        The function returns (atoms_with_dummies, effective_donor_indices) where dummy atoms
        represent haptic groups to simplify geometric operations.

        :param dummy: Element symbol used for dummy atoms representing haptic groups.
        :type dummy: str
        :return: Tuple (ASE Atoms with dummy atoms appended, list of effective donor indices).
        :rtype: tuple[ase.Atoms, list[int]]
        """
        return get_all_effective_ligand_atoms_with_effective_donor_indices(
                                                                            ligand_atoms=self.atoms,
                                                                            hapdent_idc=self.hapdent_idc,
                                                                            dummy=dummy
                                                                            )

    def get_isomers_effective_ligand_atoms_with_effective_donor_indices(self, dummy='Cu') -> tuple[ase.Atoms, list[list[int]]]:
        """
        Produce effective-ligand ASE Atoms for each geometric isomer by replacing haptic groups with dummy atoms.

        Returns a tuple (atoms_with_dummies, isomer_effective_donor_idc) where the second element is a list
        of donor-index lists for each geometrical isomer, suitable for rotation/assembly routines.

        :param dummy: Element symbol to use for dummy atoms representing haptic groups.
        :type dummy: str
        :return: (ASE Atoms for effective ligand, list of donor-index lists for each isomer).
        :rtype: tuple[ase.Atoms, list[list[int]]]
        """
        all_atoms, isomers_eff_donor_idc = get_isomers_effective_ligand_atoms_with_effective_donor_indices(
                                                    ligand_atoms=self.atoms,
                                                    geometric_isomers_hapdent_idc=self.geometric_isomers_hapdent_idc,
                                                    dummy=dummy
                                                    )
        # Doublechecking that the effective atoms are the same as the original atoms if there are no haptic interactions
        if not self.is_haptic:
            assert all_atoms == self.atoms, f'Error in isomer generation for ligand {self.unique_name}.'

        return all_atoms, isomers_eff_donor_idc

    def get_ligand_archetype_and_isomers(self) -> tuple[str, list[ase.Atoms], tuple[Union[int, tuple[int]]], float, str, float]:
        """
        Determine ligand archetype and geometric isomers, handling haptic donors via dummies.

        The routine returns:
          (archetype, real_isomers, hapdent_isomer_idc, rssd, second_archetype, weight_for_change)

        :return: Tuple containing archetype string, list of ASE Atoms for best isomers, hapdent tuples for isomers,
                 RSSD float, second-best archetype string and weight needed to change archetype.
        :rtype: tuple[str, list[ase.Atoms], tuple[Union[int, tuple[int]]], float, str, float]
        """
        eff_ligand_atoms, eff_donor_idc = self.get_all_effective_ligand_atoms_with_effective_donor_indices('Cu')
        archetype, eff_isomers, eff_isomer_donor_idc, rssd, second_archetype, weight_necessary_for_change = assign_archetype(eff_ligand_atoms, eff_donor_idc)
        # Remove Cu from the isomers to get the real ligand archetype in case of haptic ligands
        real_isomers = []
        for isomer in eff_isomers:
            real_isomer = isomer[[atom.symbol != 'Cu' for atom in isomer]]
            real_isomers.append(real_isomer)

        # Convert the effective donor indices to the real donor indices with denticities and hapticities as sublists
        hapdent_isomer_idc = []
        for eff_isomer_idc in eff_isomer_donor_idc:
            hapdent_donor_idc = []
            for eff_idx in eff_isomer_idc:
                # Look up the list index of each effective index and mirror this to the hapdent index
                list_idx = eff_donor_idc.index(eff_idx)
                hapdent_idx = self.hapdent_idc[list_idx]
                hapdent_donor_idc.append(hapdent_idx)
            hapdent_donor_idc = tuple(hapdent_donor_idc)
            # Assert that the resulting indices are correct apart from the order
            assert set(hapdent_donor_idc) == set(self.hapdent_idc), f"Error in conversion of effective donor indices to hapdent donor indices: {hapdent_isomer_idc} vs. {self.hapdent_idc}"
            hapdent_isomer_idc.append(hapdent_donor_idc)

        return archetype, real_isomers, hapdent_isomer_idc, rssd, second_archetype, weight_necessary_for_change

    def get_csv_info(self, max_entries: int=5) -> dict:
        """
        Build a flat dictionary of ligand metadata suitable for CSV export.

        Long lists (e.g. parent complex IDs) are truncated to at most max_entries entries
        and represented as comma-separated strings with an ellipsis indicator.

        :param max_entries: Maximum number of list entries to present before truncation.
        :type max_entries: int
        :return: Dictionary containing ligand global properties augmented with CSV-friendly fields.
        :rtype: dict
        """
        info = self.global_props.copy()
        info['donors'] = '-'.join(sorted(self.donor_elements))  # Add donors as a string

        # Truncate lists with potentially many entries to a maximum of `max_entries`.
        truncated_data = {
                            'csd_complex_ids': self.other_ligand_instances['parent_complex_id'],
                            'csd_metal_count': [f'{el}({count})' for el, count in self.metal_counts.items()],
                            'csd_metal_os_count': [f'{el}({count})' for el, count in self.metal_os_counts.items()],
                        }
        for key, data in truncated_data.items():
            n_data = len(data)
            data = data[:max_entries]
            data = ', '.join(data)
            if n_data > max_entries:
                data += f', ... ({n_data - max_entries} more)'
            info[key] = data

        # Sort the dictionary to make the output csv more readable, e.g. insert 'donors' after 'stoichiometry'.
        ligand_csv_props_order = list(info.keys())
        ligand_csv_props_order.insert(ligand_csv_props_order.index('stoichiometry') + 1, 'donors')
        info = sort_dict_as(d=info, order=ligand_csv_props_order, strict=True)

        return info

    def to_dict(self, include_graph: bool=True, copy: bool=False, full_global_props: bool=False) -> dict:
        """
        Serialize the ligand to a dictionary matching the MetaLig format.

        :param include_graph: If True, include the connectivity graph in the output.
        :type include_graph: bool
        :param copy: If True, return a deep copy of the dictionary.
        :type copy: bool
        :param full_global_props: If True, compute and include all expected global_props fields before serializing.
        :type full_global_props: bool
        :return: Dictionary with keys identical to MetaLig ligand entries.
        :rtype: dict
        :raises AssertionError: If returned dictionary keys differ from the expected ligand_dict_props.
        """
        if full_global_props:
            for prop in ligand_global_props:
                self.global_props[prop] = getattr(self, prop)   # Will automatically be returned if `prop` already exists in `global_props`.

        self._sort_global_props_inplace()

        d = super().to_dict(include_graph=include_graph)
        d['donor_idc'] = self.donor_idc
        d['ligand_instances'] = self.other_ligand_instances
        d['hapdent_idc'] = self.hapdent_idc
        d['geometric_isomers_hapdent_idc'] = self.geometric_isomers_hapdent_idc

        if copy:
            d = deepcopy(d)

        # Make sure that the output dictionary looks identical to the MetaLig.
        unexpected_props = set(d.keys()).symmetric_difference(ligand_dict_props)
        assert not unexpected_props, f"Properties in ligand dictionary missing or unexpected: {unexpected_props}. Expected properties: {ligand_dict_props}."

        return d

    @classmethod
    def from_dict(cls, d: dict, validity_check: bool=True) -> 'Ligand':
        """
        Construct a Ligand instance from a MetaLig-style dictionary.

        The method handles backward-compatible refactoring of deprecated formats and validates
        presence of required fields before instantiation.

        :param dict d: Dictionary matching the MetaLig ligand entry format. Example keys:

            .. code-block:: python

                {
                    'atomic_props': {...},
                    'global_props': {...},
                    'graph': {...},
                    'donor_idc': [...],
                    'ligand_instances': {...},
                    'hapdent_idc': ...,
                    'geometric_isomers_hapdent_idc': ...
                }

        :param bool validity_check: Whether to run post-construction consistency checks.
        :return: Instantiated Ligand object.
        :rtype: Ligand
        :raises ValueError: If required top-level keys or global_props['charge'] are missing.
        """
        # Check whether the dictionary is in the old or new format.
        old_props = ['n_hydrogens', 'n_atoms', 'n_protons']
        if any([prop in d for prop in old_props]):
            warnings.warn('The provided ligand database has a format that was deprecated (DARTassembler version<1.1.0). An attempt is made to convert it to the new format, but this may not always work correctly. In case of problems, please re-generate the ligand database with DARTassembler version>=1.1.0.')
            d = refactor_metalig_entry_from_v1_0_0_to_v1_1_0(d)

        different_props = ligand_dict_props.symmetric_difference(set(d.keys()))
        if different_props:
            raise ValueError(f"Input dictionary is missing or has unexpected properties: {different_props}. Expected properties: {ligand_dict_props}.")

        if not 'charge' in d['global_props']:
            raise ValueError(f"Input dictionary is missing the 'charge' property in 'global_props'.")

        # Convert the graph from a dictionary to a NetworkX graph object
        graph = graph_from_graph_dict(d['graph'])

        return cls(
            atomic_props=d["atomic_props"],
            global_props=d["global_props"],
            unique_name=d['global_props']['unique_name'],
            donor_idc=d['donor_idc'],
            graph=graph,
            ligand_instances=d['ligand_instances'],
            hapdent_idc=d['hapdent_idc'],
            geometric_isomers_hapdent_idc=d['geometric_isomers_hapdent_idc'],
            charge=d['global_props']['charge'],
            validity_check=validity_check,
        )
    
    def _check_if_molecule_valid(self) -> None:
        """
        Run ligand-specific validation checks in addition to BaseMolecule consistency tests.

        Ensures the ligand graph is connected and that donor/hapdent bookkeeping is internally consistent.

        :return: None
        :rtype: None
        :raises ValueError: If the ligand graph is not connected or donor/hapdent counts are inconsistent.
        :raises AssertionError: If recalculated hapdent_idc differs from the stored one.
        """
        super()._check_if_molecule_valid()
        if not nx.is_connected(self.graph):
            raise ValueError(f'Graph of ligand {self.unique_name} is not fully connected.')
        if not self.n_donors == self.n_denticities + self.n_haptic_atoms:
            raise ValueError(
                f'Number of donors ({self.n_donors}) does not equal number of n_denticities ({self.n_denticities}) plus n_haptic_atoms ({self.n_haptic_atoms}) in ligand {self.unique_name}.')
        if hasattr(self, 'hapdent_idc'):
            new_hapdent = get_denticities_and_hapticities_idc(graph=self.graph, donor_idc=self.donor_idc)
            assert new_hapdent == self.hapdent_idc, f'Hapdent_idc are inconsistent: old: {self.hapdent_idc}, new: {self.hapdent_idc}'

        return None
    
    def _sort_global_props_inplace(self) -> None:
        """
        Reorder self.global_props to the canonical ligand_global_props ordering for consistent serialization.

        :return: None
        :rtype: None
        """
        self.global_props = sort_dict_as(d=self.global_props, order=ligand_global_props, strict=False)
    
    def _get_rdkit_mol_from_smiles(self, sanitize: bool=False, with_metal: str=None) -> Union['rdkit.Chem.Mol',None]:
        """
        Build an RDKit molecule from the ligand SMILES, optionally connecting a metal.

        :param sanitize: If True, apply RDKit sanitization during parsing.
        :type sanitize: bool
        :param with_metal: If provided, attach this metal symbol to the ligand for SMILES generation.
        :type with_metal: str | None
        :return: RDKit Chem.Mol instance or None when SMILES is unavailable.
        :rtype: rdkit.Chem.Mol | None
        :raises ImportError: If RDKit is not installed.
        """
        try:
            from rdkit import Chem
        except (ImportError, ModuleNotFoundError):
            raise ImportError("This function requires the RDKit package to be installed.")
        smiles = self.get_smiles(with_metal=with_metal)
        if smiles is None:
            return None

        return Chem.MolFromSmiles(smiles, sanitize=sanitize)

    def _get_atomic_distance_to_original_metal(self, mode: str = 'min') -> Union[float, list, tuple]:
        """
        Compute distances from every atom in the ligand to the recorded parent metal position.

        :param mode: One of 'min', 'max', 'coordinating', 'all' selecting returned summary.
        :type mode: str
        :return: Depending on mode:
                 - 'min' -> float minimum distance,
                 - 'max' -> float maximum distance,
                 - 'coordinating' -> list[float] distances for coordinating atoms,
                 - 'all' -> tuple(min, max, list[float]).
        :rtype: float | list[float] | tuple
        """
        mode = mode.lower()

        distances = np.linalg.norm(self.atoms.positions - np.array(self.parent_metal_position), axis=1)
        if mode == 'min':
            return distances.min()
        elif mode == 'max':
            return distances.max()
        elif mode == 'coordinating':
            return distances[self.donor_idc].tolist()
        elif mode == 'all':
            return distances.min(), distances.max(), distances[self.donor_idc].tolist()

    def _get_charge_as_int(self) -> Union[int, 'nan']:
        """
        Convert stored charge to integer if numeric, otherwise return np.nan.

        :return: Integer charge or np.nan if unknown.
        :rtype: int | nan
        """
        if np.isnan(self.charge):
            return self.charge
        else:
            return int(self.charge)
        
    def _check_if_2D_symmetrical(self) -> bool:
        """
        Assess 2D donor-to-donor symmetry by attaching a pseudo-metal and comparing donor-cut graphs.

        The method generates donor-specific graphs where one donor is disconnected from a central pseudo-metal
        and compares WL-hashes to detect symmetry between donors.

        :return: True when at least two donor-specific graphs are identical (2D donor symmetry), else False.
        :rtype: bool
        """
        # todo: think about how to do this with haptic ligands
        # Explanation of algorithm: This function checks if the ligand graph is symmetrical between any of the donors. Essentially, it attaches a pseudo Hg atom to each donor and checks if the resulting ligand graph between any donors is identical. If it is, the ligand is symmetrical under these donors.

        # Get graph where all donors are connected to a pseudo Hg atom.
        graph, metal_idx = self.get_graph_with_metal(metal_symbol='Hg')

        # Make new graphs which each have one pseudo Hg atom bonding to all but one donor, i.e. remove one donor bond.
        donor_graphs = []
        donor_indices = list(graph.neighbors(metal_idx))
        for donor_idx in donor_indices:
            donor_graph = graph.copy()
            donor_graph.remove_edge(metal_idx, donor_idx)
            donor_graphs.append(donor_graph)

        # Calculate graph hashes of all graphs. If any of the graph hashes are identical, the ligand is symmetrical between at least two donors.
        graph_hashes = [get_graph_hash(donor_graph) for donor_graph in donor_graphs]
        symmetrical = len(set(graph_hashes)) < len(graph_hashes)

        return symmetrical
    
    def _get_n_beta_hydrogens(self) -> int:
        """
        Count hydrogen atoms that are exactly two bonds away from any coordinating donor (beta hydrogens).

        Uses the square of the adjacency matrix to find nodes at graph distance two while excluding alpha hydrogens.

        :return: Integer count of beta hydrogens.
        :rtype: int
        """
        # The second power of the adjacency matrix, i.e. A^2[i,j] represents the number of paths of length two from i to j. Hence, as we are only interested in Hydrogen which has a distance of two to our coordinating atoms.
        A = get_adjacency_matrix(self.graph)
        B = np.matmul(A, A)

        beta_h_indices = set()  # Using the set we avoid double counting of beta H atoms
        for donor_idx in self.donor_idc:
            for idx, element in self.graph.nodes(data=self.node_label):
                if element == "H":
                    # search for beta H while excluding alpha H
                    if B[donor_idx, idx] > 0 and A[donor_idx, idx] == 0:
                        beta_h_indices.add(idx)

        n_beta_hydrogens = len(beta_h_indices)

        return n_beta_hydrogens
    
    def _get_donors_planarity(self, with_metal: bool=True) -> float:
        """
        Compute planarity for donor atoms, optionally including the parent metal.

        For sets of fewer than three atoms the ligand is considered perfectly planar.

        :param with_metal: Include the parent metal in the planarity calculation when True.
        :type with_metal: bool
        :return: Planarity score between 0.0 and 1.0.
        :rtype: float
        """
        # Get the coordinates of the donors and the metal
        coordinates = self.get_coordinates()
        coordinates = [coordinates[donor_idx] for donor_idx in self.donor_idc]
        if with_metal:
            coordinates.append(self.parent_metal_position)

        # If there are less than 3 atoms, the molecule is planar by definition
        if len(coordinates) < 3:
            return 1.0

        planarity = get_planarity(coordinates)

        return planarity
    
    def _get_denticities_and_hapticities_idc(self) -> tuple[Union[int, tuple[int]]]:
        """
        Compute the denticity/hapticity index tuple for the ligand from its graph and donor list.

        Returns a tuple where integer entries denote monodentate donors and sub-tuples denote haptic groups.

        :return: Tuple combining ints and sub-tuples that represent denticities and hapticities.
        :rtype: tuple[Union[int, tuple[int]]]
        :raises ValueError: If an entry in computed hapdent_idc has an unexpected type.
        """
        hapdent_idc = get_denticities_and_hapticities_idc(graph=self.graph, donor_idc=self.donor_idc)

        # Convert back to atomic indices. Keep in mind that some entries are single indices and some are tuples of indices.
        atomic_hapdent_idc = []
        for idc in hapdent_idc:
            if isinstance(idc, int):    # n_donors, therefore integer
                idc = idc
                assert idc in self.donor_idc, f'Index {idc} is not in the ligand to metal indices: {self.donor_idc}'
            elif isinstance(idc, tuple):
                assert all([i in self.donor_idc for i in idc]), f'Indices {idc} is not in the ligand to metal indices: {self.donor_idc}'
            else:
                raise ValueError(f'Invalid type of hapdent_idc: {type(idc)}')
            atomic_hapdent_idc.append(idc)

        return tuple(atomic_hapdent_idc)


from functools import cached_property
import warnings
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
    check_metal_center_format, get_atomic_props_from_ase_atoms, get_ase_atoms_from_atomic_props
from DARTassembler.src.constants.chem import Element
from DARTassembler.src.metalig.utils_graph import graph_from_graph_dict, graph_to_dict_with_node_labels, view_graph, \
    get_sorted_atoms_and_indices_from_graph, get_reindexed_graph, get_graph_fragments, count_atoms_with_n_bonds, get_graph_hash, get_heavy_atoms_graph, \
    get_adjacency_matrix, assert_graph_and_coordinates_are_consistent, \
    remove_node_features_from_graph, make_multigraph_to_graph, graph_to_smiles
from DARTassembler.src.metalig.utils import get_stable_sorted_value_counts, check_equal, sort_dict_as
from DARTassembler.src.metalig.utils_molecule import get_standardized_stoichiometry_from_atoms_list, \
    unknown_rdkit_bond_orders, get_max_deviation_from_coplanarity, if_same_stoichiometries
from DARTassembler.src.metalig.geometry import assign_geometry

pseudo_metal = 'Cu'     # Pseudo metal for display in ligand xyz files and for use in the SMARTS filter.
# The set of properties that are stored in the ligand dictionary of the MetaLig database.
ligand_dict_props = {'atomic_props', 'global_props', 'graph', 'donor_idc', 'ligand_instances', 'hapdent_idc',
                     'geometric_isomers_hapdent_idc'}
# The set of properties in the ligand_instances dictionary of each ligand in the MetaLig database.
ligand_instance_props = {'ligand_name', 'parent_complex_id', 'parent_complex_charge', 'parent_metal', 'parent_metal_os'}
# The properties expected in the ligand `global_props`. Must be a list instead of a set, because the order matters for the csv output.
ligand_global_props = [
    # General properties
    'unique_name',
    'stoichiometry',
    'geometry',
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
    'geometry_rssd',
    'geometry_confidence',
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

# Helper decorator for caching properties in the Ligand() instance defined later in this script.
class cached_global_props(cached_property):
    """
    A decorator for lazy evaluation and subsequent caching of properties. In contrast to the standard cached_property, this function also looks up and stores each value in the `global_props` dictionary of the ligand instance, if the property name is in the list of expected global properties.
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
        This is a base class for a molecule, which can be either a ligand or a complex.
        :param graph: The molecular graph represented as a networkx Graph object. The nodes are indexed from 0 to n_atoms-1.
        :param atomic_props: ase.Atoms object or atomic properties dictionary in the following format:
        atomic_props = {
            'atoms': ['C', 'O', 'N', ...],
            'x': [0.0, 1.0, 2.0, ...],
            'y': [0.0, 1.0, 2.0, ...],
            'z': [0.0, 1.0, 2.0, ...],
            ... -> other properties can be added as needed, e.g. 'partial_charges': [0.0, -0.5, 0.5, ...],
        }
        :param global_props: A flat dictionary containing global properties of the molecule, such as charge, stoichiometry, etc. If None, an empty dictionary is created.
        :param validity_check: If True, makes a few consistency checks. If False, skips the validity check.
        :param node_label: The label for the nodes in the graph, which specifies the atom type.
        :param bond_label: The label for the edges in the graph, which specifies the bond type.
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

    @cached_global_props
    def n_atoms(self) -> int:
        """The number of atoms in the molecule."""
        return len(self.atomic_props['atoms'])

    @cached_global_props
    def n_hydrogens(self) -> int:
        """The number of hydrogen atoms in the molecule."""
        return sum(1 for atom in self.atomic_props['atoms'] if atom == 'H')

    @cached_global_props
    def n_protons(self) -> int:
        """The number of protons in the molecule."""
        return sum([Element(el).atomic_number for el in self.atomic_props['atoms']])

    @cached_global_props
    def n_bonds(self) -> int:
        """The number of bonds in the molecule."""
        return len(self.graph.edges)

    @cached_global_props
    def has_bond_order_attribute(self) -> bool:
        """If the molecular graph has a bond order attribute."""
        return all(self.bond_label in self.graph.edges[edge] for edge in self.graph.edges)

    @cached_global_props
    def has_unknown_bond_orders(self) -> bool:
        """If the molecular graph has any unknown bond orders."""
        if not self.has_bond_order_attribute:
            return True
        return self.count_bond_types(unknown_rdkit_bond_orders) > 0

    @cached_global_props
    def has_all_bond_orders_valid(self) -> bool:
        """If all bond orders in the molecular graph are valid and known."""
        return self.has_bond_order_attribute and not self.has_unknown_bond_orders

    @cached_global_props
    def graph_hash(self) -> str:
        """The Weisfeiler-Lehman hash of the molecular graph. This is a unique fingerprint for the graph structure, taking into account the chemical element of each atom and the connections between atoms, but not bond orders."""
        return get_graph_hash(self.graph)

    @cached_global_props
    def heavy_atoms_graph_hash(self) -> str:
        """The Weisfeiler-Lehman hash of the molecular graph with only heavy atoms (no hydrogens)."""
        heavy_graph = get_heavy_atoms_graph(self.graph, element_label=self.node_label)
        return get_graph_hash(heavy_graph)

    @cached_global_props
    def bond_order_graph_hash(self) -> Union[str, None]:
        """The Weisfeiler-Lehman hash of the molecular graph, also taking into account the bond orders. Returns None if not all bond orders are valid."""
        if not self.has_all_bond_orders_valid:
            return None
        return get_graph_hash(self.graph, node_attr=self.node_label, edge_attr=self.bond_label)

    @cached_global_props
    def stoichiometry(self) -> str:
        """The stoichiometry of the molecule."""
        return get_standardized_stoichiometry_from_atoms_list(self.atomic_props['atoms'])

    def get_reindexed_graph(self) -> nx.Graph:
        """
        Returns the reindexed graph in which the nodes are indexed from 0 to n_atoms-1. With ligands in the MetaLig this should be given by default from DARTassembler version 1.1.0 onwards.
        :return: reindexed graph
        """
        return get_reindexed_graph(self.graph)

    def get_ase_atoms(self, remove_elements: list=None, add_atoms: list=None) -> Atoms:
        """
        Get an ase.Atoms object of the ligand. Optionally remove certain chemical elements or add atoms to the molecule.
        :param remove_elements: list of chemical elements to optionally remove from the molecule, e.g. ['H'] will remove all hydrogen atoms from the output molecule.
        :param add_atoms: list of tuples of the form [(element, (x,y,z))] to add atoms to the output molecule.
        :return: ASE Atoms object of the molecule with the requested modifications.
        """
        return get_ase_atoms_from_atomic_props(atomic_props=self.atomic_props, remove_elements=remove_elements,add_atoms=add_atoms)

    def get_smiles(self) -> Union[str,None]:
        """
        Returns the SMILES string of the molecule.
        :return: SMILES string or None if not all bond orders are valid.
        """
        if not self.has_all_bond_orders_valid: # if the molecule has unknown bond orders, we cannot calculate the SMILES
            return None

        smiles = graph_to_smiles(self.graph, element_label=self.node_label, bond_label=self.bond_label)

        return smiles

    def count_C_H_bonds(self) -> int:
        """
        Returns the number of C-H bonds in the molecule.
        :return: number of C-H bonds in the molecule.
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
        Counts the number of specified bond types in the molecular graph.
        :param bond_types: list of integers of rdkit bond types
        :return: True if the molecular graph has any of the specified bond types, False otherwise
        :raises ValueError: if the molecular graph does not have a bond order attribute
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
        Returns a list of the fragment indices (unconnected components) and their elements of the molecular graph.
        :return: tuple of two lists, the first list contains the indices of the fragments, the second list contains the elements of the fragments.
        """
        return get_graph_fragments(graph=self.graph, node_label=self.node_label)

    def view_3D(self) -> None:
        """
        Opens a 3D visualization of the molecule using ASE's view function.
        :return: None
        """
        view(self.atoms)

    def view_graph(self, node_size=150) -> None:
        """
        Opens a 3D visualization of the molecular graph using ASE's view function.
        :param node_size: size of the nodes in the graph
        :return: None
        """
        view_graph(self.graph, node_label=self.node_label, node_size=node_size)

    def get_coordinates(self) -> list:
        """
        Returns the coordinates of the molecule.
        :return: Coordinates in the format [[x1, y1, z1], [x2, y2, z2], ...]
        """
        return [[self.atomic_props['x'][i], self.atomic_props['y'][i], self.atomic_props['z'][i]] for i in range(len(self.atomic_props['x']))]

    def get_interatomic_distances(self, skip_elements=None) -> tuple[float, float, np.ndarray]:
        """
        Returns the minimum, maximum and all distances between two atoms in the molecule.
        :param skip_elements: Do not include bonds with these elements in the calculation
        :return: tuple of (minimum, maximum, all) distance of atoms in Angstrom
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
        Returns the distances between all atoms in the molecule as a flat list.
        :return: list of distances between all atoms in the molecule.
        """
        distances = self.get_interatomic_distances()[2]
        unique_distances = []
        for i in range(len(distances)):
            for j in range(i+1, len(distances)):
                unique_distances.append(distances[i,j])

        return unique_distances

    def get_xyz_string(self, comment: str= '') -> str:
        """
        Returns the .xyz file format string of the molecule.
        :param comment: comment to be added to the .xyz file
        :return: .xyz file format string of the molecule
        """
        xyz = f"{len(self.atomic_props['x'])}\n"
        xyz += comment + '\n'
        for i, _ in enumerate(self.atomic_props['x']):
            xyz += f"{self.atomic_props['atoms'][i]}  {self.atomic_props['x'][i]}  {self.atomic_props['y'][i]}  {self.atomic_props['z'][i]} \n"

        # Remove trailing newline character
        xyz = xyz.rstrip('\n')

        return xyz

    def to_dict(self, include_graph: bool=True) -> dict:
        """
        Returns the molecule as a dictionary.
        :param include_graph: if True, the graph is included in the dictionary, otherwise only atomic and global properties are returned.
        :return: dictionary representation of the molecule with atomic properties, global properties and optionally the graph.
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
        Appends the molecule data to a file using the specified writer.
        :param key: key under which the molecule data will be stored in the file.
        :param writer: writer object that has a write method to write the data to the file.
        :return: None
        """
        data = {'key': key, 'value': self.to_dict()}
        writer.write(data)

        return

    def _count_atoms_with_n_bonds(self, element: Union[str, None], n_bonds: int) -> int:
        """
        Count the number of n_ligand_instances of element `element` with exactly `n_bonds` bonds.
        :param element (str, None): specification of the element, e.g. 'C'. If None, all elements are counted.
        :param n_bonds (int): count an atom if it has exactly this number of bonds
        :return (int): integer count of the n_ligand_instances
        """
        return count_atoms_with_n_bonds(graph=self.graph, element=element, n_bonds=n_bonds, graph_element_label=self.node_label)

    def _contains_only(self, atoms: Union[str, list], except_elements=None) -> bool:
        """
        Returns True if the molecule contains only elements out of the list `atoms`, otherwise False.
        :param atoms: list of elements to check for
        :param except_elements: ignore these elements in the molecule when testing
        :return:
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
        Checks whether the molecule is valid. If not, raises an error.
        :return: None
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
        Checks if the stoichiometry of the molecule is the same as the given stoichiometry.
        :param other: stoichiometry to compare to.
        :return: True if the stoichiometry is the same, False otherwise.
        """
        return if_same_stoichiometries(self.stoichiometry, other)

    def _get_rdkit_mol_from_smiles(self, sanitize: bool=False):
        """
        Returns the RDKit molecule constructed from the SMILES string of the molecule. Requires an installed RDKit package.
        :return: RDKit molecule from the given SMILES string. Returns None if the molecule has unknown bond orders and therefore no SMILES can be calculated.
        """
        from rdkit import Chem
        smiles = self.get_smiles()
        if smiles is None:
            return None

        return Chem.MolFromSmiles(smiles, sanitize=sanitize)

    def _get_planarity(self) -> float:
        """
        Calculates the planarity of all atoms in the molecule.
        :return: Planarity of the molecule as a float between 0 and 1. 0 means not planar at all (a 3D sphere), 1 means perfectly planar.
        """
        coordinates = self.get_coordinates()
        deviation = get_max_deviation_from_coplanarity(points=coordinates)  # deviation is a float that is 0 if the molecule is perfectly planar and > 0 if it is not. The higher the value, the less planar the molecule is.
        planarity = 1/ (1+ deviation)   # planarity is a float between 0 and 1. 0 means not planar at all (a sphere), 1 means perfectly planar.
        planarity = round(planarity, 10)    # round to 10 decimal places to avoid floating point errors which happen with np.linalg.svd() in different versions of numpy

        return planarity

    def _get_xyz_array(self) -> np.ndarray:
        """
        Returns the coordinates of the molecule as a numpy array in the format [[x1, y1, z1], [x2, y2, z2], ...].
        :return: numpy array of the coordinates of the molecule
        """
        return np.array([self.atomic_props['x'], self.atomic_props['y'], self.atomic_props['z']]).T

    def _ligand_naming(self, n_donors: int, ligand_list: list) -> (str, str):
        """
        Generates a unique name for the ligand based on its n_donors and the existing ligands in the ligand_list.
        :param n_donors: Number of donor atoms.
        :param ligand_list: list of existing ligands in the database.
        :return: ligand_name: unique name for the ligand, csd: CSD code of the ligand.
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
        Checks whether `inherit_global_properties` has correct input format.
        :param inherit_global_properties: input
        :return: correctly specified inherit_global_properties
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
        Removes all node features from all molecular graphs in the db except the node features specified in keep.
        :param keep: list of node features which will not be removed
        :return: None
        """
        if keep is None:
            keep = [self.node_label]

        remove_node_features_from_graph(graph=self.graph, keep=keep, inplace=True)

        return

    def _normalize_multigraph_into_graph_inplace(self) -> None:
        """
        If self.graph is a MultiGraph, it is normalized into a Graph.
        :return: None
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
        A Ligand object, usually from the MetaLig database.
        :param atomic_props: ase.Atoms object or atomic properties dictionary in the following format:
        atomic_props = {
            'atoms': ['C', 'O', 'N', ...],
            'x': [0.0, 1.0, 2.0, ...],
            'y': [0.0, 1.0, 2.0, ...],
            'z': [0.0, 1.0, 2.0, ...],
            ... -> other properties can be added as needed, e.g. 'partial_charges': [0.0, -0.5, 0.5, ...],
        }
        :param donor_idc: List of indices of donor atoms in the ligand, i.e. the atoms that are connected to the metal center.
        :param graph: The molecular graph represented as a networkx Graph object. The nodes are indexed from 0 to n_atoms-1.
        :param unique_name: Unique name for the ligand, which is used to identify it in the database, e.g. 'unq_CSD-OZIYON-02-a'.
        :param charge: Formal charge of the ligand. Can be set to np.nan for many applications if the charge is not known or not relevant.
        :param global_props: A flat dictionary containing global properties of the ligand, such as charge, stoichiometry, etc. If None, an empty dictionary is created.
        :param ligand_instances: A dictionary containing information about the ligand instances in the parent complex in the following format:
        ligand_instances = {
            'ligand_name': ['ligand1', 'ligand2', ...],
            'parent_complex_id': ['complex1', 'complex2', ...],
            'parent_complex_charge': [0, 1, ...],
            'parent_metal': ['Cu', 'Ni', ...],
            'parent_metal_os': [2, 3, ...],
        }
        :param hapdent_idc: A tuple of denticity and hapticities indices, where denticity is an integer and hapticities are tuples of indices of atoms that are part of the same haptic group:
        hapdent_idc = (0, (2,3), 6) -> a tridentate ligand with two dentic donor atoms (0 and 1) and a haptic group with atoms 2 and 3.
        :param geometric_isomers_hapdent_idc: A list of hapdent_idc, where each entry has a different order of elements in the hapdent_idc, representing geometric isomers of the ligand:
        geometric_isomers_hapdent_idc = [(0, (2,3), 6), (6, (2,3), 0)] -> two geometric isomers, here of a mer tridentate ligand. The outer donor atoms (0 and 6) are swapped to generate the two isomers.
        :param validity_check: If True, makes a few consistency checks. If False, skips the validity check.

        ==== Overview of properties ====
         +-------------------------------+----------------------+------------------------------------------------------+
        | Name                          | Type                 | Description                                          |
        +===============================+======================+======================================================+
        | atomic_props                  | Union[dict, ase.Atoms] | ASE Atoms object or dictionary of atomic properties |
        |                               |                      | with keys 'atoms', 'x', 'y', 'z', etc.               |
        +-------------------------------+----------------------+------------------------------------------------------+
        | graph                         | nx.Graph             | Molecular graph.                                   |
        +-------------------------------+----------------------+------------------------------------------------------+
        | atoms                         | ase.Atoms            | ASE Atoms object of the ligand derived from atomic_props. |
        +-------------------------------+----------------------+------------------------------------------------------+
        | global_props                  | dict                 | Flat dictionary of global ligand properties.         |
        +-------------------------------+----------------------+------------------------------------------------------+
        | node_label                    | str                  | Label for nodes in the graph, indicating atom type. |
        +-------------------------------+----------------------+------------------------------------------------------+
        | bond_label                    | str                  | Label for edges in the graph, indicating bond type. |
        +-------------------------------+----------------------+------------------------------------------------------+
        | unique_name                   | str                  | Unique identifier for the ligand.                    |
        +-------------------------------+----------------------+------------------------------------------------------+
        | stoichiometry                 | str                  | Stoichiometry of the ligand (e.g., 'C10H12N2O4').    |
        +-------------------------------+----------------------+------------------------------------------------------+
        | geometry                      | str                  | Assigned ligand geometry (e.g., '2_cis').            |
        +-------------------------------+----------------------+------------------------------------------------------+
        | charge                        | Union[int, np.nan] | Formal charge of the ligand.                     |
        +-------------------------------+----------------------+------------------------------------------------------+
        | smiles                        | Union[str, None]     | SMILES string if bond orders valid; else None.      |
        +-------------------------------+----------------------+------------------------------------------------------+
        | smiles_with_metal             | Union[str, None]     | SMILES string including pseudo-metal; else None.     |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_donors                      | int                  | Number of donor atoms in the ligand.                 |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_eff_denticities             | int                  | Effective denticity counting each haptic group as one.|
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_denticities                 | int                  | Denticity count of non-haptic donor atoms.           |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_haptic_atoms                | int                  | Total number of haptically coordinating atoms.       |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_haptic_groups               | int                  | Number of haptic groups (hapticity).                 |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_atoms                       | int                  | Total number of atoms in the ligand.                 |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_elements                    | int                  | Number of unique elements in the ligand.             |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_bonds                       | int                  | Number of bonds in the ligand’s graph.               |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_electrons                   | int                  | Total number of electrons in the ligand.           |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_protons                     | int                  | Total number of protons in the ligand.               |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_beta_hydrogens              | int                  | Number of beta hydrogen atoms.                     |
        +-------------------------------+----------------------+------------------------------------------------------+
        | molecular_weight              | float                | Molecular weight of the ligand.                      |
        +-------------------------------+----------------------+------------------------------------------------------+
        | planarity                     | float                | Overall planarity of the ligand (0 to 1.0).          |
        +-------------------------------+----------------------+------------------------------------------------------+
        | donor_planarity               | float                | Planarity of donor atoms excluding the metal.        |
        +-------------------------------+----------------------+------------------------------------------------------+
        | donor_metal_planarity         | float                | Planarity including donor atoms and metal center.    |
        +-------------------------------+----------------------+------------------------------------------------------+
        | min_interatomic_distance      | float                | Minimum distance among ligand atoms.                |
        +-------------------------------+----------------------+------------------------------------------------------+
        | max_ligand_extension          | float                | Maximum distance between any two ligand atoms.       |
        +-------------------------------+----------------------+------------------------------------------------------+
        | geometry_rssd                 | float                | RSSD value for geometry assignment.                  |
        +-------------------------------+----------------------+------------------------------------------------------+
        | geometry_confidence           | float                | Confidence score for geometry assignment.            |
        +-------------------------------+----------------------+------------------------------------------------------+
        | is_2D_symmetrical             | bool                 | True if the ligand graph is symmetrical between donors. |
        +-------------------------------+----------------------+------------------------------------------------------+
        | has_all_bond_orders_valid     | bool                 | True if all bond orders in the graph are known.      |
        +-------------------------------+----------------------+------------------------------------------------------+
        | graph_hash                    | str                  | Weisfeiler-Lehman hash of the ligand graph (no metal). |
        +-------------------------------+----------------------+------------------------------------------------------+
        | graph_hash_with_metal         | str                  | Hash of the ligand graph including pseudo-metal.     |
        +-------------------------------+----------------------+------------------------------------------------------+
        | heavy_atoms_graph_hash        | str                  | Hash of the heavy-atom-only graph (no hydrogens).    |
        +-------------------------------+----------------------+------------------------------------------------------+
        | heavy_atoms_graph_hash_with_metal | str              | Hash of heavy-atom-only graph with pseudo-metal.     |
        +-------------------------------+----------------------+------------------------------------------------------+
        | bond_order_graph_hash         | Union[str, None]     | Hash of graph including bond orders; else None.      |
        +-------------------------------+----------------------+------------------------------------------------------+
        | n_ligand_instances            | int                  | Number of instances in parent complexes.              |
        +-------------------------------+----------------------+------------------------------------------------------+
        | donor_idc                     | List[int]            | Indices of donor atoms coordinating the metal center. |
        +-------------------------------+----------------------+------------------------------------------------------+
        | other_ligand_instances        | dict                 | Ligand-instance data in parent complexes.             |
        +-------------------------------+----------------------+------------------------------------------------------+
        | hapdent_idc                   | Tuple[Union[int, Tuple[int]]] | Tuple of denticity and haptic indices.      |
        +-------------------------------+----------------------+------------------------------------------------------+
        | geometric_isomers_hapdent_idc | List[Tuple[Union[int, Tuple[int]]]] | Geometric isomer hapdent tuples.       |
        +-------------------------------+----------------------+------------------------------------------------------+
        | donor_elements                | List[str]            | Chemical symbols of donor atoms.                      |
        +-------------------------------+----------------------+------------------------------------------------------+
        | donor_positions               | np.ndarray           | 3D coordinates of donor atoms.                        |
        +-------------------------------+----------------------+------------------------------------------------------+
        | parent_complex_id             | str                  | Parent complex ID for the ligand.                     |
        +-------------------------------+----------------------+------------------------------------------------------+
        | parent_metal                  | str                  | Chemical symbol of the parent metal.                  |
        +-------------------------------+----------------------+------------------------------------------------------+
        | parent_metal_position         | np.ndarray           | 3D coordinates of the parent metal center.            |
        +-------------------------------+----------------------+------------------------------------------------------+
        | metal_counts                  | dict                 | Counts of how often ligand binds to each metal.       |
        +-------------------------------+----------------------+------------------------------------------------------+
        | metal_os_counts               | dict                 | Counts of how often ligand binds with each oxidation state. |
        +-------------------------------+----------------------+------------------------------------------------------+
        | is_haptic                     | bool                 | True if ligand has any haptic atoms.                  |
        +-------------------------------+----------------------+------------------------------------------------------+
        +-------------------------------+----------------------+------------------------------------------------------+

        ==== Overview of methods ====
        - get_ase_atoms: Returns the ligand as an ase.Atoms object.
        - get_smiles: Returns the SMILES string of the ligand, if all bond orders are valid.
        - get_graph_fragments: Returns a list of the fragment indices (unconnected components) and their elements of the molecular graph.
        - get_coordinates: Returns the coordinates of the ligand as a list of lists in the format [[x1, y1, z1], [x2, y2, z2], ...].
        - get_interatomic_distances: Returns the minimum, maximum and all distances between two atoms in the ligand.
        - get_all_interatomic_distances_flat: Returns the distances between all atoms in the ligand as a flat list.
        - get_xyz_string: Returns the .xyz file format string of the ligand.
        - to_dict: Returns the ligand as a dictionary with atomic properties, global properties and optionally the graph.
        - view_3D: Opens a 3D visualization of the ligand using ASE's view function.
        - view_graph: Opens a 3D visualization of the molecular graph using ASE's view function.
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

    @cached_global_props
    def has_confident_charge(self) -> bool:
        if 'has_confident_charge' in self.global_props:
            return self.global_props['has_confident_charge']
        else:
            return True

    @cached_global_props
    def n_ligand_instances(self) -> int:
        """The number of instances of the ligand in the parent complexes."""
        return len(self.other_ligand_instances['ligand_name'])

    @cached_global_props
    def donor_metal_planarity(self) -> float:
        """The planarity of the donor atoms including the metal center."""
        return self._get_donors_planarity(with_metal=True)

    @cached_global_props
    def donor_planarity(self) -> float:
        """The planarity of the donor atoms excluding the metal center."""
        return self._get_donors_planarity(with_metal=False)

    @cached_global_props
    def planarity(self) -> float:
        """The overall planarity of the ligand."""
        return self._get_planarity()

    @cached_global_props
    def n_donors(self):
        """The number of donor atoms in the ligand."""
        return len(self.donor_idc)

    @cached_global_props
    def donor_elements(self) -> list[str]:
        """The chemical synbols of the donor atoms in the ligand."""
        return [self.atomic_props['atoms'][i] for i in self.donor_idc]

    @cached_global_props
    def donor_positions(self) -> np.ndarray:
        """The 3D coordinates of the donor atoms in the ligand."""
        return np.array([[self.atomic_props["x"][i], self.atomic_props["y"][i], self.atomic_props["z"][i]] for i in self.donor_idc])

    @cached_global_props
    def parent_complex_id(self) -> str:
        """The parent complex ID of the ligand."""
        return self.other_ligand_instances['parent_complex_id'][0]

    @cached_global_props
    def parent_metal_position(self) -> np.ndarray:
        """The 3D coordinates of the parent metal center of the ligand."""
        return self.other_ligand_instances['parent_metal_position'][0]

    @cached_global_props
    def parent_metal(self) -> str:
        """The chemical symbol of the parent metal of the ligand."""
        return self.other_ligand_instances['parent_metal'][0]

    @cached_global_props
    def hapdent_idc(self) -> tuple:
        """The hapdent_idc of the ligand, which is a tuple containing the dentic indices and in sub-tuples the haptic indices."""
        return self._get_denticities_and_hapticities_idc()

    @cached_global_props
    def n_eff_denticities(self) -> int:
        """The effective denticity of the ligand, in which each haptic group counts as one denticity."""
        return len(self.hapdent_idc)

    @cached_global_props
    def n_denticities(self) -> int:
        """The denticity of the ligand, without considering haptic groups."""
        return sum([1 for el in self.hapdent_idc if isinstance(el, int)])

    @cached_global_props
    def n_haptic_atoms(self) -> int:
        """The total number of haptically coordinating atoms in the ligand. This is different from the hapticity, which is the number of haptic groups."""
        return sum([len(sublist) for sublist in self.hapdent_idc if isinstance(sublist, tuple)])

    @cached_global_props
    def n_haptic_groups(self) -> int:
        """The number of haptic groups in the ligand, usually called the hapticity."""
        return len([sublist for sublist in self.hapdent_idc if isinstance(sublist, tuple)])

    @cached_global_props
    def _geometry_and_geometrical_isomers(self):
        """Cache all the geometry and isomer information when required."""
        geometry, _, isomer_hapdent_idc, rssd, _, geometry_confidence = self.get_ligand_geometry_and_isomers()
        d = {
            'geometry': geometry,                                                       # str
            'geometric_isomers_hapdent_idc': isomer_hapdent_idc,                        # list of hapdent_idc
            'geometry_rssd': rssd,                                                      # float >= 0.0
            'geometry_confidence': geometry_confidence                                  # float > 1.0
        }

        # Add to global_props so that the information is saved when the ligand is written to a file
        self.global_props['geometry'] = geometry
        self.global_props['geometry_rssd'] = rssd
        self.global_props['geometry_confidence'] = geometry_confidence

        return d

    @cached_global_props
    def geometry(self) -> str:
        """The geometry of the ligand, e.g. '2_cis'."""
        return self._geometry_and_geometrical_isomers['geometry']

    @cached_global_props
    def geometric_isomers_hapdent_idc(self) -> list:
        """A list of hapdent_idc with different orders of elements, representing the geometric isomers of the ligand."""
        return self._geometry_and_geometrical_isomers['geometric_isomers_hapdent_idc']

    @cached_global_props
    def geometry_rssd(self) -> float:
        """The RSSD value of the geometry assignment of the ligand. The lower the value, the more close the real geometry is to the ideal geometry."""
        return self._geometry_and_geometrical_isomers['geometry_rssd']

    @cached_global_props
    def geometry_confidence(self) -> float:
        """The confidence of the geometry assignment of the ligand. The higher the value, the more confident the geometry assignment is."""
        return self._geometry_and_geometrical_isomers['geometry_confidence']

    @cached_global_props
    def min_interatomic_distance(self) -> float:
        """The minimum interatomic distance in the ligand."""
        return self.get_interatomic_distances()[0]

    @cached_global_props
    def max_ligand_extension(self) -> float:
        """The maximum extension of the ligand, which is the maximum distance between any two atoms in the ligand."""
        return self.get_interatomic_distances()[1]

    @cached_global_props
    def smiles(self) -> Union[str,None]:
        """The SMILES string of the ligand."""
        return self.get_smiles()

    @cached_global_props
    def metal_counts(self) -> dict:
        """A dictionary specifying how often the ligand has been observed binding to a metal center."""
        return get_stable_sorted_value_counts(self.other_ligand_instances['parent_metal'])

    @cached_global_props
    def metal_os_counts(self) -> dict:
        """A dictionary specifying how often the ligand has been observed binding to a metal center with a specific oxidation state."""
        mos_counts = [f'{el}+{mos:.0f}' if mos > 0 else f'{el}{mos:.0f}' for el, mos in zip(self.other_ligand_instances['parent_metal'], self.other_ligand_instances['parent_metal_os']) if not np.isnan(mos)]
        return get_stable_sorted_value_counts(mos_counts)

    @cached_global_props
    def smiles_with_metal(self) -> Union[str,None]:
        """The SMILES string of the ligand with a pseudo metal center (Hg) connected to it."""
        return self.get_smiles(with_metal='Hg')

    @cached_global_props
    def is_2D_symmetrical(self) -> bool:
        """If the ligand graph is symmetrical between any two donor atoms."""
        return self._check_if_2D_symmetrical()

    @cached_global_props
    def graph_hash_with_metal(self) -> str:
        """The hash of the ligand graph with a pseudo metal center (Hg) connected to it."""
        return self.get_graph_hash_with_metal(metal_symbol='Hg')

    @cached_global_props
    def heavy_atoms_graph_hash_with_metal(self) -> str:
        """The hash of the ligand graph with a pseudo metal center (Hg) connected to it, but only considering heavy atoms."""
        return self.get_heavy_atoms_graph_hash_with_metal(metal_symbol='Hg')

    @cached_global_props
    def n_beta_hydrogens(self) -> int:
        """The number of beta-Hydrogen atoms in the ligand."""
        return self._get_n_beta_hydrogens()

    @cached_global_props
    def is_haptic(self) -> bool:
        """Whether the ligand has any haptic atoms."""
        return self.n_haptic_atoms > 0

    def get_smiles(self, with_metal: str=None) -> Union[str,None]:
        """
        Returns the SMILES string of the molecule. For around 4% of ligands in the MetaLig the bond orders are unknown, in this case it returns None.
        :param with_metal: The SMILES string is generated with the specified metal symbol connected to the ligand. If None, the ligand is not connected to any metal.
        :return: SMILES string of the molecule.
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
        Get ASE molecule with metal at original metal location.
        :param metal: Element symbol of the metal to be added to the ASE molecule.
        :return: ASE Atoms object with the metal added at the original metal position.
        """
        return self.get_ase_atoms(add_atoms=[(metal, self.parent_metal_position)])

    def has_global_property_in_range(self, name, range: list = None, values: list = None) -> bool:
        """
        Checks whether the value of the given property is within the specified range or in the specified list, if either of these are given.
        :param name: name of the property in `global_props`
        :param range: List of ranges (min, max) the property should be in. Several ranges can be specified or a single range.
        :param values: list of values the property should be in.
        :return: True if the value of the given property is within the specified range or in the specified list, False otherwise
        """
        try:
            value = self.global_props[name]
        except KeyError:
            raise ValueError(f'Property {name} is not in `global_props`.')

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

    def has_specified_stoichiometry(self, elements, instruction, only_donors: bool=False) -> bool:
        """
        Checks if the ligand has the specified stoichiometry.
        :param elements: List of chemical elements the ligand should contain.
        :param instruction: Instruction for the stoichiometry. Can be any of ['must_contain_and_only_contain', 'must_at_least_contain', 'must_exclude', 'must_only_contain_in_any_amount']
        :param only_donors: If True, only the donor atoms are considered.
        :return: True if the ligand has the specified stoichiometry, False otherwise.
        """
        atoms_of_interest = [Element(el).symbol for el in elements]
        if only_donors:
            atoms = self.donor_elements
        else:
            atoms = self.atomic_props['atoms']

        if ((sorted(list(atoms)) == sorted(
                atoms_of_interest)) and instruction == "must_contain_and_only_contain") or \
                (all(elem in list(atoms) for elem in
                     atoms_of_interest) and instruction == "must_at_least_contain") or \
                ((any(elem in list(atoms) for elem in
                      atoms_of_interest) == False) and instruction == "must_exclude") or \
                ((all(elem in atoms_of_interest for elem in
                      list(atoms))) and instruction == "must_only_contain_in_any_amount"):
            atoms_present = True
        else:
            atoms_present = False

        return atoms_present

    def has_specified_metal_centers(self, metal_centers) -> bool:
        """
        Checks if the ligand has the specified metal centers. The metal centers are checked against the original metal center of the ligand.
        :param metal_centers: List of metal centers to check for.
        :return: True if the ligand has the specified metal centers, False otherwise.
        """
        for metal_center in metal_centers:
            check_metal_center_format(metal_center)

        parent_metals = list(self.metal_counts.keys())
        parent_mos = list(self.metal_os_counts.keys())
        all_parent_metal_centers = parent_metals + parent_mos
        has_metal_centers = any(metal in all_parent_metal_centers for metal in metal_centers)

        return has_metal_centers

    def has_specified_smarts(self, smarts, should_contain, include_metal=None) -> bool:
        """
        Checks if the ligand contains the specified SMARTS pattern. If the ligand has no valid SMILES string, it always fails the check.
        :param smarts: SMARTS pattern to check for.
        :param should_contain: If True, the ligand should contain the SMARTS pattern. If False, the ligand should not contain the SMARTS pattern.
        :param include_metal: If True, the metal center is included in the SMARTS pattern. If False, the metal center is not included in the SMARTS pattern.
        :return: True if the ligand contains the specified SMARTS pattern, False otherwise.
        """
        include_metal = False if include_metal is None else include_metal
        smiles = self.smiles_with_metal if include_metal else self.smiles
        if smiles is None:  # If the ligand has no valid SMILES string, always fail the ligand.
            return False

        has_pattern = has_smarts_pattern(smarts=smarts, smiles=smiles)
        return has_pattern == should_contain

    def get_graph_with_metal(self, metal_symbol: Union[str, None]) -> tuple[nx.Graph, int]:
        """
        Returns the graph of the ligand but with the specified metal center connected to the coordinating atoms. The metal is connected to the coordinating atoms with a bond type of 1.
        :param metal_symbol: Chemical symbol of the metal center.
        :return: Graph of the ligand including the metal center and the index of the metal node in the graph.
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
        Returns the graph hash of the ligand including the metal and the bonds to the metal.
        :return: graph hash in the format '9cfe1644c35cf7f9ef3b747b268cd586'
        """
        graph_with_metal, _ = self.get_graph_with_metal(metal_symbol=metal_symbol)
        return get_graph_hash(graph_with_metal, node_attr=self.node_label)

    def get_heavy_atoms_graph_hash_with_metal(self, metal_symbol) -> str:
        """
        Returns the graph hash of the ligand including the metal and the bonds to the metal, but ignoring hydrogen atoms.
        :return: graph hash in the format '9cfe1644c35cf7f9ef3b747b268cd586'
        """
        graph_with_metal, _ = self.get_graph_with_metal(metal_symbol=metal_symbol)
        heavy_graph_with_metal = get_heavy_atoms_graph(graph_with_metal)
        return get_graph_hash(heavy_graph_with_metal, node_attr=self.node_label)

    def get_xyz_string(self, comment: str='', with_metal: bool=True) -> str:
        """
        Returns a .xyz formatted string of the ligand structure.
        :param comment: comment for the xyz file.
        :param with_metal: If True, the metal center in its original position is included in the xyz file.
        :return: String in the format of a .xyz file.
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
        Returns an ase.Atoms object containing all atoms in the ligand plus the dummy atoms plus a list of effective donor indices of this ase.Atoms object.
        :param dummy: Element symbol of the dummy atom.
        :return: Tuple of ase.Atoms object and list of effective donor indices.
        """
        return get_all_effective_ligand_atoms_with_effective_donor_indices(
                                                                            ligand_atoms=self.atoms,
                                                                            hapdent_idc=self.hapdent_idc,
                                                                            dummy=dummy
                                                                            )

    def get_isomers_effective_ligand_atoms_with_effective_donor_indices(self, dummy='Cu') -> tuple[ase.Atoms, list[list[int]]]:
        """
        This function solves the issue of making a ligand with haptic interactions into an effective ligand without haptic interactions. It will return two properties: an ase.Atoms() object of the effective ligand atoms and a list of lists, in which the outer list is each geometrical isomer and the inner list is all donor indices, either of the dentic donors or the dummy atoms of each haptic group. For a ligand without haptic interactions, these outputs are identical to `self.atoms` and `self.geometric_isomers_hapdent_idc`.
        For a ligand with haptic interaction, for each haptic group, a dummy atom is appended to the end of a copy of `self.atoms`, where the element symbol of the dummy atom is specified by the `dummy` parameter and the position is the mean of all atoms of the haptic group. Also, for a ligand with haptic interactions, the `self.geometric_isomers_hapdent_idc` is transformed so that the inner lists contain for dentic donors as before the index of this donor atom in the ase.Atoms object, while for haptic groups, the index of the dummy atom of this haptic group in the ase.Atoms object is used. For the ligand rotation in the Isomer(), this means that this ligand can be treated like a ligand without haptic interactions. However, one needs to keep in mind to remove the dummy atoms after they are not needed anymore, and then to not use the effective donor indices returned here, but the original `self.geometric_isomers_hapdent_idc` instead which still contains the haptic information.
        :param dummy: Element symbol of the dummy atom.
        :return: Tuple of effective ligand atoms and donor indices for each isomer, ready for ligand rotation.
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

    def get_ligand_geometry_and_isomers(self) -> tuple[str, list[ase.Atoms], tuple[Union[int, tuple[int]]], float, str, float]:
        """
        Returns the ligand geometry, its geometrical isomers and other related information. Handles haptic donors by replacing each haptic group with a single dummy atom.
        :return: Tuple of:
        - The assigned geometry, e.g. '2_cis' (str)
        - List of ASE Atoms objects of the best isomers
        - List of hapdent tuples for each isomer
        - The root sum of squared differences (RSSD) of the assigned geometry (float)
        - The second-best geometry (str)
        - The weight necessary for a change in geometry (float)
        """
        eff_ligand_atoms, eff_donor_idc = self.get_all_effective_ligand_atoms_with_effective_donor_indices('Cu')
        geometry, eff_isomers, eff_isomer_donor_idc, rssd, second_geometry, weight_necessary_for_change = assign_geometry(eff_ligand_atoms, eff_donor_idc)
        # Remove Cu from the isomers to get the real ligand geometry in case of haptic ligands
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

        return geometry, real_isomers, hapdent_isomer_idc, rssd, second_geometry, weight_necessary_for_change

    def get_csv_info(self, max_entries: int=5) -> dict:
        """
        Returns a flat dictionary with information about the ligand, intended for a csv file. The properties are the `global_props` plus the CSD complex IDs, metal counts, and metal oxidation state counts are included, truncated to a maximum of `max_entries` entries.
        :param max_entries: Maximum number of entries for lists with potentially many entries.
        :return: Dictionary with information about the ligand, ready to be written to a csv file or similar.
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
        Returns a dictionary representation of the ligand, including all properties and the graph, identical to the one in the MetaLig.
        :param include_graph: If True, the graph is included in the dictionary. If False, the graph is not included.
        :param copy: If True, the dictionary is deep-copied before returning. If False, the original dictionary is returned.
        :param full_global_props: If True, calculate all expected properties of the `global_props` that are not yet calculated.
        :return: Dictionary representation of the ligand.
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
        Loads the ligand from a dictionary in the format of the MetaLig.
        :param d: Dictionary containing the ligand properties.
        :param validity_check: If True, make a few consistency checks to ensure that the ligand is valid. If False, skip the validity check.
        :return: Ligand object
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
        Makes a few consistency checks to ensure that the ligand is valid.
        :return: None
        """
        super()._check_if_molecule_valid()
        if not nx.is_connected(self.graph):
            raise ValueError(f'Graph of ligand {self.unique_name} is not fully connected.')
        if not self.n_donors == self.n_denticities + self.n_haptic_atoms:
            raise ValueError(
                f'Number of donors ({self.n_donors}) does not equal number of n_denticities ({self.n_denticities}) plus n_haptic_atoms ({self.n_haptic_atoms}) in ligand {self.unique_name}.')

        return None
    
    def _sort_global_props_inplace(self) -> None:
        """
        Sorts the global properties of the ligand to increase readability and consistency.
        :return: None
        """
        self.global_props = sort_dict_as(d=self.global_props, order=ligand_global_props, strict=False)
    
    def _get_rdkit_mol_from_smiles(self, sanitize: bool=False, with_metal: str=None) -> Union['rdkit.Chem.Mol',None]:
        """
        Returns the RDKit molecule constructed from the SMILES string of the molecule. Requires an installed RDKit package.
        :param sanitize: If True, the molecule is sanitized. If False, the molecule is not sanitized.
        :param with_metal: If not None, the ligand graph is connected to the metal with the specified symbol. If None, the ligand is not connected to any metal.
        :return: RDKit molecule from the given SMILES string or None if the molecule has unknown bond orders.
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
        Returns the distance of the ligand from the metal. Can also be used to get the maximum distance of an atom in the ligand to the metal or the distances of all coordinating elements.
        :param mode: any of ['min', 'max', 'coordinating', 'all']
        :return: Returns the specified distance(s), depending on `mode`.
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
        Returns the charge of the ligand as an integer or as np.nan if the charge is not specified.
        :return: Charge of the ligand as an integer or np.nan if the charge is not specified.
        """
        if np.isnan(self.charge):
            return self.charge
        else:
            return int(self.charge)
        
    def _check_if_2D_symmetrical(self) -> bool:
        """
        Checks if the ligand graph is symmetrical between any two donors. Essentially, this checks whether the ligand graph is symmetrical under "flipping" the ligand for generating geometric isomers. However, this does not check for 3D symmetry. Often, planar ligands are 3D symmetrical if they are 2D symmetrical, but the more bulky the ligand, the more likely it is that the ligand is not 3D symmetrical even if it is 2D symmetrical.
        This function is easy to imagine for bidentate ligands, but it also works for tridentate ligands: e.g. for planar tridentate ligands, the ligand graph might be symmetrical between the outer two donors, but different for the middle donor. This will be picked up, because the function checks if the graph looks symmetrical for any two donors.
        :return: True if the ligand graph is symmetrical between donors, False otherwise.
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
        Calculates the number of beta-Hydrogen atoms. alpha-Hydrogen is ignored. beta-Hydrogen is defined as a H atom which is exactly two bonds away from a coordinating atom, i.e. three bonds away from the metal.
        :return: Number of beta-Hydrogen atoms
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
        Calculates the planarity of donors and metal center. Empirically it was found (for tridentate and tetradentate ligands) that ligands with a planarity of 0.61 or higher can be considered planar, while ligands with a planarity below 0.61 can be considered non-planar.
        :param with_metal: If True, the original metal center is included in the calculation. If False, the metal center is not included.
        :return: Planarity of the molecule as a float between 0 and 1.0. 0 means not planar at all (a sphere), 1.0 means perfectly planar.
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
        Returns a tuple with the donor indices of the ligand in which haptic groups are in sub-tuples.
        :return: Tuple of indices with haptic groups in sub-tuples.
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
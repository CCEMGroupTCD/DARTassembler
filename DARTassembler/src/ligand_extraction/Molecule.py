# standard Python packages
import warnings
from functools import cached_property
import ase
import networkx as nx
import numpy as np
from copy import deepcopy
from typing import Union
import pysmiles
from DARTassembler.src.metalig.refactor_v1_0_0 import refactor_metalig_entry_from_v1_0_0_to_v1_1_0
# some special functions which are required
from DARTassembler.src.ligand_extraction.composition import Composition
from DARTassembler.src.ligand_extraction.utilities_Molecule import get_planarity, get_denticities_and_hapticities_idc, get_isomers_effective_ligand_atoms_with_effective_donor_indices, get_all_effective_ligand_atoms_with_effective_donor_indices, format_hapdent_idc, has_smarts_pattern
from ase.visualize import view
import re
import pandas as pd

# collection of molecule objects of other packages
from ase import Atoms
from rdkit import Chem
from pysmiles import read_smiles

from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.constants.constants import metals_in_pse
from DARTassembler.src.ligand_extraction.bond_orders import graph_to_smiles
# importing own scripts
from DARTassembler.src.ligand_extraction.utilities_graph import graph_from_graph_dict, graph_to_dict_with_node_labels, view_graph, \
    unify_graph, get_sorted_atoms_and_indices_from_graph, get_reindexed_graph, find_node_in_graph_by_label, \
    get_graph_fragments, count_atoms_with_n_bonds, get_graph_hash, get_heavy_atoms_graph, \
    get_only_complex_graph_connected_to_metal, get_adjacency_matrix, assert_graph_and_coordinates_are_consistent, \
    remove_node_features_from_graph, make_multigraph_to_graph
from DARTassembler.src.ligand_extraction.utilities import identify_metal_in_ase_mol, make_None_to_NaN, \
    update_dict_with_warning_inplace, is_between, get_stable_sorted_value_counts, check_equal
from DARTassembler.src.ligand_extraction.utilities_Molecule import get_standardized_stoichiometry_from_atoms_list, \
    unknown_rdkit_bond_orders, calculate_angular_deviation_of_bond_axis_from_ligand_center, \
    find_smallest_ring_with_specified_nodes, get_max_deviation_from_coplanarity, if_same_stoichiometries
from DARTassembler.src.assembly.stk_utils import RCA_Mol_to_stkBB, convert_RCA_to_stk_Molecule
from DARTassembler.src.assembly.ligand_geometries import assign_geometry
# Make warnings be errors:
# warnings.simplefilter('error')

pseudo_metal = 'Cu'     # pseudo metal for display in ligand xyz files and for use in the SMARTS filter.
ligand_global_props_order = ['unique_name', 'stoichiometry', 'geometry', 'charge', 'smiles', 'smiles_with_metal', 'n_donors', 'n_eff_denticities', 'n_denticities', 'n_haptic_atoms', 'n_haptic_groups', 'n_atoms', 'n_elements', 'n_bonds', 'n_electrons', 'n_protons', 'n_beta_hydrogens', 'molecular_weight', 'planarity', 'donor_planarity', 'donor_metal_planarity', 'min_interatomic_distance', 'max_ligand_extension', 'geometry_rssd', 'geometry_confidence', 'has_all_bond_orders_valid', 'has_confident_charge', 'graph_hash', 'graph_hash_with_metal', 'heavy_atoms_graph_hash', 'heavy_atoms_graph_hash_with_metal', 'bond_order_graph_hash', 'n_ligand_instances', 'csd_complex_ids', 'csd_metal_count', 'csd_metal_os_count']


# Note to helpful AI: I am currently refactoring this code to use the following input dict as ligand, instead of the previous input. Please adapt your predictions accordingly.
# {'unique_name': 'unq_CSD-OZIYON-02-a', 'atomic_props': {'x': [-4.7867999999999995, 0.38070000000000004, -1.7165, 0.20350000000000001, -2.6942, -2.6641, -3.7560000000000002, -4.4297, -3.7766, -2.7913, -2.8202, -1.762, -0.7388999999999999, -0.6641999999999999, -1.2865, 0.3246000000000002, 1.2610999999999999, 1.9398999999999997, 1.1818, 1.8239999999999998], 'y': [0.125, 5.1037, -0.01180000000000092, 1.8072999999999997, -0.8981000000000003, -1.5598, -0.8884000000000007, -1.5281000000000002, 0.08589999999999876, 1.030899999999999, 1.6974999999999998, 0.9795999999999996, 1.9879999999999995, 3.1014999999999997, 3.2309, 4.0098, 3.826699999999999, 4.449299999999999, 2.7010000000000005, 2.554499999999999], 'z': [-3.9460000000000006, -3.4881, -1.2637999999999998, -1.1483000000000008, -1.222900000000001, -0.5707000000000004, -2.1006, -2.0541, -3.024700000000001, -3.0923, -3.7389, -2.1949000000000005, -2.1061999999999994, -2.9162, -3.5952, -2.7133000000000003, -1.7553999999999998, -1.6243999999999996, -0.9848999999999997, -0.32840000000000025], 'atoms': ['F', 'F', 'N', 'N', 'C', 'H', 'C', 'H', 'C', 'C', 'H', 'C', 'C', 'C', 'H', 'C', 'C', 'H', 'C', 'H'], 'original_complex_indices': [1, 2, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]}, 'global_props': {'charge': 0, 'stoichiometry': 'C10H6F2N2', 'molecular_weight': 192.16880632599998, 'n_donors': 2, 'n_atoms': 20, 'n_elements': 4, 'n_bonds': 21, 'n_electrons': 98, 'n_protons': 98, 'n_hydrogens': 6, 'n_C_H_bonds': 6, 'n_ligand_instances': 1, 'has_haptic_interactions': False, 'has_beta_hydrogens': True, 'has_all_bond_orders_valid': True, 'has_bond_orders': True, 'has_unknown_bond_orders': False, 'parent_complex_id': 'OZIYON', 'parent_metal': 'Ir', 'parent_metal_os': 3.0, 'min_interatomic_distance': 0.9291298778965196, 'max_ligand_extension': 8.745613586821682, 'graph_hash': 'd810e651de6b310aabd5ca7060829beb', 'graph_hash_with_metal': '9cfe1644c35cf7f9ef3b747b268cd586', 'heavy_atoms_graph_hash': 'f76078eb3bbe68614cc779c42ff70282', 'heavy_atoms_graph_hash_with_metal': '8d31df32a8d11ecf0b01db06d7cba93f', 'bond_order_graph_hash': '59194cf8052a23ab8b2e41804249930e'}, 'graph': {'graph': {0: {8: {'bond_type': 1}}, 1: {15: {'bond_type': 1}}, 2: {11: {'bond_type': 12}}, 3: {18: {'bond_type': 12}}, 4: {6: {'bond_type': 12}}, 5: {4: {'bond_type': 1}}, 6: {8: {'bond_type': 12}}, 7: {6: {'bond_type': 1}}, 8: {9: {'bond_type': 12}}, 9: {11: {'bond_type': 12}}, 10: {9: {'bond_type': 1}}, 11: {12: {'bond_type': 1}}, 12: {13: {'bond_type': 12}}, 13: {15: {'bond_type': 12}}, 14: {13: {'bond_type': 1}}, 15: {16: {'bond_type': 12}}, 16: {18: {'bond_type': 12}}, 17: {16: {'bond_type': 1}}, 18: {19: {'bond_type': 1}}, 19: {18: {'bond_type': 1}}}, 'node_attributes': {0: {'node_label': 'F'}, 1: {'node_label': 'F'}, 2: {'node_label': 'N'}, 3: {'node_label': 'N'}, 4: {'node_label': 'C'}, 5: {'node_label': 'H'}, 6: {'node_label': 'C'}, 7: {'node_label': 'H'}, 8: {'node_label': 'C'}, 9: {'node_label': 'C'}, 10: {'node_label': 'H'}, 11: {'node_label': 'C'}, 12: {'node_label': 'C'}, 13: {'node_label': 'C'}, 14: {'node_label': 'H'}, 15: {'node_label': 'C'}, 16: {'node_label': 'C'}, 17: {'node_label': 'H'}, 18: {'node_label': 'C'}, 19: {'node_label': 'H'}}}, 'donor_idc': [2, 3], 'ligand_instances': {'ligand_name': ['CSD-OZIYON-02-a'], 'parent_complex_id': ['OZIYON'], 'parent_complex_charge': [0], 'parent_metal': ['Ir'], 'parent_metal_os': [3.0]}, 'parent_metal_position': [0.0, 0.0, 0.0]}


class RCA_Molecule(object):
    """
    This is a base class for a molecule, which can be either a ligand or a complex.
    """

    def __init__(self,
                 # mol: Atoms = None,
                 atomic_props: dict = None,
                 global_props: dict = None,
                 graph = None,
                 # has_ligands=True,
                 # reindex_graph: bool = False,
                 # warnings: list = None,
                 # other_props: dict = None,
                 validity_check = True,
                 node_label: str = 'node_label',
                 bond_label: str = 'bond_type',
                 # **kwargs
                 ):
        # if warnings is None:
        #     warnings = []
        # if other_props is None:
        #     other_props = {}
        if atomic_props is None:
            atomic_props = {}
        if global_props is None:
            global_props = {}

        # self.warnings = warnings
        self.atomic_props = atomic_props
        self.global_props = global_props
        self.mol = self.get_ase_molecule()
        self.node_label = node_label  # node label in the graph specifying the atom type
        self.bond_label = bond_label  # edge label in the graph specifying the bond type
        self.graph = get_reindexed_graph(graph)




        # todo: refactor so that `comment` is not in the atomic_props anymore but in the global_props
        # if 'comment' in self.atomic_props:
        #     self.global_props['comment'] = ''.join(self.atomic_props['comment'])
        #     del self.atomic_props['comment']

        # Generate mol from atomic_props if possible and no mol given yet

        # if has_ligands:
        #     # if we expect ligands, we can set up an empty ligand list
        #     self.ligands = []


        # If these keys are found in global_props, we treat them as pre‐computed values:
        # precomputed_candidates = [
        #     # 'n_atoms',
        #     'n_hydrogens',
        #     'n_protons',
        #     'n_bonds',
        #     'has_bond_order_attribute',
        #     'has_unknown_bond_orders',
        #     'has_good_bond_orders',
        #     'graph_hash',
        #     'heavy_atoms_graph_hash',
        #     'bond_order_graph_hash',
        #     # 'hash',
        #     'stoichiometry'
        # ]
        # for key in precomputed_candidates:
        #     if key in self.global_props:
        #         setattr(self, f'{key}', self.global_props[key])

        # If these keys are in other_props, also use them as precomputed:
        # for key in precomputed_candidates:
        #     if key in other_props:
        #         setattr(self, f'{key}', other_props[key])

        # Set properties from the MetaLig database as attributes
        # self.set_other_props_as_properties(other_props=other_props)

        # self.add_additional_molecule_information_to_global_props()

        if validity_check:
            self.validity_check_created_molecule()

    @cached_property
    def n_atoms(self) -> int:
        return len(self.atomic_props['atoms'])

    @cached_property
    def n_hydrogens(self) -> int:
        return sum(1 for atom in self.atomic_props['atoms'] if atom == 'H')

    @cached_property
    def n_protons(self) -> int:
        return self.get_n_protons()

    @cached_property
    def n_bonds(self) -> int:
        return len(self.graph.edges)

    @cached_property
    def has_bond_order_attribute(self) -> bool:
        return self.check_for_bond_order_attribute()

    @cached_property
    def has_unknown_bond_orders(self) -> bool:
        return self.check_for_unknown_bond_orders()

    @cached_property
    def has_good_bond_orders(self) -> bool:
        return self.check_for_good_bond_orders()

    @cached_property
    def graph_hash(self) -> str:
        return get_graph_hash(self.graph)

    @cached_property
    def heavy_atoms_graph_hash(self) -> str:
        return self.get_heavy_atoms_graph_hash()

    @cached_property
    def bond_order_graph_hash(self) -> str:
        return self.get_bond_order_graph_hash()

    # @cached_property
    # def hash(self) -> str:
    #     return self.get_hash()

    @cached_property
    def stoichiometry(self) -> str:
        try:
            stoichiometry = self.global_props['stoichiometry']
        except KeyError:
            stoichiometry = self.get_standardized_stoichiometry()
            self.global_props['stoichiometry'] = stoichiometry
        return stoichiometry

    def validity_check_created_molecule(self) -> None:
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

        return

    # def add_warning(self, warning: str) -> None:
    #     """
    #     Adds a warning to the molecule.
    #     :param warning: warning to add
    #     :return: None
    #     """
    #     self.warnings.append(warning)
    #
    #     return

    def if_same_stoichiometry(self, other_stoi: str) -> bool:
        """
        Checks if the stoichiometry of the molecule is the same as the given stoichiometry.
        :param other_stoi: stoichiometry to compare to
        :return: True if the stoichiometry is the same, False otherwise
        """
        return if_same_stoichiometries(self.stoichiometry, other_stoi)

    def get_reindexed_graph(self) -> nx.Graph:
        """
        Returns the reindexed graph in which the nodes are indexed from 0 to n_atoms-1.
        :return: reindexed graph
        """
        reindexed_graph = get_reindexed_graph(self.graph)

        try:
            ligand_to_metal = self.ligand_to_metal
        except AttributeError:  # Molecule is not a ligand
            ligand_to_metal = None

        assert_graph_and_coordinates_are_consistent(
                                                    graph=reindexed_graph,
                                                    atoms=self.atomic_props['atoms'],
                                                    graph_hash=self.graph_hash,
                                                    ligand_to_metal=ligand_to_metal,
                                                    node_label=self.node_label
                                                    )
        return reindexed_graph

    def check_if_graph_and_coordinates_are_consistent(self) -> None:
        """
        Checks if graph and 3D coordinates are consistent.
        :return: None
        """
        try:
            ligand_to_metal = self.ligand_to_metal
        except AttributeError:  # Molecule is not a ligand
            ligand_to_metal = None
        assert_graph_and_coordinates_are_consistent(
                                                        graph=self.graph,
                                                        graph_hash=self.graph_hash,
                                                        atoms=self.atomic_props['atoms'],
                                                        ligand_to_metal=ligand_to_metal,
                                                        node_label=self.node_label
                                                    )

        return

    # def get_mol_from_input(self, mol):
    #     """
    #     Returns the molecule from the input.
    #     :param mol:
    #     :return:
    #     """
    #     if mol is None:
    #         coord_list_3D = [[self.atomic_props[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
    #                          enumerate(self.atomic_props["x"])]
    #         atom_list = self.atomic_props["atoms"]
    #         mol = Atoms(atom_list, positions=coord_list_3D)
    #
    #     return mol

    def get_ase_molecule(self, remove_elements: list=[], add_atoms: list=[]) -> Atoms:
        """
        Get ASE molecule from atomic properties.
        :param remove_elements: list of elements to remove from the molecule
        :param add_atoms: list of tuples of the form [(element, (x,y,z))] to add to the returned ase molecule
        """
        # Check input add_atoms:
        for el, coords in add_atoms:
            assert len(coords) == 3, f"Coordinates for element {el} are not of length 3: {coords}"
            assert isinstance(el, str), f"Element {el} is not a string"

        coord_list_3D = [[self.atomic_props[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
                         enumerate(self.atomic_props["x"])]
        atom_list = self.atomic_props["atoms"]

        # Remove specified elements
        if remove_elements:
            for i, (el, coords) in enumerate(zip(atom_list, coord_list_3D)):
                if el in remove_elements:
                    del atom_list[i]
                    del coord_list_3D[i]

        # Add specified elements
        for el, coords in add_atoms:
            atom_list.append(el)
            coord_list_3D.append(coords)

        mol = Atoms(atom_list, positions=coord_list_3D)

        return mol

    def set_other_props_as_properties(self, other_props):
        for prop, value in other_props.items():
            if hasattr(self, prop) and not check_equal(val1=self.__getattribute__(prop), val2=value):
                raise ValueError(f'Property {prop} is already set to {self.__getattribute__(prop)} and cannot be set to {value}.')
            self.__setattr__(prop, value)

        return

    def add_additional_molecule_information_to_global_props(self):
        if not all(key in self.global_props for key in ['n_atoms', 'n_elements', 'molecular_weight', 'n_C_H_bonds']):
            n_elements = len(list(np.unique(self.atomic_props['atoms'])))
            mol_weight = self.mol.get_masses().sum()
            n_C_H_bonds = self.count_C_H_bonds()
            info = {
                        'n_atoms': self.n_atoms,
                        'n_elements': n_elements,
                        'molecular_weight': mol_weight,
                        'n_C_H_bonds': n_C_H_bonds
                        }

            update_dict_with_warning_inplace(
                                                dict_to_update=self.global_props,
                                                dict_with_information=info
                                            )

        return

    def get_smiles(self) -> Union[str,None]:
        """
        Returns the SMILES string of the molecule.
        @return: SMILES string of the molecule.
        """
        if not self.check_for_good_bond_orders(): # if the molecule has unknown bond orders, we cannot calculate the SMILES
            return None

        smiles = graph_to_smiles(self.graph, element_label=self.node_label, bond_label=self.bond_label)

        return smiles

    def get_rdkit_mol_from_smiles(self, sanitize: bool=False) -> Union[Chem.Mol,None]:
        """
        Returns the RDKit molecule constructed from the SMILES string of the molecule.
        @return: RDKit molecule from the given SMILES string. Returns None if the molecule has unknown bond orders and therefore no SMILES can be calculated.
        """
        smiles = self.get_smiles()
        if smiles is None:
            return None

        return Chem.MolFromSmiles(smiles, sanitize=sanitize)

    def count_C_H_bonds(self) -> int:
        """
        Returns the number of C-H bonds in the molecule.
        @return: number of C-H bonds in the molecule.
        """
        n_bonds = 0
        atoms = self.graph.nodes(data=True)
        for idx1, idx2 in self.graph.edges():
            el1 = atoms[idx1][self.node_label]
            el2 = atoms[idx2][self.node_label]
            if sorted([el1, el2]) == ['C', 'H']:
                n_bonds += 1

        return n_bonds

    def get_n_protons(self) -> int:
        n_protons = sum([DART_Element(el).atomic_number for el in self.atomic_props['atoms']])
        return n_protons

    def has_bond_type(self, bond_types: list) -> bool:
        """
        Checks whether the molecular graph has any of the specified bond types.
        @param bond_types: list of integers of rdkit bond types
        @return: True if the molecular graph has any of the specified bond types, False otherwise
        """
        if not self.check_for_bond_order_attribute():
            return False

        for idx1, idx2, bond_dict in self.graph.edges(data=True):
            bond_type = bond_dict[self.bond_label]
            if bond_type in bond_types:
                return True

        return False

    def count_bond_types(self, bond_types: list) -> Union[int, float]:
        """
        Counts the number of specified bond types in the molecular graph.
        @param bond_types: list of integers of rdkit bond types
        @return: True if the molecular graph has any of the specified bond types, False otherwise
        """
        if not self.check_for_bond_order_attribute():
            return np.nan

        n = 0
        for idx1, idx2, bond_dict in self.graph.edges(data=True):
            bond_type = bond_dict[self.bond_label]
            if bond_type in bond_types:
                n += 1

        return n

    def check_for_bond_order_attribute(self) -> bool:
        """
        Checks if the graph has bond orders for all bonds. Note, these bond orders can be unknown to rdkit and therefore be not useful. If you want to test just for good bond orders, use `self.check_for_good_bond_orders()`.
        """
        # try:
        #     # Do not compute again if already computed
        #     contains_bond_orders = self.has_bond_order_attribute
        # except AttributeError:
        bond_orders_present = [self.bond_label in self.graph.edges[edge] for edge in self.graph.edges]
        contains_bond_orders = all(bond_orders_present)

        if not contains_bond_orders and any(bond_orders_present):
            warnings.warn('Not all bonds in the molecule have bond orders. Some bonds have bond orders, some do not. Return False in check `self.check_for_bond_order_attribute()`.')

        return contains_bond_orders

    def check_for_unknown_bond_orders(self) -> bool:
        """
        Checks whether the molecular graph has any bond orders that cannot be understood by smiles.
        @return:
        """
        # try:
        #     # Do not compute again if already computed
        #     return self.has_unknown_bond_orders
        # except AttributeError:
        #     pass

        if not self.check_for_bond_order_attribute():
            return True

        return self.has_bond_type(unknown_rdkit_bond_orders)

    def count_unknown_bond_orders(self) -> Union[int, float]:
        """
        Counts the number of bond orders which are specified as unknown.
        @return:
        """
        if not self.check_for_bond_order_attribute():
            return np.nan

        return self.count_bond_types(unknown_rdkit_bond_orders)

    def check_for_good_bond_orders(self) -> bool:
        """
        Checks whether the molecular graph has any bond orders that cannot be understood by smiles.
        @return: True if the molecular graph has any bond orders that cannot be understood by smiles, False otherwise
        """
        # try:
        #     good_bond_orders_present = self.has_good_bond_orders
        # except AttributeError:
        good_bond_orders_present = self.check_for_bond_order_attribute() and not self.check_for_unknown_bond_orders()

        return good_bond_orders_present

    def has_specified_atomic_neighbors(self, atom, neighbors) -> bool:
        """
        Checks if there's at least one node of the specified atom type that is connected
        to the specified neighbors in at least their given count, but can also be connected
        to other atom types.

        Parameters:
        - atom (str): The atom label to search for.
        - neighbors (list): A list of neighbor atom labels to check against the atom's neighbors.
        - node_label (str, optional): The node attribute key that contains the atom label. Default is 'node_label'.

        Returns:
        - bool: True if at least one specified atom is connected to all specified neighbors in at least their given count, False otherwise.
        """
        neighbor_counts = pd.Series(neighbors).value_counts().to_dict()
        for node, data in self.graph.nodes(data=True):
            if data.get(self.node_label) == atom:  # Check if the node is the specified atom
                # Count the occurrence of each neighbor atom type
                atomic_neighbors = [self.graph.nodes[neighbor].get(self.node_label) for neighbor in
                                    self.graph.neighbors(node)]
                atomic_neighbor_counts = pd.Series(atomic_neighbors).value_counts().to_dict()

                # Check if the specified neighbors are present in the required count
                all_present = True
                for required_neighbor, required_count in neighbor_counts.items():
                    if atomic_neighbor_counts.get(required_neighbor, 0) < required_count:
                        all_present = False
                        break  # A required neighbor is not present in sufficient quantity, no need to check further

                if all_present:
                    return True  # Found a node that meets the criteria

        return False  # No nodes met the criteria

    @classmethod
    def make_from_atomic_properties(cls,
                                    atomic_props_mol: dict,
                                    global_props_mol: dict,
                                    graph=None,
                                    # graph_creating_strategy: str = "default",
                                    # **kwargs
                                    ):
        """
        A more convenient creation method, as in general the atomic properties already imply the information for
        the ase mol
        """

        # coord_list_3D = [[atomic_props_mol[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
        #                  enumerate(atomic_props_mol["x"])]
        # atom_list = atomic_props_mol["atoms"]

        return cls(#mol=Atoms(atom_list, positions=coord_list_3D),
                   atomic_props=atomic_props_mol,
                   global_props=global_props_mol,
                   graph=graph,
                   # graph_creating_strategy=graph_creating_strategy,
                   # **kwargs
                   )

    def get_graph_fragments(self) -> tuple[list, list]:
        """
        Returns a list of the fragment indices (unconnected components) and their elements of the molecular graph.
        """
        indices_fragments, element_fragments = get_graph_fragments(graph=self.graph, atom_label=self.node_label)

        return indices_fragments, element_fragments

    def get_atomic_positions(self, atomic_indices) -> list:
        """
        Returns the atomic positions of the atoms with the given indices.
        """
        positions = []
        for i in atomic_indices:
            x, y, z = self.atomic_props['x'][i], self.atomic_props['y'][i], self.atomic_props['z'][i]
            positions.append([x, y, z])

        return positions

    def view_3d(self):
        view(self.mol)

    # Graph stuff
    # def make_graph(self, graph_creating_strategy: str, **kwargs):
    #     """
    #     This method also allows to overwrite graphs via the command
    #     self.graph = self.make_graph(...)
    #     which will set self.graph to the newly created graph
    #     """
    #     from DARTassembler.src.ligand_extraction.GraphCreation import GraphCreation
    #
    #     return GraphCreation(
    #         graph_creating_strategy=graph_creating_strategy,
    #         molecule=self.mol,
    #         atomic_props=self.atomic_props,
    #         **kwargs
    #     ).G

    def view_graph(self, node_size=150):
        """
        simple plot of the molecule as a graph, only connectivity, no distances
        """
        view_graph(self.graph, node_label=self.node_label, node_size=node_size)

    @cached_property
    def has_graph_hash(self):
        return False if self.graph_hash is None else True

    @cached_property
    def min_interatomic_distance(self):
        try:
            min = self.global_props['min_interatomic_distance']
        except KeyError:
            min, _, _ = self.get_atomic_distances_between_atoms()
            self.global_props['min_interatomic_distance'] = min
        return min

    @cached_property
    def max_ligand_extension(self):
        try:
            max = self.global_props['max_ligand_extension']
        except KeyError:
            _, max, _ = self.get_atomic_distances_between_atoms()
            self.global_props['max_ligand_extension'] = max
        return max

    # @cached_property
    # def has_hash(self):
    #     return False if self.hash is None else True
    #
    # @cached_property
    # def hash(self):
    #     return self.get_hash()

    # def __hash__(self):
    #     if not self.has_hash():
    #         self.get_hash()
    #     return self.hash

    # def get_graph_hash(self):
    #     self.graph_hash = get_graph_hash(self.graph)
    #     return self.graph_hash

    def get_bond_order_graph_hash(self):
        if self.check_for_good_bond_orders():
            graph_hash = get_graph_hash(self.graph, node_attr=self.node_label, edge_attr=self.bond_label)
        else:
            graph_hash = None

        return graph_hash

    def get_heavy_atoms_graph_hash(self):

        heavy_graph = get_heavy_atoms_graph(self.graph, element_label=self.node_label)
        graph_hash = get_graph_hash(heavy_graph)

        return graph_hash

    # def get_hash(self):
    #     # this hash is deterministic, unlike pythons inbuilt hash()
    #     self.hash = int(hashlib.md5(self.graph_hash.encode(encoding='UTF-8', errors='strict')).hexdigest(), 16)
    #
    #     return self.hash

    def get_standardized_stoichiometry(self) -> str:
        """
        Returns a string with the stoichiometry. Metals come first, then C, then H, then other elements. If an element exists only once, the count of 1 is not written.
        :return: stoichiometry (str)
        """
        formula = get_standardized_stoichiometry_from_atoms_list(self.get_elements_list())
        return formula


    def get_elements_list(self) -> list:
        """
        Returns a list of the elements in the molecule.
        :return: list of elements
        """
        return self.atomic_props['atoms']


    def calculate_planarity(self) -> float:
        """
        Calculates the planarity of all atoms in the molecule.
        @return: Planarity of the molecule as a float between 0 and 1. 0 means not planar at all (a 3D sphere), 1 means perfectly planar.
        """
        coordinates = self.get_coordinates_list()
        deviation = get_max_deviation_from_coplanarity(points=coordinates)  # deviation is a float that is 0 if the molecule is perfectly planar and > 0 if it is not. The higher the value, the less planar the molecule is.
        planarity = 1/ (1+ deviation)   # planarity is a float between 0 and 1. 0 means not planar at all (a sphere), 1 means perfectly planar.
        planarity = round(planarity, 10)    # round to 10 decimal places to avoid floating point errors which happen with np.linalg.svd() in different versions of numpy

        return planarity

    def get_coordinates_list(self) -> list:
        """
        Returns the coordinates of the ligand without the metal.
        @return: Coordinates of the ligand without the metal in format [[x1, y1, z1], [x2, y2, z2], ...]
        """
        coords = [[self.atomic_props['x'][i], self.atomic_props['y'][i], self.atomic_props['z'][i]] for i in range(len(self.atomic_props['x']))]

        return coords

    def get_atomic_distances_between_atoms(self, skip_elements: Union[str, list]=[]) -> float:
        """
        Returns the minimum, maximum and all distances between two atoms in the molecule.
        @param skip_elements: Do not include bonds with these elements in the calculation
        @return: tuple of (minimum, maximum, all) distance of atoms in Angstrom
        """
        if isinstance(skip_elements, str):
            skip_elements = [skip_elements]

        atoms = self.mol.get_chemical_symbols()
        valid = np.array([not (el in skip_elements) for el in atoms])

        distances = self.mol.get_all_distances()
        distances = distances[valid,:][:,valid]

        if len(distances) > 1:
            min_dist = np.where(distances > 0, distances, np.inf).min()
            max_dist = np.where(distances > 0, distances, -np.inf).max()
        else:
            min_dist = 0
            max_dist = 0

        return min_dist, max_dist, distances

    def get_all_inter_atomic_distances_as_list(self):
        """
        Returns the distances between all atoms in the molecule.
        @return: list of distances between all atoms in the molecule
        """
        distances = self.get_atomic_distances_between_atoms()[2]
        unique_distances = []
        for i in range(len(distances)):
            for j in range(i+1, len(distances)):
                unique_distances.append(distances[i,j])

        return unique_distances

    # def __eq__(self, other):
    #     if not self.stoichiometry == other.stoichiometry:
    #         return False
    #
    #     return graphs_are_equal(self.graph, other.graph)
    #
    # def __ne__(self, other):
    #     return not self.__eq__(other)

    # Finally we want to be able to turn our complex into an .xyz file
    def get_xyz_file_format_string(self, comment: str='') -> str:
        """
        returns a string that can be written into an .xyz file
        """
        str_ = f"{len(self.atomic_props['x'])}\n"
        str_ += comment + '\n'
        for i, _ in enumerate(self.atomic_props['x']):
            str_ += f"{self.atomic_props['atoms'][i]}  {self.atomic_props['x'][i]}  {self.atomic_props['y'][i]}  {self.atomic_props['z'][i]} \n"

        return str_

    def get_xyz_as_array(self) -> np.ndarray:
        """
        returns a numpy array of the xyz coordinates. Shape: (n_atoms, 3)
        """
        return np.array([self.atomic_props['x'], self.atomic_props['y'], self.atomic_props['z']]).T

    def print_to_xyz(self, path: str):
        if not path.endswith(".xyz"):
            raise ValueError("The path must end with .xyz")

        with open(path, "w+") as file:
            file.write(self.get_xyz_file_format_string())

    # helper method for the de-assembly
    def ligand_naming(self, denticity: int, ligand_list) -> (str, str):
        try:
            lig_key = f'CSD-{self.global_props['parent_complex_id']}'
            csd = self.global_props["CSD_code"]
        except KeyError:
            lig_key = f'CSD-{self.global_props['parent_complex_id']}'
            csd = self.global_props["CSD_code"]

        from DARTassembler.src.constants.constants import mini_alphabet
        j = 0
        while True:
            ligand_name = f'{lig_key}-0{denticity}-{mini_alphabet[j]}'
            if ligand_name not in [lig.name for lig in ligand_list]:
                break
            else:
                j += 1

        return ligand_name, csd

    def check_input_inherit_global_properties(self, inherit_global_properties: list) -> list:
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

    def remove_node_features_from_molecular_graphs_inplace(self, keep: list = None):
        """
        Removes all node features from all molecular graphs in the db except the node features specified in keep.
        :param keep: list of node features which will not be removed
        :return: None
        """
        if keep is None:
            keep = [self.node_label]

        remove_node_features_from_graph(graph=self.graph, keep=keep, inplace=True)

        return

    def normalize_multigraph_into_graph_inplace(self):
        """
        If self.graph is a MultiGraph, it is normalized into a Graph.
        :return: None
        """
        self.graph = make_multigraph_to_graph(self.graph)

        return

    def write_to_mol_dict(self):
        return {
            "atomic_props": self.atomic_props,
            "global_props": self.global_props,
            "graph_dict": graph_to_dict_with_node_labels(self.graph),
        }

    def append_to_file(self, key: str, writer):
        data = {'key': key, 'value': self.write_to_mol_dict()}
        writer.write(data)

        return

    @classmethod
    def read_from_mol_dict(cls,
                           dict_: dict,
                           graph_creating_strategy: str = "default",
                           **kwargs
                           ):
        assert {"atomic_props", "global_props", "graph_dict"}.issubset(set(dict_.keys()))
        if dict_["graph_dict"] is None:
            return cls.make_from_atomic_properties(atomic_props_mol=dict_["atomic_props"],
                                                   global_props_mol=dict_["global_props"],
                                                   # graph_creating_strategy=graph_creating_strategy,
                                                   **kwargs
                                                   )
        elif isinstance(dict_["graph_dict"], dict):
            return cls.make_from_atomic_properties(atomic_props_mol=dict_["atomic_props"],
                                                   global_props_mol=dict_["global_props"],
                                                   graph=graph_from_graph_dict(dict_["graph_dict"])
                                                   )
        elif isinstance(dict_["graph_dict"], nx.MultiGraph) or isinstance(dict_["graph_dict"], nx.Graph):
            return cls.make_from_atomic_properties(atomic_props_mol=dict_["atomic_props"],
                                                   global_props_mol=dict_["global_props"],
                                                   graph=unify_graph(G=dict_["graph_dict"])
                                                   )
        else:
            print("Unreadable graph format")
            return cls.make_from_atomic_properties(atomic_props_mol=dict_["atomic_props"],
                                                   global_props_mol=dict_["global_props"],
                                                   # graph_creating_strategy=graph_creating_strategy,
                                                   **kwargs
                                                   )

    def count_atoms_with_n_bonds(self, element: Union[str, None], n_bonds: int, graph_element_label: str='node_label') -> int:
        """
        Count the number of occurrences of element `element` with exactly `n_bonds` bonds.
        @param element (str, None): specification of the element, e.g. 'C'. If None, all elements are counted.
        @param n_bonds (int): count an atom if it has exactly this number of bonds
        @param graph_element_label: the label of the element string in the graph attributes. Only necessary if element is not None.
        @return (int): integer count of the occurrences
        """
        return count_atoms_with_n_bonds(graph=self.graph, element=element, n_bonds=n_bonds, graph_element_label=graph_element_label)


    def contains_only(self, atoms: Union[str, list], except_elements: list=[]) -> bool:
        """
        Returns True if the molecule contains only elements out of the list `atoms`, otherwise False.
        @param atoms: list of elements to check for
        @param except_elements: ignore these elements in the molecule when testing
        @return:
        """
        if isinstance(atoms, str):
            atoms = [atoms]

        own_atoms = [atom for atom in self.atomic_props['atoms'] if not atom in except_elements]
        contains_only_atoms = all(np.isin(own_atoms, atoms))

        return contains_only_atoms

    def to_stk_mol(self):
        return convert_RCA_to_stk_Molecule(self)

    def to_pymatMol(self):
        from pymatgen.core.structure import Molecule as PyMatMol
        return PyMatMol(species=self.atomic_props["atoms"],
                        coords=[[self.atomic_props[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
                                enumerate(self.atomic_props["x"])])

class RCA_Ligand(RCA_Molecule):
    """
    This class is used to represent ligands in the MetaLig database.
    """

    def __init__(self,
                 atomic_props: dict,
                 ligand_to_metal: list,
                 graph: nx.Graph,
                 unique_name: str,
                 global_props: dict = None,
                 ligand_instances: dict = None,
                 hapdent_idc: list = None,
                 geometric_isomers_hapdent_idc: list = None,
                 # other_props=None,
                 validity_check=True,
                 ):

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
        if not {'ligand_name', 'parent_complex_id', 'parent_complex_charge', 'parent_metal', 'parent_metal_os'}.issubset(set(ligand_instances.keys())):  # Check if all necessary keys are present
            raise ValueError('The dictionary `ligand_instances` must contain the keys `ligand_name`, `parent_complex_id`, `parent_complex_charge`, `parent_metal`, `parent_metal_os`.')

        super().__init__(#mol=Atoms(atom_list, positions=coord_list_3D),
                         atomic_props=atomic_props,
                         global_props=global_props,
                         graph=graph,
                         # has_ligands=False,
                         # other_props=other_props,
                         validity_check=validity_check,
                         # **kwargs
                         )
        coord_list_3D = [[atomic_props[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
                         enumerate(atomic_props["x"])]
        atom_list = atomic_props["atoms"]
        self.coordinates = {i: [at, coord_list_3D[i]] for i, at in enumerate(atom_list)}

        # This attribute transforms the indices of atoms in the atomic properties (starting with 0 until n-1) to the corresponding indices of self.graph. This is necessary because self.graph has indices not starting at 0 and going to n-1, but it has kept the original indices from the complex. This is a historic issue that has not been solved yet.
        self.atomic_index_to_graph_index = {atm_idx: graph_idx for atm_idx, graph_idx in enumerate(sorted(self.graph.nodes))}
        self.graph_index_to_atomic_index = {graph_idx: atm_idx for atm_idx, graph_idx in self.atomic_index_to_graph_index.items()}

        # self.denticity = denticity # integer denticity, -1 for unconnected
        # self.name = name # str name of ligand
        self.unique_name = unique_name
        self.pred_charge = self.global_props['charge']
        self.other_ligand_instances = ligand_instances

        # Saving the hapdent_idc as json converts the tuples to lists, so we need to convert them back to tuples
        if hapdent_idc is not None:
            self.hapdent_idc = format_hapdent_idc(hapdent_idc)
        if geometric_isomers_hapdent_idc is not None:
            self.geometric_isomers_hapdent_idc = [format_hapdent_idc(isomer_idc) for isomer_idc in geometric_isomers_hapdent_idc]

        # the indices and elements where the ligands was bound to the metal
        self.ligand_to_metal = ligand_to_metal
        self.local_elements = self.get_local_elements()
        self.pred_charge_is_confident = self.global_props['has_confident_charge']
        self.occurrences = self.global_props['n_ligand_instances']

        # Calculate planarity of the ligand and its donors.
        self.donor_metal_planarity = self.calculate_donors_planarity(with_metal=True)
        self.donor_planarity = self.calculate_donors_planarity(with_metal=False)
        self.planarity = self.calculate_planarity()
        self.global_props['planarity'] = self.planarity
        self.global_props['donor_planarity'] = self.donor_planarity
        self.global_props['donor_metal_planarity'] = self.donor_metal_planarity

        # Calculate denticity and hapticity of the ligand
        self.denticity = len(self.ligand_to_metal)
        assert self.denticity == self.n_denticities + self.n_haptic_atoms, f'Number of donors ({self.denticity}) does not equal number of n_denticities ({self.n_denticities}) plus n_haptic_atoms ({self.n_haptic_atoms}) in ligand {self.unique_name}.'

        # Mention some properties so it's certain they are stored in global_props and computed if they don't exist yet.
        self.n_denticities
        self.n_haptic_atoms
        self.n_eff_denticities
        self.n_haptic_groups
        self.geometry
        self.stoichiometry
        self.n_beta_hydrogens
        self.smiles
        self.smiles_with_metal
        self.min_interatomic_distance
        self.max_ligand_extension


        # self.was_connected_to_metal = len(self.local_elements) > 0
        #
        # if "csd_code" in kwargs.keys():
        #     self.csd_code = kwargs['csd_code']
        #
        # if "unique_name" in kwargs:
        #     self.unique_name = kwargs["unique_name"]
        #
        # if 'original_metal' in kwargs.keys():
        #     self.original_metal = kwargs['original_metal']
        # if 'original_metal_position' in kwargs.keys():
        #     self.original_metal_position = kwargs['original_metal_position']
        # try:
        #     # Only the original ligands had this property, it's deleted for unique ligands and cannot be used.
        #     self.original_metal_symbol = DART_Element(self.original_metal).symbol
        # except (ValueError, AttributeError):
        #     pass
        #
        # if 'original_metal_os' in kwargs.keys():
        #     self.original_metal_os = kwargs['original_metal_os']

        # if not hasattr(self, 'is_centrosymmetric') or not hasattr(self, 'centrosymmetry_ang_dev'):
        #     self.is_centrosymmetric, self.centrosymmetry_ang_dev = self.check_if_centrosymmetric(return_ang_dev=True)
        # if not hasattr(self, 'graph_hash_with_metal'):
        #     # The graph_hash_with_metal is not calculated with the real original metal, but with a pseudo metal which is always the same, so that only the connections to the metal matter, but different original metals won't give different hashes. This means ligands are considered the same independent of the original metal under this hash.
        #     self.graph_hash_with_metal = self.get_graph_hash_with_metal(metal_symbol='Hg')
        # if not hasattr(self, 'heavy_atoms_graph_hash_with_metal'):
        #     self.heavy_atoms_graph_hash_with_metal = self.get_heavy_atoms_graph_hash_with_metal(metal_symbol='Hg')
        # if not hasattr(self, 'has_betaH'):
        #     self.has_betaH = self.betaH_check()
        # if not hasattr(self, 'has_neighboring_coordinating_atoms'):
        #     self.has_neighboring_coordinating_atoms = self.check_for_neighboring_coordinating_atoms()
        # if not hasattr(self, 'stats'):
        #     self.stats = self.get_ligand_stats()

        assert nx.is_connected(self.graph), f'Graph of ligand with name {self.unique_name} is not fully connected.'
        self._sort_global_props_inplace()

    @cached_property
    def original_complex_id(self):
        return self.other_ligand_instances['parent_complex_id'][0]

    @cached_property
    def original_metal_position(self):
        return self.other_ligand_instances['parent_metal_position'][0]

    @cached_property
    def parent_metal(self):
        return self.other_ligand_instances['parent_metal'][0]

    @cached_property
    def hapdent_idc(self):
        return self.get_denticities_and_hapticities_idc()

    @cached_property
    def n_eff_denticities(self):
        try:
            n_eff_denticities = self.global_props['n_eff_denticities']
        except KeyError:
            n_eff_denticities = len(self.hapdent_idc)   # effective ligand coordination number
            self.global_props['n_eff_denticities'] = n_eff_denticities
        return n_eff_denticities

    @cached_property
    def n_denticities(self):
        try:
            n_denticities = self.global_props['n_denticities']
        except KeyError:
            n_denticities = sum([1 for el in self.hapdent_idc if isinstance(el, int)])
            self.global_props['n_denticities'] = n_denticities
        return n_denticities

    @cached_property
    def n_haptic_atoms(self):
        try:
            n_haptic_atoms = self.global_props['n_haptic_atoms']
        except KeyError:
            n_haptic_atoms = sum([len(sublist) for sublist in self.hapdent_idc if isinstance(sublist, tuple)])
            self.global_props['n_haptic_atoms'] = n_haptic_atoms
        return n_haptic_atoms

    @cached_property
    def n_haptic_groups(self):
        try:
            n_haptic_groups = self.global_props['n_haptic_groups']
        except KeyError:
            n_haptic_groups = len([sublist for sublist in self.hapdent_idc if isinstance(sublist, tuple)])
            self.global_props['n_haptic_groups'] = n_haptic_groups
        return n_haptic_groups

    @cached_property
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

    @cached_property
    def geometry(self):
        try:
            return self.global_props['geometry']
        except KeyError:
            return self._geometry_and_geometrical_isomers['geometry']

    @cached_property
    def geometric_isomers_hapdent_idc(self):
        return self._geometry_and_geometrical_isomers['geometric_isomers_hapdent_idc']

    @cached_property
    def geometry_rssd(self):
        try:
            return self.global_props['geometry_rssd']
        except KeyError:
            return self._geometry_and_geometrical_isomers['geometry_rssd']

    @cached_property
    def geometry_confidence(self):
        try:
            return self.global_props['geometry_confidence']
        except KeyError:
            return self._geometry_and_geometrical_isomers['geometry_confidence']

    @cached_property
    def smiles(self):
        try:
            smiles = self.global_props['smiles']
        except KeyError:
            smiles = self.get_smiles()
            self.global_props['smiles'] = smiles
        return smiles

    @cached_property
    def count_metals(self):
        # Sort by count so that the most common metal is first
        return get_stable_sorted_value_counts(self.other_ligand_instances['parent_metal'])

    @cached_property
    def mos_counts(self):
        mos_counts = [f'{el}+{mos:.0f}' if mos > 0 else f'{el}{mos:.0f}' for el, mos in zip(self.other_ligand_instances['parent_metal'], self.other_ligand_instances['parent_metal_os']) if not np.isnan(mos)]
        # Sort by count so that the most common metal-os combination is first
        return get_stable_sorted_value_counts(mos_counts)

    @cached_property
    def smiles_with_metal(self):
        try:
            smiles = self.global_props['smiles_with_metal']
        except KeyError:
            smiles = self.get_smiles(with_metal='Hg')
            self.global_props['smiles_with_metal'] = smiles
        return smiles

    @cached_property
    def is_centrosymmetric(self):
        return self.check_if_centrosymmetric()

    @cached_property
    def centrosymmetry_ang_dev(self):
        return self.calculate_angular_deviation_from_centrosymmetry()

    @cached_property
    def graph_hash_with_metal(self):
        return self.get_graph_hash_with_metal(metal_symbol='Hg')

    @cached_property
    def heavy_atoms_graph_hash_with_metal(self):
        return self.get_heavy_atoms_graph_hash_with_metal(metal_symbol='Hg')

    @cached_property
    def has_betaH(self):
        return self.n_beta_hydrogens > 0

    @cached_property
    def n_beta_hydrogens(self):
        try:
            n = self.global_props['n_beta_hydrogens']
        except KeyError:
            n = self.get_n_beta_hydrogens()
            self.global_props['n_beta_hydrogens'] = n
        return n

    @cached_property
    def has_neighboring_coordinating_atoms(self):
        return self.n_haptic_atoms > 0

    def get_smiles(self, with_metal: str=None) -> Union[str,None]:
        """
        Returns the SMILES string of the molecule.
        @param with_metal: If not None, the ligand graph is connected to the metal with the specified symbol.
        @return: SMILES string of the molecule.
        """
        if not self.check_for_good_bond_orders(): # if the molecule has unknown bond orders, we cannot calculate the SMILES
            return None

        if with_metal is None:
            graph = self.graph
        else:
            if not DART_Element(with_metal).is_metal:
                raise ValueError(f'Invalid input for with_metal: {with_metal}. Must be a metal symbol.')

            graph = self.get_graph_with_metal(metal_symbol=with_metal)

        smiles = graph_to_smiles(graph, element_label=self.node_label, bond_label=self.bond_label)

        return smiles

    def get_rdkit_mol_from_smiles(self, sanitize: bool=False, with_metal: str=None) -> Union[Chem.Mol,None]:
        """
        Returns the RDKit molecule constructed from the SMILES string of the molecule.
        @param sanitize: If True, the molecule is sanitized. If False, the molecule is not sanitized.
        @param with_metal: If not None, the ligand graph is connected to the metal with the specified symbol. If 'original', the original metal symbol is used. If None, the ligand is not connected to any metal.
        @return: RDKit molecule from the given SMILES string. Returns None if the molecule has unknown bond orders and therefore no SMILES can be calculated.
        """
        smiles = self.get_smiles(with_metal=with_metal)
        if smiles is None:
            return None

        return Chem.MolFromSmiles(smiles, sanitize=sanitize)

    # def get_ligand_stats(self) -> dict:
    #     """
    #     Returns a dictionary with potentially interesting ligand statistics.
    #     @return: Dictionary with stats
    #     """
    #
    #     min_distance_to_metal, max_distance_to_metal, coordinating_atom_distances = self.get_atomic_distance_to_original_metal(mode='all')
    #     min_dist, max_dist, _ = self.get_atomic_distances_between_atoms()
    #     max_dist_per_atoms = max_dist / len(self.atomic_props['atoms'])
    #     stats = {
    #                 'min_distance_to_metal': min_distance_to_metal,
    #                 'max_distance_to_metal': max_distance_to_metal,
    #                 'coordinating_atom_distances_to_metal': coordinating_atom_distances,
    #                 'min_atomic_distance': min_dist,
    #                 'max_atomic_distance': max_dist,
    #                 'max_dist_per_atoms': max_dist_per_atoms
    #             }

        return stats

    def _sort_global_props_inplace(self) -> None:
        real_props = list(self.global_props.keys())
        ideal_order = ligand_global_props_order
        sorted_keys = sorted(real_props, key=lambda x: ideal_order.index(x) if x in ideal_order else len(ideal_order))
        sorted_props = {key: self.global_props[key] for key in sorted_keys}
        assert set(real_props) == set(sorted_props.keys()), f'Sorting failed. Real props: {real_props}, sorted props: {sorted_props.keys()}'
        self.global_props = sorted_props

        return

    def calculate_angular_deviation_from_centrosymmetry(self) -> float:
        """
        Calculate the angular deviation from centrosymmetry. Only defined for monodentate ligands.
        :return: Angular deviation from centrosymmetry for monodentates, np.nan for others.
        """
        return calculate_angular_deviation_of_bond_axis_from_ligand_center(atomic_props=self.atomic_props, donor_indices=self.ligand_to_metal, metal_position=self.original_metal_position)

    def check_if_centrosymmetric(self, threshold: float=8, return_ang_dev=False) -> Union[bool, tuple[bool, float]]:
        """
        Check if the ligand is centrosymmetric. Only defined for monodentate ligands.
        @return: True if the ligand is monodentate and centrosymmetric, False otherwise.
        """
        if self.denticity != 1:
            if return_ang_dev:
                return False, np.nan
            else:
                return False

        ang_dev = self.calculate_angular_deviation_from_centrosymmetry()
        centrosymmetric = bool(ang_dev < threshold)

        if return_ang_dev:
            return centrosymmetric, ang_dev
        else:
            return centrosymmetric

    def check_bidentate_planarity(self, tol: bool=0.7, return_deviation: bool=False) -> Union[bool, tuple[bool, float]]:
        """
        Checks if the ligand is bidentate planar. Only defined for bidentate ligands.
        @param tol: Maximum tolerance to be considered planar. Default is 0.7.
        @param return_deviation: If True, the calculated deviation is returned as well
        @return: True if the ligand is bidentate planar, False otherwise
        """
        if self.denticity != 2:
            if return_deviation:
                return False, np.nan
            else:
                return False

        graph, metal_idx = self.get_graph_with_metal(metal_symbol='Hg', return_metal_index=True)

        nodes = [metal_idx] + [self.atomic_index_to_graph_index[idx] for idx in self.ligand_to_metal]
        assert sorted([graph.nodes(data=self.node_label)[idx] for idx in nodes]) == sorted(['Hg'] + self.local_elements), nodes

        ring_indices = find_smallest_ring_with_specified_nodes(graph=graph, nodes=nodes)
        if ring_indices is None:
            is_coplanar, max_dist = False, np.nan
        else:
            ring_indices = [self.graph_index_to_atomic_index[idx] for idx in ring_indices if not idx == metal_idx]
            all_coordinates = {idx: coord for idx, coord in enumerate(self.get_coordinates_list())}
            coordinates = [self.original_metal_position] + [all_coordinates[idx] for idx in ring_indices]

            deviation = get_max_deviation_from_coplanarity(points=coordinates)
            is_coplanar = deviation < tol

        if return_deviation:
            return is_coplanar, deviation
        else:
            return is_coplanar

    def is_2D_symmetrical(self) -> bool:
        """
        Checks if the ligand graph is symmetrical between donors. Essentially, this checks whether the ligand graph is symmetrical under "flipping" the ligand for generating geometric isomers. However, this does not check for 3D symmetry. Often, planar ligands are 3D symmetrical if they are 2D symmetrical, but the more bulky the ligand, the more likely it is that the ligand is not 3D symmetrical even if it is 2D symmetrical.
        This function is easy to imagine for bidentate ligands, but it also works for tridentate ligands: e.g. for planar tridentate ligands, the ligand graph might be symmetrical between the outer two donors, but different for the middle donor. This will be picked up, because the function checks if the graph looks symmetrical for any two donors.
        :return: True if the ligand graph is symmetrical between donors, False otherwise.
        """
        # Explanation of algorithm: This function checks if the ligands graph is symmetrical between any of the donors. Essentially, it attaches a pseudo Hg atom to each donor and checks if the resulting ligand graph between any donors is identical. If it is, the ligand is symmetrical under these donors.

        # Get graph where all donors are connected to a pseudo Hg atom.
        graph, metal_idx = self.get_graph_with_metal(metal_symbol='Hg', return_metal_index=True)

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

    def is_mer_tridentate(self) -> bool:
        return (self.denticity == 3) and (self.if_donors_planar(with_metal=True)) and (not self.has_neighboring_coordinating_atoms)

    def is_mer_tetradentate(self) -> bool:
        return (self.denticity == 4) and (self.if_donors_planar(with_metal=True)) and (not self.has_neighboring_coordinating_atoms)

    def count_atoms_with_n_bonds(self, element: Union[str, None], n_bonds: int, remember_metal: bool=False) -> int:
        """
        Count the number of occurrences of element `element` with exactly `n_bonds` bonds.
        @param element (str, None): specification of the element, e.g. 'C'. If None, all elements are counted.
        @param n_bonds (int): count an atom if it has exactly this number of bonds
        @param graph_element_label (str): the label of the element string in the graph attributes. Only necessary if element is not None.
        @param remember_metal (bool): If the original bonds to the metal should be considered or not. Only relevant for ligands.
        @return (int): integer count of the occurrences
        """
        try:
            graph = self.graph_with_metal if remember_metal else self.graph
        except AttributeError:
            graph = self.graph

        n = count_atoms_with_n_bonds(graph=graph, element=element, n_bonds=n_bonds, graph_element_label=self.node_label)

        return n

    def get_atomic_distance(self, mode: str= 'min'):
        """
        Returns the distance of atoms to each other.
        @param mode: any of ['min', 'max', 'coordinating', 'all']
        @return: Returns the specified distance(s)
        """
        distances = np.linalg.norm(self.mol.positions - np.array(self.original_metal_position), axis=1)

        mode = mode.lower()
        if mode == 'min':
            return distances.min()
        elif mode == 'max':
            return distances.max()
        elif mode == 'coordinating':
            return distances[self.ligand_to_metal].tolist()
        elif mode == 'all':
            return distances.min(), distances.max(), distances[self.ligand_to_metal].tolist()

    def is_in_global_props(self, name, range: list = None, values: list = None) -> bool:
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

        # Always exclude NaN values.
        # if np.isnan(value):
        #     return False

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
        atoms_of_interest = [DART_Element(el).symbol for el in elements]
        if only_donors:
            atoms = self.local_elements
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
        return any([metal in list(self.count_metals.keys()) for metal in metal_centers])

    def has_specified_smarts(self, smarts, should_contain, include_metal=None) -> bool:
        """
        Checks if the ligand contains the specified SMARTS pattern. If the ligand has no valid SMILES string, it always fails the check.
        :param smarts: SMARTS pattern to check for.
        :param should_contain: If True, the ligand should contain the SMARTS pattern. If False, the ligand should not contain the SMARTS pattern.
        :param include_metal: If True, the metal center is included in the SMARTS pattern. If False, the metal center is not included in the SMARTS pattern.
        :return:
        """
        include_metal = False if include_metal is None else include_metal
        smiles = self.smiles_with_metal if include_metal else self.smiles
        if smiles is None:  # If the ligand has no valid SMILES string, always fail the ligand.
            return False

        has_pattern = has_smarts_pattern(smarts=smarts, smiles=smiles)
        return has_pattern == should_contain

    def get_atomic_distance_to_original_metal(self, mode: str= 'min'):
        """
        Returns the distance of the ligand from the metal. Can also be used to get the maximum distance of an atom in the ligand to the metal or the distances of all coordinating elements.
        @param mode: any of ['min', 'max', 'coordinating', 'all']
        @return: Returns the specified distance(s)
        """
        distances = np.linalg.norm(self.mol.positions - np.array(self.original_metal_position), axis=1)

        mode = mode.lower()
        if mode == 'min':
            return distances.min()
        elif mode == 'max':
            return distances.max()
        elif mode == 'coordinating':
            return distances[self.ligand_to_metal].tolist()
        elif mode == 'all':
            return distances.min(), distances.max(), distances[self.ligand_to_metal].tolist()

    def get_coordinates_with_metal(self):
        """
        Returns the coordinates of the ligand with the metal.
        @return: Coordinates of the ligand with the metal
        """
        return [self.original_metal_position] + self.get_coordinates_list()

    def get_graph_with_metal(self, metal_symbol: Union[str, None]=None, return_metal_index: bool=False) -> Union[nx.Graph, tuple[nx.Graph, int]]:
        """
        Returns the graph of the ligand but with the specified metal center connected to the coordinating atoms. The metal is connected to the coordinating atoms with a bond type of 1.
        @param metal_symbol: Symbol of the metal.
        @return: Graph of the ligand with the metal. If return_metal_index is True, the metal index is also returned.
        """
        graph_with_metal = nx.Graph(self.graph)   # unfreeze graph

        # Add metal node and bonds of coordinating atoms
        metal_idx = max(self.graph.nodes) + 1
        graph_with_metal.add_node(metal_idx, node_label=metal_symbol)       # hardcode: name 'node_label' is used for the element string
        for atom_idx in self.ligand_to_metal:
            # Indices of graph and atomic indices don't match
            graph_idx = self.atomic_index_to_graph_index[atom_idx]
            graph_with_metal.add_edge(metal_idx, graph_idx, bond_type=1)    # hardcode: name 'bond_type' is used for the bond type

        if return_metal_index:
            return graph_with_metal, metal_idx
        else:
            return graph_with_metal

    def get_graph_hash_with_metal(self, metal_symbol) -> str:
        """
        Returns the graph hash of the ligand including the metal and the bonds to the metal.
        @return: graph hash
        """
        graph_with_metal = self.get_graph_with_metal(metal_symbol=metal_symbol)
        return get_graph_hash(graph_with_metal, node_attr=self.node_label)

    def get_heavy_atoms_graph_hash_with_metal(self, metal_symbol) -> str:
        """
        Returns the graph hash of the ligand including the metal and the bonds to the metal. Only heavy atoms are considered.
        @return: graph hash
        """
        graph_with_metal = self.get_graph_with_metal(metal_symbol=metal_symbol)
        heavy_graph_with_metal = get_heavy_atoms_graph(graph_with_metal)
        return get_graph_hash(heavy_graph_with_metal, node_attr=self.node_label)

    def get_local_elements(self) -> list:
        """
        Calculates the elements in the first coordination sphere from `self.ligand_to_metal`.
        :return: list of elements in first coordination sphere
        """
        return [self.atomic_props["atoms"][i] for i in self.ligand_to_metal]

    def get_donor_positions(self) -> np.array:
        """
        Returns the positions of the donor atoms.
        :return: np.array of donor positions
        """
        return np.array([[self.atomic_props["x"][i], self.atomic_props["y"][i], self.atomic_props["z"][i]]
                         for i in self.ligand_to_metal])

    def get_charge_as_int(self) -> Union[float, int]:
        """
        Returns the charge of the ligand as an integer.
        """
        if not hasattr(self, 'pred_charge'):
            raise AttributeError('The charge of the ligand has not been predicted yet.')
        elif np.isnan(self.pred_charge):
            return self.pred_charge
        else:
            return int(self.pred_charge)

    def get_xyz_file_format_string(self, comment: str='', with_metal:bool=True) -> str:
        """
        Returns a string that can be written into an .xyz file.
        @comment: comment for the xyz file.
        @param with_metal: If True, the metal center in it's original position is included in the xyz file, otherwise not.
        """
        if comment is None: # default comment specifying important properties of the ligand
            donors = '-'.join(self.local_elements)
            comment = f'Ligand ID: {self.unique_name}  ===  Stoichiometry: {self.stoichiometry}  ===  Charge: {self.get_charge_as_int()}  ===  Donors: {donors}'

        n_ligand_atoms = len(self.atomic_props['x'])
        if with_metal:
            str_ = f"{n_ligand_atoms+1}\n" # +1 for the metal
            str_ += comment + '\n'
            str_ += f"{pseudo_metal}  {self.original_metal_position[0]}  {self.original_metal_position[1]}  {self.original_metal_position[2]} \n"     # metal atom
        else:
            str_ = f"{n_ligand_atoms}\n"
            str_ += comment + '\n'

        # Add ligand atoms
        for i, _ in enumerate(self.atomic_props['x']):
            str_ += f"{self.atomic_props['atoms'][i]}  {self.atomic_props['x'][i]}  {self.atomic_props['y'][i]}  {self.atomic_props['z'][i]} \n"

        return str_

    def get_n_beta_hydrogens(self) -> int:
        """
        Calculates the number of beta-Hydrogen atoms. Alpha H is ignored. Beta H is defined as a H atom which is exactly two bonds away from a coordinating atom, i.e. three bonds away from the metal.
        @return: Number of beta Hydrogen atoms
        """
        graph = self.get_reindexed_graph()     # historical issue
        A = get_adjacency_matrix(graph)
        # The second power of the adjacency matrix, i.e. A^2[i,j] represents the number of paths of length two from i to j. Hence, as we are only interested in Hydrogen which has a distance of two to our coordinating atoms.
        B = np.matmul(A, A)

        beta_h_indices = set()  # Using the set we avoid double counting of beta H atoms
        for donor_idx in self.ligand_to_metal:
            for idx, element in graph.nodes(data=self.node_label):
                if element == "H":
                    # search for beta H while excluding alpha H
                    if B[donor_idx, idx] > 0 and A[donor_idx, idx] == 0:
                        beta_h_indices.add(idx)

        n_beta_hydrogens = len(beta_h_indices)

        return n_beta_hydrogens

    def betaH_check(self) -> bool:
        """
        Checks if the ligand has beta Hydrogen. Alpha H is ignored. Beta H is defined as a H atom which is exactly two bonds away from a coordinating atom, i.e. three bonds away from the metal.
        @return: True if beta Hydrogen is present, False otherwise
        """
        graph = self.get_reindexed_graph()     # historical issue
        A = get_adjacency_matrix(graph)

        # The second power of the adjacency matrix, i.e. A^2[i,j] represents the number of paths of length two
        # from i to j. Hence, as we are only interested in Hydrogen which has a distance of two to our coordinating atoms
        B = np.matmul(A, A)

        for functional_index in self.ligand_to_metal:
            for index, atom_symbol in graph.nodes(data=self.node_label):
                # search for beta Hydrogen while excluding alpha Hydrogen
                if B[functional_index, index] > 0 and A[functional_index, index] == 0 and atom_symbol == "H":
                    # self.view_3d()
                    return True

        return False

    def check_for_neighboring_coordinating_atoms(self) -> bool:
        """
        Checks if any of the coordinating atoms are neighbors.
        @return: True if two coordinating atoms are neighbors, False otherwise
        """
        for i in self.ligand_to_metal:
            for j in self.ligand_to_metal:
                # Indices of graph and atomic indices don't match historically
                graph_index_i = self.atomic_index_to_graph_index[i]
                graph_index_j = self.atomic_index_to_graph_index[j]
                if i != j and self.graph.has_edge(graph_index_i, graph_index_j):
                    return True

        return False

    def get_denticities_and_hapticities_idc(self) -> tuple[Union[int, tuple[int]]]:
        """
        Returns a tuple with the donor indicesof the ligand in which haptic groups are in sub-tuples.
        :return: Tuple of indices with haptic groups in sub-tuples.
        """
        graph_indices = [self.atomic_index_to_graph_index[i] for i in self.ligand_to_metal]
        graph = self.graph
        hapdent_idc = get_denticities_and_hapticities_idc(graph=graph, donor_idc=graph_indices)

        # Convert back to atomic indices. Keep in mind that some entries are single indices and some are tuples of indices.
        atomic_hapdent_idc = []
        for idc in hapdent_idc:
            if isinstance(idc, int):    # denticity, therefore integer
                idc = self.graph_index_to_atomic_index[idc]
                assert idc in self.ligand_to_metal, f'Index {idc} is not in the ligand to metal indices: {self.ligand_to_metal}'
            elif isinstance(idc, tuple):
                idc = tuple([self.graph_index_to_atomic_index[i] for i in idc])
                assert all([i in self.ligand_to_metal for i in idc]), f'Indices {idc} is not in the ligand to metal indices: {self.ligand_to_metal}'
            else:
                raise ValueError(f'Invalid type of hapdent_idc: {type(idc)}')
            atomic_hapdent_idc.append(idc)

        return tuple(atomic_hapdent_idc)

    def get_effective_donor_atoms(self, dummy='Cu', only_haptic=False) -> ase.Atoms:
        """
        Returns an ase.Atoms with all effective donor atoms: if an atom is non-haptic, it is returned as is. If a group of donor atoms is haptic, the entire haptic group is replaced by a single dummy atom at the mean position of the haptic group.
        :param only_haptic: If True, only dummy atoms representing haptic groups are returned. If False, all effective donor atoms are returned, the dummy atoms and the normal non-haptic atoms.
        :param dummy: Element symbol of the dummy atom.
        :return: ase.Atoms object with effective donor atoms of the same length as hapdent_idc.
        """
        return get_effective_donor_atoms(
                                            hapdent_idc=self.hapdent_idc,
                                            ligand_atoms=self.mol,
                                            dummy=dummy,
                                            only_haptic=only_haptic,
                                            )

    def get_all_effective_ligand_atoms_with_effective_donor_indices(self, dummy='Cu') -> tuple[ase.Atoms, list[int]]:
        """
        Returns an ase.Atoms object containing all atoms in the ligand plus the dummy atoms plus a list of effective donor indices of this ase.Atoms object.
        :param dummy: Element symbol of the dummy atom.
        :return: Tuple of ase.Atoms object and list of effective donor indices.
        """
        return get_all_effective_ligand_atoms_with_effective_donor_indices(
                                                                            ligand_atoms=self.mol,
                                                                            hapdent_idc=self.hapdent_idc,
                                                                            dummy='Cu'
                                                                            )

    def get_isomers_effective_ligand_atoms_with_effective_donor_indices(self, dummy='Cu') -> tuple[ase.Atoms, list[list[int]]]:
        """
        For each geometric isomer of this ligand, returns the effective donor atoms in which a dummy donor atom is placed at the mean position of each haptic group. Also returns the effective donor indices of each isomer. In total, the resulting ase.Atoms and indices can be used like any other ligand without haptic interactions.
        :param dummy: Element symbol of the dummy atom.
        :return: Tuple of ase.Atoms object and list of effective donor indices for each isomer.
        """
        eff_atoms, isomers_eff_donor_idc = get_isomers_effective_ligand_atoms_with_effective_donor_indices(
                                                    ligand_atoms=self.mol,
                                                    geometric_isomers_hapdent_idc=self.geometric_isomers_hapdent_idc,
                                                    dummy=dummy
                                                    )
        # Doublechecking that the effective atoms are the same as the original atoms if there are no haptic interactions
        if not self.has_neighboring_coordinating_atoms:
            assert eff_atoms == self.mol, f'Error in isomer generation for ligand {self.unique_name}.'

        return eff_atoms, isomers_eff_donor_idc

    # def sort_atomic_props_to_have_coordinating_atoms_first(self):
    #     """
    #     Sorts the atomic properties such that the atoms that are coordinating are first in the list.
    #     Attention: This does not change anything else, only the atomic properties. The graph and the ligand_to_metal are still the same!
    #     """
    #     for key in self.atomic_props.keys():
    #         self.atomic_props[key] = [self.atomic_props[key][i] for i in self.ligand_to_metal] + [
    #             self.atomic_props[key][i] for i in range(len(self.atomic_props[key])) if i not in self.ligand_to_metal]
    #
    #     return

    def get_assembly_dict(self):
        """
        :return: {index: list, type: list, xyz_str: str}
        """
        dict_ = {}
        dict_["index"] = self.ligand_to_metal
        dict_["type"] = [self.atomic_props["atoms"][i] for i in self.ligand_to_metal]
        dict_["str"] = self.get_xyz_file_format_string()

        return dict_

    # def add_atom(self, symbol: str, coordinates: list[float]):
    #     if len(coordinates) != 3:
    #         print("Wrong number of coordinates specified")
    #         raise ValueError
    #
    #     self.atomic_props["atoms"].append(symbol)
    #     self.atomic_props["x"].append(coordinates[0])
    #     self.atomic_props["y"].append(coordinates[1])
    #     self.atomic_props["z"].append(coordinates[2])
    #
    #     # now we need to update the self.mol, graph wont be updated
    #     self.print_to_xyz(path="../tmp/tmp.xyz")
    #     self.mol = io.read("../tmp/tmp.xyz")

    # def functional_atom_check(self, atoms_of_interest: [str, list[str]] = None):
    #     if atoms_of_interest is None:
    #         atoms_of_interest = ["N", "O"]
    #     elif isinstance(atoms_of_interest, str):
    #         atoms_of_interest = [atoms_of_interest]
    #     return set(self.get_assembly_dict()["type"]).issubset(set(atoms_of_interest))

    # def planar_check(self, eps=2):
    #     """
    #     Checks if the ligand has planar donor atoms. Check is done only for denticity 3 and 4, otherwise False is returned.
    #     :param eps: durch try'n'error obtained
    #     eps für (d=4) -> 1
    #     :return:
    #     """
    #     # Deprecated because this method has a lot of issues. For the tridentate, the calculation is ok but not perfect and for the tetradentate it does not include the metal.
    #     raise DeprecationWarning("This method is deprecated. Use if_donors_planar() instead.")
    #
    #     if self.denticity <= 2:     # ligands with 2 or less coordinating atoms are per definition planar
    #         return True
    #     elif self.denticity >= 5:   # for 5 or more coordinating atoms, the planarity is not yet implemented
    #         return False
    #
    #     functional_coords = [[self.atomic_props[key][i] for key in ["x", "y", "z"]] for i in self.ligand_to_metal]
    #     assert len(functional_coords) == self.denticity, f"Error in Planar Check for ligand {self.name}"
    #
    #     if self.denticity == 3:
    #         c1, c2, c3 = Point3D(functional_coords[0]), Point3D(functional_coords[1]), Point3D(functional_coords[2])
    #         E = Plane(c1, c2, Point3D([0, 0, 0]))
    #         if round(E.distance(c3)) < eps:
    #             return True
    #
    #     if self.denticity == 4:
    #         c1, c2, c3, c4 = Point3D(functional_coords[0]), Point3D(functional_coords[1]), Point3D(
    #             functional_coords[2]), Point3D(functional_coords[3])
    #         E = Plane(c1, c2, c3)
    #         if round(E.distance(c4)) < eps:
    #             return True
    #
    #     return False

    def calculate_donors_planarity(self, with_metal: bool=True) -> float:
        """
        Calculates the planarity of donors and metal center.
        @param with_metal: If True, the original metal center is included in the calculation. If False, the metal center is not included.
        @return: Planarity of the molecule as a float between 0 and 1. 0 means not planar at all (a sphere), 1 means perfectly planar.
        """
        # Get the coordinates of the donors and the metal
        coordinates = self.get_coordinates_list()
        coordinates = [coordinates[donor_idx] for donor_idx in self.ligand_to_metal]
        if with_metal:
            coordinates.append(self.original_metal_position)

        # If there are less than 3 atoms, the molecule is planar by definition
        if len(coordinates) < 3:
            return 1.0

        planarity = get_planarity(coordinates)

        return planarity

    def if_donors_planar(self, threshold: float=0.61, with_metal: bool=True) -> bool:
        """
        Checks if the donors are planar.
        @param threshold: Threshold for the planarity. If the planarity is above the threshold, the donors are considered planar. By default, the threshold is set to 0.61 which has empirically been found to be a good value for tridentate and tetradentate ligands.
        @param with_original_metal: If True, the original metal center is included in the calculation. If False, the metal center is not included.
        @return: True if the donors are planar, False otherwise.
        """
        planarity = self.calculate_donors_planarity(with_metal=with_metal)
        is_planar = planarity >= threshold
        is_planar = bool(is_planar)

        return is_planar

    def get_ligand_geometry_and_isomers(self) -> tuple[str, list[ase.Atoms], tuple[Union[int, tuple[int]]], float, str, float]:
        """
        Returns the geometry of the ligand, the best isomers and some other potentially interesting information. This function is the go-to function for the ligand geometry assignment, even able to handle ligands with haptic interactions.
        :return: Tuple of:
        - The assigned geometry
        - List of ASE Atoms objects of the best isomers
        - List of hapdent tuples for each isomer
        - The root sum of squared differences (RSSD) of the assigned geometry
        - The second best geometry
        - The weight necessary for a change in geometry
        """
        eff_ligand_atoms, eff_donor_idc = self.get_all_effective_ligand_atoms_with_effective_donor_indices('Cu')
        geometry, eff_isomers, eff_isomer_donor_idc, rssd, second_geometry, weight_necessary_for_change = assign_geometry(eff_ligand_atoms, eff_donor_idc)
        # Remove Cu from the isomers to get the real ligand geometry
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
            # immutable_hapdent_idc = [tuple(idc) if isinstance(idc, tuple) else idc for idc in hapdent_donor_idc]
            # immutable_original_hapdent_idc = [tuple(idc) if isinstance(idc, tuple) else idc for idc in self.hapdent_idc]
            assert set(hapdent_donor_idc) == set(self.hapdent_idc), f"Error in conversion of effective donor indices to hapdent donor indices: {hapdent_isomer_idc} vs. {self.hapdent_idc}"
            hapdent_isomer_idc.append(hapdent_donor_idc)

        return geometry, real_isomers, hapdent_isomer_idc, rssd, second_geometry, weight_necessary_for_change

    def get_ligand_output_info(self, max_entries=5) -> dict:
        self._sort_global_props_inplace()
        info = self.global_props.copy()

        truncate_data = {
                            'csd_complex_ids': self.other_ligand_instances['parent_complex_id'],
                            'csd_metal_count': [f'{el}({count})' for el, count in self.count_metals.items()],
                            'csd_metal_os_count': [f'{el}({count})' for el, count in self.mos_counts.items()],
                        }
        for key, data in truncate_data.items():
            n_data = len(data)
            data = data[:max_entries]
            data = ', '.join(data)
            if n_data > max_entries:
                data += f', ... ({n_data - max_entries} more)'
            info[key] = data

        return info

    def write_to_mol_dict(self, include_graph_dict: bool=True) -> dict:
        self._sort_global_props_inplace()
        d = {
                'atomic_props': self.atomic_props,
                'global_props': self.global_props,
        }
        if include_graph_dict:
            d['graph'] = graph_to_dict_with_node_labels(self.graph)

        d['donor_idc'] = self.ligand_to_metal
        d['ligand_instances'] = self.other_ligand_instances
        d['hapdent_idc'] = self.hapdent_idc
        d['geometric_isomers_hapdent_idc'] = self.geometric_isomers_hapdent_idc

        # output = ['warnings', 'atomic_props', 'global_props', 'n_atoms', 'n_hydrogens', 'n_protons', 'graph_hash', 'n_bonds', 'has_bond_order_attribute', 'has_unknown_bond_orders', 'has_good_bond_orders', 'heavy_atoms_graph_hash', 'bond_order_graph_hash', 'stoichiometry', 'original_complex_id', 'local_elements', 'was_connected_to_metal', 'original_metal', 'original_metal_position', 'original_metal_symbol', 'original_metal_os', 'is_centrosymmetric', 'centrosymmetry_ang_dev', 'graph_hash_with_metal', 'heavy_atoms_graph_hash_with_metal', 'has_betaH', 'has_neighboring_coordinating_atoms', 'stats', 'unique_name', 'pred_charge', 'pred_charge_is_confident', 'all_ligand_names', 'identical_ligand_info', 'occurrences', 'same_graph_denticities', 'count_metals', 'n_same_graph_denticities', 'n_metals', 'n_same_graphs', 'has_unconnected_ligands', 'all_ligands_metals', 'same_graph_charges', 'n_pred_charges', 'common_graph_with_diff_n_hydrogens', 'n_electrons', 'odd_n_electron_count', 'has_warnings', 'denticity', 'name', 'ligand_to_metal']
        #
        # # do_not_output_automatically = ['mol', 'node_label', 'rdkit_mol', 'graph', 'coordinates', 'hash', 'csd_code', 'graph_with_metal', 'atomic_index_to_graph_index', 'graph_index_to_atomic_index']
        # for prop in output:
        #     # if not prop in do_not_output_automatically:
        #     d[prop] = getattr(self, prop)

        return d

    @classmethod
    def read_from_mol_dict(cls, dict_: dict):
        """
        Reads the ligand from a provided dictionary.
        """
        # In the old format, these properties were in the dict_ dictionary, instead of in the global_props as in the new format.
        old_props = ['n_hydrogens', 'n_atoms', 'n_protons']
        if all([prop in dict_ for prop in old_props]):
            # The input dictionary is in the old format. Convert it to the new format.
            dict_ = refactor_metalig_entry_from_v1_0_0_to_v1_1_0(dict_)

        necessary_props = {'atomic_props', 'global_props', 'graph', 'donor_idc', 'ligand_instances', 'hapdent_idc', 'geometric_isomers_hapdent_idc'}
        different_props = necessary_props.symmetric_difference(set(dict_.keys()))
        assert not different_props, f"Missing or unexpected properties in the ligand input dictionary: {different_props}"

        # Add default properties for properties not in the json. This is useful for properties which have been introduced in later versions of the code.
        # optional_props = {'warnings': []}
        # for prop, default in optional_props.items():
        #     if not prop in dict_:
        #         dict_[prop] = default

        # other_props = {key: val for key, val in dict_.items() if not key in necessary_props}

        # # Optionally add graph if it is present in the dictionary
        # if 'graph_dict' in dict_ and not (dict_['graph_dict'] is None):
        graph = graph_from_graph_dict(dict_['graph'])
        # else:
        #     graph = None

        return cls(
            atomic_props=dict_["atomic_props"],
            global_props=dict_["global_props"],
            # denticity=dict_["denticity"],
            unique_name=dict_['global_props']['unique_name'],
            ligand_to_metal=dict_['donor_idc'],
            graph=graph,
            ligand_instances=dict_['ligand_instances'],
            hapdent_idc=dict_['hapdent_idc'],
            geometric_isomers_hapdent_idc=dict_['geometric_isomers_hapdent_idc'],
            # warnings=dict_['warnings'],
            # other_props=other_props,
            validity_check=True,
        )

    def to_stk_bb(self):
        """
        this is really only designed for ligands as a normal RCA_Molecule doesn't have the required properties
        """
        return RCA_Mol_to_stkBB(self)

    def get_ase_molecule_with_metal(self, metal: str=None) -> Atoms:
        """
        Get ASE molecule with metal at original metal location. If no metal is specified, the original metal is used.
        """
        return self.get_ase_molecule(add_atoms=[(metal, self.original_metal_position)])

    def get_xtb_descriptors(self):
        """
        Return xtb descriptors for the molecule. Needs xtb installed.
        """
        from dev.src11_machine_learning.utils.utilities_ML import get_xtb_descriptors

        try:
            charge = self.pred_charge
            self.xtb_descriptors = get_xtb_descriptors(self.get_xyz_file_format_string(), charge=charge,
                                                       n_unpaired=None)
        except AttributeError:
            self.xtb_descriptors = {
                                'ionization_potential': np.nan,
                                'electron_affinity': np.nan,
                                'HOMO': np.nan,
                                'LUMO': np.nan,
                                'Dipole_x': np.nan,
                                'Dipole_y': np.nan,
                                'Dipole_z': np.nan
                                }

        return self.xtb_descriptors

    def generate_descriptors(self, get_3D_descriptors: bool = True) -> dict:
        """
        Generates descriptors for the ligand. Needs pymatgen and mordred installed.
        """
        from src11_machine_learning.utils.utilities_ML import get_element_descriptors
        descriptors = {}

        # Global descriptors
        descriptors['denticiy'] = self.denticity
        # descriptors['planar'] = self.planar_check()
        descriptors['molecular_weight'] = self.global_props['molecular_weight']
        descriptors['n_atoms'] = self.n_atoms
        descriptors['n_bonds'] = self.n_bonds

        # Coordinated atoms descriptors
        coords_descriptors = [get_element_descriptors(el) for el in self.local_elements]
        coords_descriptors = {key: np.mean([d[key] for d in coords_descriptors]) for key in
                              coords_descriptors[0].keys()}
        descriptors.update(coords_descriptors)

        # Graph and 3D descriptors
        import mordred
        mol = Chem.MolFromSmiles(self.get_smiles())
        calc = mordred.Calculator(descriptors, ignore_3D=not get_3D_descriptors)
        mordred_descriptors = calc.pandas(mol).to_dict(orient='records')[0]
        descriptors.update(mordred_descriptors)

        return descriptors



class RCA_Complex(RCA_Molecule):

    def __init__(self,
                 mol: Atoms = None,
                 atomic_props: dict = None,
                 global_props: dict = None,
                 pymat_mol=None,
                 graph=None,
                 graph_creating_strategy="default",
                 has_ligands=True,
                 reindex_graph: bool = False,
                 other_props={},
                 **kwargs):

            super().__init__(
                            mol=mol,
                            atomic_props=atomic_props,
                            global_props=global_props,
                            pymat_mol=pymat_mol,
                            graph=graph,
                            graph_creating_strategy=graph_creating_strategy,
                            has_ligands=has_ligands,
                            reindex_graph=reindex_graph,
                            other_props=other_props,
                            **kwargs
                             )
            self.check_if_graph_and_coordinates_are_consistent()
            self.shift_metal_to_origin()

            self.metal = identify_metal_in_ase_mol(self.mol)
            self.metal_atomic_number = DART_Element(self.metal).atomic_number
            self.metal_oxi_state = make_None_to_NaN(self.global_props['metal_oxi_state'])
            self.charge = make_None_to_NaN(self.global_props['charge'])
            self.metal_idx = self.atomic_props['atoms'].index(self.metal)
            self.metal_position = (self.atomic_props['x'][self.metal_idx], self.atomic_props['y'][self.metal_idx], self.atomic_props['z'][self.metal_idx])
            self.fully_connected = nx.is_connected(self.graph)

            self.donor_indices = sorted(self.graph.neighbors(self.metal_idx))
            self.donor_elements = [self.atomic_props['atoms'][idx] for idx in self.donor_indices]
            self.n_donors = len(self.donor_indices)

            self.metal_position = (self.atomic_props['x'][self.metal_idx], self.atomic_props['y'][self.metal_idx], self.atomic_props['z'][self.metal_idx])
            self.fully_connected = nx.is_connected(self.graph)
            self.mol_id = self.global_props['CSD_code']

            # Set parameters for octahedral complexes
            self.is_octahedral = self.check_octahedral()
            mean_distance, sd_distance, max_angular_deviation = self.calculate_distortion_parameters()
            self.global_props['donors_mean_dist'] = mean_distance
            self.global_props['donors_sd_dist'] = sd_distance
            self.global_props['oct_max_ang_dev'] = max_angular_deviation

            self.add_additional_complex_information_to_global_props()

            return

    def has_metal_os(self):
        return ~np.isnan(self.metal_oxi_state)

    def calculate_distortion_parameters(self) -> tuple[float, float, float]:
        """
        Function to calculate distortion parameters of the complex
        :return: sd_distance (standard deviation of neighbour distances from the metal atom) and max_angular_deviation (maximum deviation from 90 or 180 degrees in bond angles)
        """
        # Calculate neighbour distances and positions
        metal_pos = np.array(self.metal_position)
        neighbours = self.donor_indices
        neighbour_positions = []
        neighbour_distances = []
        for neighbour_idx in neighbours:
            neighbour_pos = np.array([self.atomic_props['x'][neighbour_idx], self.atomic_props['y'][neighbour_idx], self.atomic_props['z'][neighbour_idx]])
            neighbour_positions.append(neighbour_pos)
            neighbour_distances.append(np.linalg.norm(neighbour_pos - metal_pos))

        # Calculate standard deviation of neighbour distances
        if len(neighbour_distances) > 0:
            mean_distance = float(np.mean(neighbour_distances))
            sd_distance = float(np.std(neighbour_distances))
        else:
            mean_distance = np.nan
            sd_distance = np.nan

        # Calculate maximum angular deviation for octahedral complexes
        if len(neighbours) == 6:
            max_angular_deviation = 0
            for i in range(len(neighbour_positions)):
                for j in range(i+1, len(neighbour_positions)):
                    # Calculate the vector from the metal to the neighbours
                    vector_i = neighbour_positions[i] - metal_pos
                    vector_j = neighbour_positions[j] - metal_pos

                    # Calculate the angle between these vectors
                    cos_angle = np.dot(vector_i, vector_j) / (np.linalg.norm(vector_i) * np.linalg.norm(vector_j))
                    cos_angle = np.clip(cos_angle, -1, 1)  # Ensure cos_angle is within the valid range [-1, 1]
                    angle = np.arccos(cos_angle) * 180 / np.pi  # Convert to degrees

                    # Check the angular deviation from 90 or 180 degrees
                    angular_deviation = min(abs(angle - 90), abs(angle - 180))
                    max_angular_deviation = max(max_angular_deviation, angular_deviation)
        else:
            max_angular_deviation = np.nan

        return mean_distance, sd_distance, max_angular_deviation

    def check_octahedral(self, angular_threshold=40) -> bool:
        """
        Function to check if the complex is octahedral or not
        The conditions for a complex to be octahedral are:
            - The complex must have 6 donor atoms
            - The maximum deviation from 90 or 180 degrees in bond angles should be less than `angular_threshold`
        :return: True if the complex is octahedral, False otherwise
        """
        if len(self.donor_indices) != 6:
            return False

        _, _, max_angular_deviation = self.calculate_distortion_parameters()
        is_oct = bool(max_angular_deviation > angular_threshold)    # bool() for json serialisation

        return is_oct

    def count_ligands_with_stoichiometry(self, atoms: list, only_connected=False):
        if isinstance(atoms, str):
            atoms = [atoms]
        atoms = sorted(atoms)

        n = 0
        for lig in self.ligands:
            if not only_connected or lig.was_connected_to_metal:
                lig_atoms = sorted(lig.atomic_props['atoms'])
                if atoms == lig_atoms:
                    n += 1

        return n

    def count_n_unconnected_ligands(self, max_n_atoms: np.inf) -> int:
        """
        Counts the number of unconnected ligands with a maximum number of atoms of `max_n_atoms`.
        @param max_n_atoms: Maximum number of atoms of ligands to count
        @return: Number of ligands with this maximum number of atoms
        """
        n = sum([not lig.was_connected_to_metal and len(lig.atomic_props['atoms']) <= max_n_atoms for lig in self.ligands])
        return n


    def count_atoms_in_ligands(self, atoms: Union[list, str], only_if_connected_to_metal: bool=False, per_ligand:bool=False) -> Union[list, int]:
        """
        Returns the number of occurrences of the specified elements in the complex.
        @param atoms: list of atoms to count the total occurrences
        @param only_if_connected_to_metal: If True, ignore unconnected ligands
        @param per_ligand:If True, returns a list of occupancies for each ligand. Otherwise returns the total for all complexes.
        @return: Number of occurrences of specified elements, either as list per ligand or the total as integer
        """
        if isinstance(atoms, str):
            atoms = [atoms]

        n = []
        for lig in self.ligands:
            if only_if_connected_to_metal and lig.was_connected_to_metal:
                n.append(sum([el in atoms for el in lig.atomic_props['atoms']]))

        if not per_ligand:
            n = sum(n)

        return n

    def get_only_complex_graph_connected_to_metal(self, atom_label: object = 'node_label') -> nx.Graph:
        """
        Returns the graph of only the metal complex without unconnected ligands.
        """
        complex_graph = get_only_complex_graph_connected_to_metal(graph=self.graph, atom_label=atom_label, metal=self.metal)
        return complex_graph

    def has_fragment(self, frag: Union[str, list]) -> bool:
        """
        Checks whether the complex graph has an unconnected fragment with the elements specified in frag.
        @param frag: Either a string specifying a single atom or a list of strings specifying several atoms. E.g. frag='Cl' checks if any unconnected chloride exists. frag='[H, H, O]' checks whether any unconnected water exists. The order in the list doesn't matter but the number of occurrences of elements does matter.
        """
        if isinstance(frag, str):
            frag = [frag]

        fragments = [sorted(comp) for comp in nx.connected_components(self.graph)]
        el_fragments = [[self.atomic_props['atoms'][i] for i in fragment] for fragment in fragments]

        has_unconnected_fragment = False
        frag = sorted(frag)
        for fragment in el_fragments:
            fragment = sorted(fragment)
            if fragment == frag:
                has_unconnected_fragment = True
                break

        return has_unconnected_fragment

    def count_ligands_containing_only(self, atoms: Union[str, list], denticity_range: list=[], n_atoms_range: list=[], except_elements=[]) -> int:
        """
        Count how many ligands contain only atoms specified in `atoms`, except for elements in `except elements`.
        @param atoms: atoms to check for.
        @param denticity_range: if specified, count ligands only if it has a denticity in this range (inclusive)
        @param n_atoms_range: if specified, count ligands only if the number of atoms which are not excluded by `except_element` is in this range (inclusive)
        @param except_elements: Ignore these elements in the ligands when checking the presence of the atoms in the ligand and when checking n_atoms_range.
        @return: Integer
        """
        n = 0
        for lig in self.ligands:
            correct_denticity = is_between(lig.denticity, denticity_range)
            correct_n_atoms = is_between(len([a for a in lig.atomic_props['atoms'] if not a in except_elements]), n_atoms_range)
            if correct_denticity and correct_n_atoms and lig.contains_only(atoms, except_elements=except_elements):
                n += 1

        return n


    def count_coordinating_atoms_with_distance_to_metal_greater_than(self, distance, element=None, max_n_atoms=np.inf) -> int:
        """
        Returns the number of coordinating atoms of type `element` with a distance to the metal greater than `distance`.
        @param max_n_atoms: Include only ligands with a number of atoms up to this number.
        """
        n = 0
        for lig in self.ligands:
            for idx, el in enumerate(lig.local_elements):
                correct_atom = el == element or element is None
                if correct_atom and len(lig.atomic_props['atoms']) <= max_n_atoms:
                    if lig.get_atomic_distance_to_original_metal(mode='coordinating')[idx] > distance:
                        n += 1

        return n

    def has_consistent_stoichiometry_with_CSD(self, print_different_elements=False) -> bool:
        try:
            csd_stoichiometry = self.global_props['CSD_stoichiometry']
            pattern = r'[A-Z][a-z]?\d*'
            csd_stoichiometry = ' '.join(re.findall(pattern, csd_stoichiometry))
            csd_comp = Composition(csd_stoichiometry)
        except KeyError:
            warnings.warn('Global property `CSD_stoichiometry` not found. Skip check for consistent stoichiometry with xyz.')
            return True
        except:
            warnings.warn(f'CSD_stoichiometry could not be parsed for complex {self.mol_id}. Skip check for consistent stoichiometry with xyz.')
            return True

        xyz_comp = Composition(self.stoichiometry)
        consistent = csd_comp.almost_equals(xyz_comp)

        if not consistent and print_different_elements:
            xyz_elements = [DART_Element(el).symbol for el in xyz_comp.elements]
            csd_elements = [DART_Element(el).symbol for el in csd_comp.elements]
            different_elements = set(xyz_elements).symmetric_difference(csd_elements)
            if different_elements:
                print(f'Differing elements in stoichiometries of xyz and CSD in complex {self.mol_id}: {list(different_elements)}')
            diff_dict = {el: int(xyz_comp.elements[el]-csd_comp.elements[el]) for el in set(csd_elements+xyz_elements)}
            diff_dict = {key: val for key, val in diff_dict.items() if val != 0}
            print(f'Different elements counts for complex {self.mol_id}: {diff_dict}')

        return consistent

    def has_consistent_stoichiometry_with_smiles(self, smiles: str, ignore_element_count: bool=False, print_warnings: bool=True, only_complex: bool=False) -> bool:
        try:
            mol_graph = read_smiles(smiles, explicit_hydrogen=True, zero_order_bonds=False)

            if only_complex:
                _, frag_atoms = get_graph_fragments(graph=mol_graph, atom_label='element')
                atoms = [a for a in frag_atoms if self.metal in a][0]
            else:
                atoms, _ = get_sorted_atoms_and_indices_from_graph(mol_graph, atom_label='element')

            comp = Composition(''.join(atoms))
        except KeyError:
            if print_warnings:
                warnings.warn('Global property `smiles` not found. Skip check for consistent stoichiometry with xyz.')
            return True
        except:
            if print_warnings:
                warnings.warn(f'Smiles could not be parsed for complex {self.mol_id}. Skip check for consistent stoichiometry with xyz.')
            return True

        if only_complex:
            _, frag_atoms = self.get_graph_fragments()
            atoms = [a for a in frag_atoms if self.metal in a][0]
            xyz_comp = Composition(''.join(atoms))
        else:
            xyz_comp = Composition(self.stoichiometry)

        if ignore_element_count:
            # Check only if the same elements occur, not if the elements have the same count.
            xyz_elements = [DART_Element(el).symbol for el in xyz_comp.elements]
            csd_elements = [DART_Element(el).symbol for el in comp.elements]
            different_elements = set(xyz_elements).symmetric_difference(csd_elements)
            consistent = len(different_elements) == 0
        else:
            consistent = comp.almost_equals(xyz_comp)

        return consistent


    def complex_is_biggest_fragment(self, allow_complexes_greater_than: int = np.nan) -> bool:
        """
        Checks whether the complex (fragment with the transition metal) is the fragment with the most atoms. This is useful to check whether the transition metal might not belong to an actual complex but to a counterion.
        :param allow_complexes_greater_than: Always return true for complexes with a higher number of atoms than this parameter
        """
        _, el_fragments = self.get_graph_fragments()

        max_n_other_atoms = 0
        n_complex_atoms = 0
        for atoms in el_fragments:
            is_complex = any([DART_Element(atom).atomic_number in metals_in_pse for atom in atoms])
            if is_complex:
                assert n_complex_atoms == 0, f'There seem to be complexes with more than one transition metal in complex {self.mol_id}.'
                n_complex_atoms = len(atoms)
            else:
                if len(atoms) > max_n_other_atoms:
                    max_n_other_atoms = len(atoms)
        assert n_complex_atoms != 0, f'There seem to be no transition metals in complex {self.mol_id}.'

        if max_n_other_atoms > n_complex_atoms and n_complex_atoms <= allow_complexes_greater_than:
            is_biggest_fragment = False
            # print(f'N complex atoms in {self.mol_id}: {n_complex_atoms}')
        else:
            is_biggest_fragment = True

        return is_biggest_fragment

    def add_additional_complex_information_to_global_props(self):
        info = {}
        if 'partial_charge' in self.atomic_props:
            metal_q = self.atomic_props['partial_charge'][self.atomic_props['atoms'].index(self.metal)]
            info['metal_partial_charge'] = metal_q

        update_dict_with_warning_inplace(
                                            dict_to_update=self.global_props,
                                            dict_with_information=info
                                        )

        return info

    def shift_metal_to_origin(self):
        """
        Actually, this method is outdated not required.
        However, the idea was to shift the molecule to the origin, so that the metal is right in the origin
        which is the first part - the shift

        and afterwards, rotate the the metal, so that one of the functional atoms is algined with the x-axis
        """
        #
        # get shift vector
        metal_symb = identify_metal_in_ase_mol(self.mol)
        metal_index = self.atomic_props["atoms"].index(metal_symb)

        # do the shifting
        for key in ["x", "y", "z"]:
            self.atomic_props[key] = [value - self.atomic_props[key][metal_index] for value in self.atomic_props[key]]

        """
        # Now we rotate. First we identify the element we want to rotate on
        metal_node = find_node_in_graph_by_label(G=self.graph, label_to_find=metal_symb, expected_hits=1)
        rotation_element_index = list(self.graph.neighbors(metal_node))[0]
        # Next, we obtain the rotation matrix
        rotation_element_vector = np.array([self.atomic_props[key][rotation_element_index] for key in ["x", "y", "z"]])
            # now we can idenftify the desired vector, on which we'd like to rotate the vector
        desired_rotation = np.array([np.linalg.norm(rotation_element_vector), 0, 0])
        # and we can thus find the rotation matrix
        A = R.align_vectors(rotation_element_vector.reshape(-1, 1).T, desired_rotation.reshape(-1, 1).T)[0].as_matrix()
        for index, _ in enumerate(self.atomic_props["x"]):
            location_vec_for_index = np.array([self.atomic_props[key][index] for key in ["x", "y", "z"]])
            # now we rotate
            v_new = location_vec_for_index.dot(A)
            # and we modify the atomic properties
            for i, key in enumerate(["x", "y", "z"]):
                self.atomic_props[key][index] = v_new.tolist()[i] if abs(v_new.tolist()[i]) > 0.001 else 0
        """
        return

    def get_output_info(self) -> dict:

        info = {
            'Complex ID': self.mol_id,
            'Stoichiometry': self.stoichiometry,
            'Metal': self.metal,
            'Metal OS': self.metal_oxi_state,
            'Charge': self.charge,
            'Number of Atoms': self.n_atoms,
            'Fully Connected': self.fully_connected,
            'Coordinating Atoms': self.donor_elements,
            'Number of Donors': self.n_donors,
            'Octahedral': self.is_octahedral,
            # 'Planarity': self.calculate_planarity(),
            'Haptic': any(lig.has_neighboring_coordinating_atoms for lig in self.ligands),
            'Beta-Hydrogen': any(lig.has_betaH for lig in self.ligands),
            'Ligand Unique Names': [lig.unique_name for lig in self.ligands],
            **self.global_props,
            }

        return info

    def de_assemble(self, inherit_global_properties: list = ['CSD_code']):
        """
        now only graph based, makeslife waaay easier

        All that is based on the following assumption:
        the nodes are denoted as integers and we shall assume that these integers
        correspond to their index in the atomic properties, i.e. the position
        so, for example for node 1, the x-coordinate can be found by
        self.atomic_props["x"][1] or more general self.atomic_props["x"][node]
        Ich glaube das haelt auf jeden Fall fuer alle selber erzeugten Graphen
        den Dummen Testparamater muss ich mitschleppen, um das gegen meine alte, ultra behinderte Grapherstellung testen
        zu koennen
        """
        if not hasattr(self, "ligands"):
            self.ligands = []

        inherit_global_properties = self.check_input_inherit_global_properties(inherit_global_properties)

        atoms, idc = get_sorted_atoms_and_indices_from_graph(self.graph)
        if 'atoms' in self.atomic_props:
            # if not atoms == self.atomic_props['atoms']:
            #     breakpoint()

            assert atoms == self.atomic_props['atoms'], 'Order of atoms in graph and in atomic_props doesn\'t match.'

        # first we gather some information around the metal in the initial graph
        graph = deepcopy(self.graph)
        metal_in_complex = identify_metal_in_ase_mol(self.mol)
        metal_node = find_node_in_graph_by_label(G=graph, label_to_find=metal_in_complex, expected_hits=1)
        metal_neighbors = list(graph.neighbors(metal_node))  # all neighbor nodes of the metal
        metal_neighbor_elements = [graph.nodes[i]['node_label'] for i in metal_neighbors]

        for i, el in zip(idc, atoms):
            #graph.nodes[i]['orig_idx'] = i
            #graph.nodes[i]['metal_neighbor'] = i in metal_neighbors
            assert graph.nodes[i]['node_label'] == el, 'atom and index don\'t match.'

        # next we create the ripped graph
        ripped_graph = deepcopy(graph)
        ripped_graph.remove_node(metal_node)

        conn_components = [sorted(comp) for comp in
                           nx.connected_components(ripped_graph)]  # sorting of components very important
        conn_components = [comp for comp in
                           sorted(conn_components,
                                  key=str)]  # important: sort by string of list, makes order of ligands unique

        for component in conn_components:
            # if this set is empty, the ligand has no connection to the metal
            functional_atom_indices = sorted(list(set(component).intersection(set(metal_neighbors))))
            assert max(component) <= len(
                ripped_graph), 'Indices dont make sense. Most likely this is an implementation error where due to the deletion of the metal atom the indices of original and ripped graph dont match.'

            denticity = len(functional_atom_indices)
            if denticity == 0:
                # because in the assembly 0 is the placeholder for the reactant, whereas -1 means this is just a remainder in the .xyz, not
                # connected to the metal at all
                denticity = -1

            # with that it is insanely easy to determine the atomic properties of the ligand
            assert component == sorted(
                component), 'The list of ligand indices is not sorted, but that is assumed in many parts of this project.'
            ligand_atomic_props = {key_: [el for i, el in enumerate(item_) if i in component] for key_, item_ in
                                   self.atomic_props.items()}
            ligand_atomic_props['original_complex_indices'] = component

            ligand_graph = deepcopy(ripped_graph.subgraph(component))
            atoms_lig, idc_lig = get_sorted_atoms_and_indices_from_graph(ligand_graph)

            # problem: the functional_atom_indices are the indices in the full original metal, rather than the ligand only
            # so we have to convert them to the index in the ligand_atomic_props
            local_indices = [component.index(ind) for ind in functional_atom_indices]
            local_elements = [ligand_atomic_props['atoms'][i] for i in local_indices]

            # Doublechecking
            if local_indices != []:  # otherwise error for unconnected ligands
                assert max(local_indices) < len(
                    ligand_graph), 'local_indices make no sense, an index is greater than the number of elements.'
            assert all(
                el in metal_neighbor_elements for el in local_elements), 'Inconsistent elements of the metal neighbors.'
            assert local_indices == sorted(local_indices), 'local_indices is not sorted but should be.'
            assert atoms_lig == ligand_atomic_props['atoms'], 'elements of graph and atomic_props not consistent'
            assert [atoms[i] for i in component] == ligand_atomic_props[
                'atoms'], 'ligand atoms not consistent with original complex atoms.'

            ligand_name, csd = self.ligand_naming(denticity, self.ligands)

            ligand_global_props = {prop: self.global_props[prop] for prop in inherit_global_properties}
            new_lig = RCA_Ligand(denticity=denticity,
                                 ligand_to_metal=local_indices,
                                 atomic_props=ligand_atomic_props,
                                 name=ligand_name,
                                 graph=ligand_graph,
                                 global_props=ligand_global_props,
                                 original_metal=DART_Element(metal_in_complex).atomic_number,
                                 original_metal_position=self.metal_position,
                                 original_metal_os=self.metal_oxi_state
                                 )

            self.ligands.append(new_lig)

        if not self.ligands:
            print(f'WARNING: Complex {self.global_props["CSD_code"]} has no ligands extracted.')

    def get_smiles(self, only_core_complex: bool=False) -> str:
        full_smiles = super().get_smiles()
        smiles = full_smiles
        if only_core_complex:
            smiles = [sm for sm in smiles.split('.') if self.metal in sm]
            assert len(smiles) == 1, 'There should be exactly one SMILES string containing the metal.'
            smiles = smiles[0]

        # Assert
        smiles_graph = pysmiles.read_smiles(smiles, explicit_hydrogen=True)
        smiles_hash = get_graph_hash(smiles_graph, node_attr='element')
        core_graph = self.get_only_complex_graph_connected_to_metal()
        core_graph_hash = get_graph_hash(core_graph)
        if not smiles_hash == core_graph_hash and len(core_graph.nodes) < 9999999999:
            view_graph(core_graph, save_path='/Users/timosommer/Downloads/core_graph.png')
            view_graph(smiles_graph, save_path='/Users/timosommer/Downloads/smiles_graph.png', node_label='element')

            a = 1

        # assert smiles_hash == core_graph_hash, 'SMILES and graph hash do not match.'

        return smiles

    def generate_descriptors_of_complex_graph(self, only_core_complex: bool=True, xtb=False):
        """
        Generates descriptors for the complex and its ligands.
        """
        from dev.src11_machine_learning.utils.utilities_ML import get_element_descriptors
        from dev.src11_machine_learning.dataset_preparation.descriptors import RAC
        descriptors = {}

        # Complex descriptors
        descriptors['n_ligands'] = sum(1 for lig in self.ligands if lig.denticity > 0 or not only_core_complex)

        # Metal center descriptors
        metal_descriptors = get_element_descriptors(self.metal)
        metal_descriptors = {f'metal_{key}': val for key, val in metal_descriptors.items()}
        descriptors.update(metal_descriptors)

        # Coordinating_atom descriptors
        coords_descriptors = [get_element_descriptors(el) for el in self.donor_elements]
        coords_descriptors = {key: np.mean([d[key] for d in coords_descriptors]) for key in coords_descriptors[0].keys()}
        coords_descriptors = {f'coord_{key}': val for key, val in coords_descriptors.items()}
        descriptors.update(coords_descriptors)

        # Graph descriptors
        # RDKit descriptors
        # smiles = self.get_smiles(only_core_complex=only_core_complex)
        # rdkit_descriptors = RDKit_2D([smiles]).compute_2Drdkit().to_dict(orient='records')[0]
        # descriptors.update(rdkit_descriptors)
        # Own RAC descriptors
        graph = self.get_only_complex_graph_connected_to_metal() if only_core_complex else self.graph
        graph_descriptors, labels = RAC().molecule_autocorrelation(mol=graph, return_labels=True)
        graph_descriptors = {f'own_graph_{key}': val for key, val in zip(labels, graph_descriptors)}
        descriptors.update(graph_descriptors)


        # XTB descriptors of ligands
        if xtb:
            xtb_descriptors = []
            for lig in self.ligands:
                if lig.denticity > 0 or not only_core_complex:
                    xtb_descriptors.append(lig.get_xtb_descriptors())


        return descriptors

    def write_to_mol_dict(self, include_graph_dict: bool=True, include_ligands: bool=True) -> dict:
        d = {}

        # Manually initialize special fields
        if include_graph_dict:
            d['graph_dict'] = graph_to_dict_with_node_labels(self.graph)
        d['ligands'] = [lig.write_to_mol_dict() for lig in self.ligands]

        # Do not output write these fields to the output dictionary, mostly because they are not json serializable
        do_not_output_automatically = ['ligands', 'node_label', 'mol', 'smiles', 'rdkit_mol', 'graph', 'coordinates', 'hash', 'csd_code', 'graph_with_metal', 'atomic_index_to_graph_index', 'graph_index_to_atomic_index']

        for prop, val in vars(self).items():
            if not prop in do_not_output_automatically:
                d[prop] = val

        return d

    @classmethod
    def read_from_mol_dict(cls, dict_: dict, **kwargs):
        """
        Reads the ligand from a provided dictionary.
        """
        necessary_props = ["atomic_props", "global_props", "graph_dict"]
        assert set(necessary_props).issubset(set(dict_.keys())), f'Any of the necessary keys {necessary_props} is not present.'

        if 'ligands' in dict_:
            dict_['ligands'] = [RCA_Ligand.read_from_mol_dict(lig) for lig in dict_['ligands']]
            # This sounds stupid but this is because otherwise the RCA_molecule class sets up self.ligands = [] which then collides when the actual ligands which are read in here are added because for safety the code checks that it doesn't overwrite anything.
            has_ligands = False
        else:
            has_ligands = True

        if 'total_q' in dict_:
            dict_['charge'] = dict_['total_q']
            del dict_['total_q']

        other_props = {key: val for key, val in dict_.items() if not key in necessary_props}

        # Optionally add graph if it is present in the dictionary
        if 'graph_dict' in dict_ and not (dict_['graph_dict'] is None):
            graph = graph_from_graph_dict(dict_['graph_dict'])
        else:
            graph = None

        return cls(
            atomic_props=dict_["atomic_props"],
            global_props=dict_["global_props"],
            graph=graph,
            has_ligands=has_ligands,
            other_props=other_props,
            **kwargs
        )






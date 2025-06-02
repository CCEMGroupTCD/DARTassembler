#########################################################################################
# This file contains the classes and methods that are used to process the input data    #
# and generate the assembled transition metal complexes                                 #
#########################################################################################
from numpy import ndarray

from DARTassembler.src.assembly.ligand_geometries import try_all_geometrical_isomer_possibilities
from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.DataBase import RCA_Ligand, LigandDB
from scipy.optimize import differential_evolution, linear_sum_assignment, dual_annealing, brute
from DARTassembler.src.constants import Periodic_Table as PerTab
from typing import Dict, Any, List, Optional, Tuple, Union
import warnings
from scipy.spatial.transform import Rotation as R
from ase.visualize import view
import plotly.express as px
from scipy.spatial.distance import cdist
from sklearn.cluster import MeanShift, estimate_bandwidth, dbscan
import pandas as pd
from ase import Atoms
import numpy as np
import itertools
import plotly.graph_objects as go
import ase
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt


from DARTassembler.src.ligand_extraction.io_custom import load_json
from DARTassembler.src.ligand_extraction.utilities_graph import view_graph, graph_to_dict_with_node_labels, \
    graph_from_graph_dict
from DARTassembler.src.ligand_extraction.Molecule import RCA_Molecule
import networkx as nx
from copy import deepcopy


def get_merged_graph_from_ligands_and_metal_centers(
                                                    ligands: list[RCA_Ligand],
                                                    metal_centers: list[list[ase.Atom]]
                                                    ) -> tuple[nx.Graph, list, list]:
    """
    Merges the graphs from the ligands into one graph. The metal is added as a node with index 0 and connected to the donor atoms of the ligands.
    # todo
        This function does not yet work perfectly for multidentate bridging ligands. The donor atoms of the ligands are not correctly connected to the metal center. E.g. for a bidentate bridging atom with two metal donors, all donor atoms are connected to each metal of the two metal centers. This needs to be either fixed or documented.
    :param ligands: List of RCA_Ligand objects
    :param metal_centers: List of list of connected metal centers for each ligand
    :return: Tuple of the merged graph of the complex, the indices of the ligand atoms and the indices of the ligand donor atoms
    """
    ligand_graphs = [deepcopy(lig.graph) for lig in ligands]
    unique_metal_centers = get_all_unique_metal_centers(metal_centers)

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
        graph.add_nodes_from(ligand_graph.nodes(data=True))     # add ligand nodes
        graph.add_edges_from(ligand_graph.edges())              # add ligand edges

    # Connect the metal centers to the ligands
    ligand_donor_indices = [[] for _ in ligands]
    for i, (ligand, ligand_metal_centers, ligand_graph) in enumerate(zip(ligands, metal_centers, ligand_graphs)):
        for metal_center in ligand_metal_centers:
            unique_metal_center_idx = [i for i, atom in enumerate(unique_metal_centers) if are_atoms_equal(atom, metal_center)][0]
            for atomic_donor_idx in ligand.ligand_to_metal:
                assert ligand.atomic_props['atoms'][atomic_donor_idx] in ligand.local_elements, f"Atom {ligand.atomic_props['atoms'][atomic_donor_idx]} is not a donor atom of ligand."
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
    for ligand in ligands:
        all_atomic_elements += ligand.atomic_props['atoms']
    all_graph_elements = [graph.nodes[node]['node_label'] for node in sorted(graph.nodes)]
    assert all_graph_elements == all_atomic_elements, f"The graph elements do not match the atomic elements: {all_graph_elements} vs {all_atomic_elements}!"

    atomic_donor_elements = sorted([el for lig in ligands for el in lig.local_elements])
    graph_donor_elements = sorted([graph.nodes[node]['node_label'] for idc in ligand_donor_indices for node in sorted(graph.nodes) if node in idc])
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

def get_rotated_ligands(ligands: List[RCA_Ligand], target_vectors: List[List[float]],
                        ligand_origins: List[List[float]]) -> list[list[Atoms]]:
    rotated_ligands = []
    for ligand, target_vectors, origin in zip(ligands, target_vectors, ligand_origins):
        # Extract the geometry and donor atoms of the ligand
        atoms, donor_atoms = ligand.get_isomers_effective_ligand_atoms_with_effective_donor_indices()
        # Cast the target vectors to numpy arrays
        target_vectors = [np.array(v) for v in target_vectors]
        # Align the donor atoms of the ligand to the target vectors
        ligand_isomers, donor_atoms_ordered, rssd = try_all_geometrical_isomer_possibilities(atoms=atoms,
                                                                                             donor_idc=donor_atoms[0],
                                                                                             target_vectors=target_vectors)
        # Remove the dummy atom from the haptic ligands
        ligand_isomers = remove_haptic_dummy_atom(atoms_list=ligand_isomers, dummy_atom="Cu",
                                                  donor_atoms_idc=ligand.hapdent_idc)

        # Append the rotated ligands to the list
        rotated_ligands.append(ligand_isomers)

    return rotated_ligands

def remove_haptic_dummy_atom(atoms_list: List[Atoms], dummy_atom: str, donor_atoms_idc: Tuple[Tuple[int]]):
    """
    Removes the dummy atom from the generated isomers.
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
            dummy_idc = [i for i, atom in enumerate(atoms) if atom.symbol == dummy_atom]
            dummy_idc.sort(
                reverse=True)  # This is important so that the larger index is removed first so as not to change the index of the other atoms
            for dummy_idx in dummy_idc:
                atoms.pop(dummy_idx)
        return atoms_list

def get_all_unique_metal_centers(metal_centers) -> List[ase.Atom]:
    """
    Get a list of all unique metal centers.
    :return: List of ase.Atom objects
    """
    unique_metal_centers = [metal_centers[0][0]]    # initialize the list with the first metal center
    for metal_list in metal_centers:
        for metal in metal_list:
            metal_in_list = any([are_atoms_equal(metal, m) for m in unique_metal_centers])
            if not metal_in_list:
                unique_metal_centers.append(metal)

    return unique_metal_centers

def are_atoms_equal(atom1: ase.Atom, atom2: ase.Atom) -> bool:
    """
    Compares two atoms to see if they are equal.
    """
    return atom1.symbol == atom2.symbol and np.allclose(atom1.position, atom2.position)

class LigandSpec:
    """
    Represents a ligand entry in the batch file's geometry section
    It can handle different numbers of vectors dynamically and unique ligand specific input options
    LigandSpec --> Ligand Specification
    """

    def __init__(self, name: str, data: Dict[str, Any]):
        """
        Initializes a LigandSpec object
        :param name: The ligand's identifier (e.g., 'ligand_1')
        :param data: Dictionary containing ligand properties
        :raises ValueError: If required keys are missing
        """
        self.name = name
        self.data = data
        self.origin = self._get_vector(key="origin", required=True)
        self.vectors = self._extract_vectors()
        self.ligand_db_path = data.get("ligand_db", None)
        self.temp_dent = data.get("temp_dent", None)  # TODO: WARNING this is a temporary fix and will need to be removed in the future
        self.ligand_db = LigandDB.load_from_json(path=self.ligand_db_path, n_max=30)
        self.update_geometry = data.get("update_geometry", None)

        self.effective_ligand_coordination_number = len(self.vectors)

    def _extract_vectors(self) -> Dict[str, List[float]]:
        """
        Extracts all vector_X keys dynamically
        :return: Dictionary of vector names to their respective coordinates
        """
        vectors = {key: self._get_vector(key=key) for key in self.data if key.startswith("vector")}
        return vectors

    def _get_vector(self, key: str, required: bool = False) -> Optional[List[float]]:
        """
        Retrieves a vector from the ligand data
        :param key: Key to retrieve
        :param required: If True, raises an error if the key is missing
        :return: List of float coordinates if found, otherwise None
        :raises ValueError: If the required key is missing or malformed
        """
        key_count = sum(1 for k in self.data if k == key)  # Count occurrences of key

        if key_count > 1:
            raise ValueError(f"Fatal Error: Ligand '{self.name}' must have exactly one '{key}' key, found {key_count}.")

        value = self.data.get(key)
        if required and value is None:
            raise ValueError(f"Fatal Error: '{key}' is required for ligand '{self.name}'.")

        if value is not None:
            if not isinstance(value, list) or len(value) != 3:
                raise ValueError(f"Fatal Error: '{key}' in ligand '{self.name}' must be a 3-element list.")
            return [float(v) for v in value]

        return None


class MetalSpec:
    def __init__(self, name: str, data: Dict[str, Any]):
        """
        Initializes a MetalSpec object. This object is used to store input instructions concerning a metal atom
        :param name: The metal's identifier (e.g., 'metal_1')
        :param data: Dictionary containing metal properties
        :raises ValueError: If required keys are missing
        """
        self.name = name
        self.data = data
        self.metal_type = self._get_metal_type()
        self.metal_oxidation_state = self._get_metal_oxidation_state()
        self.coord = self._get_origin()

    def _get_origin(self) -> List[float]:
        """
        Retrieves the metal's coordinate entry
        :return: The metal coordinates as a list of floats
        :raises ValueError: If the coordinate format is incorrect
        """
        key_word = "origin"
        if key_word not in self.data or not isinstance(self.data[key_word], list) or len(self.data[key_word]) != 3:
            raise ValueError(f"Fatal Error: {key_word} must be a 3-element list for metal '{self.name}'")
        return [float(x) for x in self.data[key_word]]

    def _get_metal_type(self) -> str:
        """
        Retrieves the metal type from the data
        :return: str i.e. Ru, Mn, Fe, etc
        """
        # check that metal exists on the periodic table
        metal_str = str(self.data.get("metal_type", ""))
        if metal_str not in PerTab.all_atomic_symbols:
            raise ValueError(f"Fatal Error: Metal '{metal_str}' not found in the periodic table.")
        return metal_str

    def _get_metal_oxidation_state(self) -> int:
        """
        Retrieves the metal oxidation state from the data
        :return: The metals oxidation state i.e. +1, +2, +3, etc
        """
        try:
            oxidation_state = int(self.data["metal_oxidation_state"])
        except (KeyError, ValueError, TypeError):
            raise ValueError(f"Fatal Error: Invalid oxidation state for metal '{self.name}'.")
        return oxidation_state


class BatchInput:
    # Define the required and optional keys for the batch input
    REQUIRED_KEYS = {"name", "random_seed", "max_num_complexes", "total_charge", "geometry", "isomer_instruction"}
    OPTIONAL_KEYS = {"description", "metadata"}

    def __init__(self, batch: Dict[str, Any]):
        """
        Parses and validates a single batch entry from the assembly YAML input.
        Initializes ligand and metal specifications and computes total oxidation state.
        """

        self.batch = batch
        self._validate_batch()

        # Extract and validate the required keys
        self.batch_name = self._get_and_cast(key="name", expected_type=str)
        self.random_seed = self._get_and_cast(key="random_seed", expected_type=float)
        self.max_num_complexes = self._get_and_cast(key="max_num_complexes", expected_type=int)
        self.total_charge = self._get_and_cast(key="total_charge", expected_type=int)
        self.geometry = self._get_and_cast(key="geometry", expected_type=list)  # Expecting a list of dictionaries

        # Metals and ligands
        self.ligands = []
        self.metals = []
        self._process_geometry()

        # Calculate the total metal oxidation state by summing the oxidation states of all metals
        self.total_metal_oxidation_state = int(sum([metal.metal_oxidation_state for metal in self.metals]))

    def _validate_batch(self) -> None:
        """
        Validates the provided batch dictionary, ensuring all required keys exist
        :raises ValueError: If a required key is missing
        """
        # Ensure the batch is a dictionary
        if not isinstance(self.batch, dict):
            raise ValueError("Fatal Error: Batch input must be a dictionary.")

        # Check for missing and extra keys
        missing_keys = self.REQUIRED_KEYS - self.batch.keys()
        extra_keys = self.batch.keys() - self.REQUIRED_KEYS - self.OPTIONAL_KEYS

        if missing_keys:
            raise ValueError(f"Fatal Error: Missing required keys in input file: {', '.join(missing_keys)}")

        if extra_keys:
            print(f"Warning: Unrecognized keys found in batch: [{', '.join(extra_keys)}] These will be ignored")

    def _get_and_cast(self, key: str, expected_type: type) -> Any:
        """
        Safely retrieve a value from the batch dictionary and convert it to the expected type
        :param key: The key to retrieve from the batch dictionary
        :param expected_type: The type to which the value should be converted
        :return: The value converted to the expected type
        """
        value = self.batch.get(key)
        if not isinstance(value, expected_type):
            try:
                value = expected_type(value)
            except Exception:
                raise ValueError(f"Fatal Error: '{key}' must be of type {expected_type.__name__}")
        return value

    @staticmethod
    def _validate_metal(metal_data: Dict[str, Any]) -> List[float]:
        """
        Validate the metal's coordinate entry
        :param metal_data: Dictionary containing metal position data
        :return: The metal coordinates as a list of floats
        :raises ValueError: If the coordinate format is incorrect
        """
        if "coord" not in metal_data or not isinstance(metal_data["coord"], list) or len(metal_data["coord"]) != 3:
            raise ValueError("Fatal Error: 'metal' must contain a 'coord' key with a 3-element list")
        return [float(x) for x in metal_data["coord"]]

    def _process_geometry(self) -> None:
        """
        Processes the 'geometry' key to extract metal and ligand entries
        """

        for entry in self.geometry:
            if not isinstance(entry, dict) or len(entry) != 1:
                raise ValueError("Fatal Error: 'geometry' must be a list of dictionaries")

            key, value = next(iter(entry.items()))

            if key.startswith("metal"):
                self.metals.append(MetalSpec(key, value))
            elif key.startswith("ligand"):
                self.ligands.append(LigandSpec(key, value))
            else:
                raise ValueError(f"Fatal Error: Unexpected key '{key}' in geometry")


class AssembledIsomer(RCA_Molecule):
    def __init__(self, ligands: Dict[int, RCA_Ligand], target_vectors: List[Dict[Any, List[float]]], ligand_origins: List[List[float]], metal_origins: List[List[float]],
                 metal_types: List[str], monometallic: bool = False):
        """
        Generates novel transition metal complexes from ligands and metals
        :param ligands:         Dictionary of ligand objects
        :param target_vectors:  List of dictionaries containing target vectors for each ligand
        :param ligand_origins:  List of ligand origins
        :param metal_origins:   List of metal origins
        :param metal_types:     List of metal types
        :param monometallic:    Boolean flag for monometallic complexes
        """
        # Define the class variables
        self.ligands = ligands.values()  # for each ligand, get the ligand object
        self.target_vectors = target_vectors  # List of target vectors i.e [[0, 0, 0], [0, 0, 1], ...]
        self.ligand_origins = ligand_origins  # List of ligand_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_origins = metal_origins  # List of metal_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_types = metal_types  # List of metal types i.e ['Ru', 'Mn', ...]
        self.mono_metallic = monometallic  # Boolean flag for monometallic complexes


        # Assemble the complex
        self.assembled_complexes = self._assemble_complex_new_new()

        self._tmc_validity_checks()

    def _tmc_validity_checks(self) -> None:
        """Some short checks specifically for transition metal complexes."""
        # Doublecheck if all the metals are really metals. Don't raise an error in case it's intentional.
        for metal_idx in self.metal_idc:
            is_metal = DART_Element(self.atomic_props['atoms'][metal_idx]).is_metal
            if not is_metal:
                warnings.warn(f"Metal center is not a metal.")

        return

    def to_dict(self):
        """
        Converts the AssembledIsomer object to a dictionary.
        :return: Dictionary representation of the AssembledIsomer object
        """
        data =  {
            "atomic_props": self.atomic_props,
            "graph": graph_to_dict_with_node_labels(self.graph),
            "global_props": self.global_props,
            "metal_idc": self.metal_idc,
            "donor_idc": self.donor_idc,
            "ligand_idc": self.ligand_idc,
            "ligand_info": self.ligand_info,
        }

        return data

    def get_metal_center_atoms(self) -> ase.Atoms:
        """
        Get the metal atoms of the complex.
        :return: ASE Atoms object containing the metal atoms
        """
        atoms = ase.Atoms()
        for metal_idx in self.metal_idc:
            atoms += self.mol[metal_idx]
        return atoms

    @classmethod
    def from_json(cls, filepath) -> 'AssembledIsomer':
        """
        Loads an AssembledIsomer object from a JSON file.
        :param filepath: Path to the JSON file
        :return: AssembledIsomer object
        """
        data = load_json(filepath)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'AssembledIsomer':
        """
        Creates an AssembledIsomer object from a dictionary in the correct format.
        :param data: Dictionary containing the AssembledIsomer data
        :return: AssembledIsomer object
        """
        data['graph'] = graph_from_graph_dict(data['graph'])
        return cls(**data)

    @classmethod
    def from_ligands_and_metal_centers(
                                        cls,
                                        ligands: List[RCA_Ligand],
                                        target_vectors: List[List[float]],
                                        metal_centers: Union[List[List[Union[str, List[float]]]], str],
                                        ligand_origins: List[List[float]] = None,
                                        ) -> Tuple[List['AssembledIsomer'], List[str]]:
        """
        Generates novel transition metal complexes from ligands and metal centers.
        :param ligands: List of RCA_Ligand objects from the MetaLig database.
        :param target_vectors: List of target vectors for each ligand.
        :param metal_centers: List of tuple with element and position for each metal center. If a string is provided, it is assumed to be the chemical element of a mono-metallic complex at the origin.
        :param ligand_origins: List of the origin for each ligand.

        Example usage for assembling a bi-metallic complex with three monodentate ligands, one of them bridging:
        target_vectors = [
                            [[1, 0, 0]],
                            [[0, 0, 1]],
                            [[-1, 0, 0]],
                         ]
        ligand_origins = [
                            [0, 0, 0],
                            [0, 0, 0],
                            [0, 0, 0]
                            ] # if all are at the origin like here, the `ligand_origins` can also be omitted
        ru = ['Ru', [1, 0, 0]]
        fe = ['Fe', [-1, 0, 0]]
        metal_centers = [
                            [ru],       # metal center for the first ligand
                            [ru, fe],   # metal centers for the second, bridging ligand
                            [fe]        # metal center for the third ligand
                        ]
        isomers, warnings = AssembledIsomer.from_ligands_and_metal_centers(
                                    ligands=ligands,
                                    target_vectors=target_vectors,
                                    ligand_origins=ligand_origins,
                                    metal_centers=metal_centers
                                    )
        """
        # Handle default values
        if ligand_origins is None:
            ligand_origins = [[0.0, 0.0, 0.0] for _ in ligands]
        if isinstance(metal_centers, str):
            # If the metal center is provided as a chemical element, it's a mono-metallic complex at the origin
            metal_centers = [[ase.Atom(symbol=metal_centers, position=[0, 0, 0])] for _ in ligands]
        else:
            # If the metal center is provided as a list of elements and positions, convert to ASE Atoms objects
            metal_centers = [[ase.Atom(symbol=metal[0], position=metal[1]) for metal in metal_list] for metal_list in metal_centers]

        # Check input format
        all_same_length = len(ligands) == len(target_vectors) == len(ligand_origins) == len(
            metal_centers)
        if not all_same_length:
            raise ValueError(
                "The input of ligands, target vectors, ligand origins and metal centers must have the same length.")

        rotated_ligands = get_rotated_ligands(ligands=ligands, target_vectors=target_vectors, ligand_origins=ligand_origins)

    def _assemble_complex_new_new(self) -> List[Atoms]:
        # --- Step 1: Partition ligands based on the shape of the target vectors ---
        ligand_lists = self._assign_ligands_to_vectors(
            ligands=list(self.ligands),
            vectors=self.target_vectors
        )

        # --- Step 2: Loop through the swapped ligand permutations and generate the isomers ---
        all_isomers = []
        for ligand_swapped_combo in ligand_lists:
            aligned_ligands = []  # This will contain one complexes ligands [l1, l2, [l3-i1, l3-i2], l4, ...]
            for idx, (ligand, ligand_target_vectors, origin) in enumerate(zip(ligand_swapped_combo, self.target_vectors, self.ligand_origins)):
                # Retrieve the ligand's geometry and donor atom indices.
                geometry, donor_atoms = ligand.get_isomers_effective_ligand_atoms_with_effective_donor_indices()

                # Set the tags for the ligand atoms
                # geometry.set_tags([ligand.elcn for _ in range(len(geometry))])
                geometry.new_array('multi_tags', np.full((len(geometry), 2), [ligand.elcn, idx + 1], dtype=int))

                # Convert target vector dictionary values to numpy arrays.
                target_vectors = [np.array(v) for v in ligand_target_vectors.values()]

                # Align the donor atoms of the ligand to the target vectors
                # todo need to revise. This code currently does not depend on the donor atoms but this functionality needs to be made more robust
                # todo currently too many isomers are being generated for pentagonal, trigonal, etc.
                ligand_isomers, donor_atoms_ordered, rssd = try_all_geometrical_isomer_possibilities(atoms=geometry,
                                                                                                     donor_idc=donor_atoms[0],
                                                                                                     target_vectors=target_vectors)

                # Translate the ligand to its correct location in the complex
                for isomer in ligand_isomers:
                    isomer.set_positions(isomer.get_positions() + np.array(origin))

                # Remove dummy atoms (e.g., "Cu") from haptic ligands.
                cleaned_isomers = self._remove_haptic_dummy_atom(
                    atoms_list=ligand_isomers,
                    dummy_atom="Cu",
                    donor_atoms_idc=ligand.hapdent_idc
                )
                # Append the rotated ligands to the list
                aligned_ligands.append(cleaned_isomers)

            all_isomers.append(aligned_ligands)

        # --- Step 3: Assemble the complexes using the ligand permutations ---
        all_isomers = self._gen_all_isomers(all_isomers)

        # --- Step 4: optimize the rotation of each mono-coordinating ligand around their respective coordination axis simultaneously ---
        optimizer = AxialOpt(complexes=all_isomers, target_vectors=self.target_vectors, ligand_origins=self.ligand_origins)
        all_isomers = optimizer.opt_mono_rotation()
        all_isomers = ReduceIsomers(all_isomers, metal_centres=self.metal_origins, rssd_threshold=2, debug=True).get_unique_isomers()

        return all_isomers

    def _flatten(self, nested_list: list):
        """
        Helper function to flatten a nested list
        :param nested_list:
        :return a generator object that can be used to iterate over the flattened list
        """
        for item in nested_list:
            if isinstance(item, list):
                yield from self._flatten(item)
            else:
                yield item

    def _assign_ligands_to_vectors(self, ligands: List[RCA_Ligand], vectors: List[Dict[str, List[float]]]) -> List[List[RCA_Ligand]]:
        """
        Assigns groups of ligands to vector entries based on the number of keys in each vector dictionary.
        For each unique key count (in order of first appearance in the vectors list), this function
        partitions the ligand list into contiguous groups. Each group’s size is determined by the number
        of vector entries that have that key count. Then, every vector entry is mapped to the corresponding
        ligand group.

        For example, given:
          ligands = ["L1", "L2", "L3", "L4"]
          vectors = [
              {'a': [0, 1]},                      # key count 1
              {'b': [1, 0]},                      # key count 1
              {'c': [1, 1], 'd': [2, 2]},           # key count 2
              {'e': [3, 3], 'f': [4, 4], 'g': [5, 5]}  # key count 3
          ]
        The unique key counts are 1, 2, and 3 with frequencies 2, 1, and 1 respectively. Thus,
        the ligand groups would be:
          Group for key count 1: ["L1", "L2"]
          Group for key count 2: ["L3"]
          Group for key count 3: ["L4"]
        and the output would be:
          [['L1', 'L2'], ['L1', 'L2'], ['L3'], ['L4']]

        :param: ligands (list): A list of ligand identifiers.
        :param: vectors (list of dict): A list of dictionaries, each representing a vector with arbitrary keys.

        :return: list: A list of ligand groups corresponding to each vector entry.

        :raise: ValueError: If there are not enough ligands to assign to all vector groups.
        """
        # Step 1: Determine the order of unique key counts and the frequency of each.
        unique_key_counts = []  # Order in which each unique key count first appears.
        frequency_by_key_count = {}  # Frequency of each key count.

        for vec in vectors:
            key_count = len(vec)
            if key_count not in frequency_by_key_count:
                unique_key_counts.append(key_count)
                frequency_by_key_count[key_count] = 0
            frequency_by_key_count[key_count] += 1

        # Verify that we have enough ligands.
        total_required_ligands = sum(frequency_by_key_count.values())
        if total_required_ligands > len(ligands):
            raise ValueError("Not enough ligands provided to assign to all vector entries.")

        # Step 2: Partition the ligand list into groups based on the frequency of each unique key count.
        ligand_groups = {}  # Maps key count to its corresponding ligand group.
        start_index = 0
        for key_count in unique_key_counts:
            group_size = frequency_by_key_count[key_count]
            ligand_groups[key_count] = ligands[start_index:start_index + group_size]
            start_index += group_size

        # Step 3: Create the output mapping: for each vector entry, select the ligand group that corresponds
        # to its number of keys.
        _list = [ligand_groups[len(vec)] for vec in vectors]

        results = []
        current_permutation = []
        used = set()

        def backtrack(index: int):
            # When all positions have been assigned, store the permutation.
            if index == len(_list):
                results.append(current_permutation.copy())
                return

            # Iterate over allowed ligands for the current position.
            for ligand in _list[index]:
                if ligand in used:
                    continue  # Skip if already used.
                # Choose this ligand.
                used.add(ligand)
                current_permutation.append(ligand)
                backtrack(index + 1)
                # Backtrack: remove the ligand and mark it as available.
                current_permutation.pop()
                used.remove(ligand)

        backtrack(0)
        return results

    def _add_metals(self, ligand_structure: Atoms):
        """
        Adds the metals to the complex
        """
        for metal_type, metal_origin in zip(self.metal_types, self.metal_origins):
            metal = Atoms(symbols=metal_type, positions=[metal_origin])
            metal.new_array('multi_tags', np.full((len(metal), 2), [0, 0], dtype=int))
            ligand_structure += metal
        return ligand_structure

    def _gen_all_isomers(self, ligands: List[Any]):
        """
        Generate all possible isomers from a list of ligands which have multiple isomers
        :param ligands: list: [[ligand1_isomer1, Ligand1_isomer2], [ligand2_isomer1, ligand2_isomer2], [ligand3_isomer1, ligand3_isomer2, ligand3_isomer3], ...]
        :return: list of ase objects
        """
        isomers = []
        for ligand_lists in ligands:
            # Generate all combinations; each combination is a tuple with one isomer per ligand.
            combinations = list(itertools.product(*ligand_lists))
            for combo in combinations:
                combined = Atoms()  # Start with an empty Atoms object.
                for ligand in combo:  # Iterate over the ligands in the combination.
                    combined += ligand  # combining Atoms objects.
                combined = self._add_metals(ligand_structure=combined)  # Add the metals to the complex
                isomers.append(combined)  # Store all the new isomers

        # Warnings for each isomer. If an isomer has no issues, the note should be ''. If an isomer is excluded because of clashing ligands or symmetrical ligands, the note should be `clashing_ligands' or `symmetrical_ligands`. todo @Cian: Please add corresponding warnings here after implementing the checks.
        warnings = ['' for _ in isomers]

        return isomers, warnings


class AxialOpt:
    def __init__(self, complexes: List[Atoms], target_vectors: List[Dict[str, List[float]]], ligand_origins: List[List[float]]):
        """
        This class will take a list of ASE Atoms objects and optimize mono-coordinating ligands around their coordination axis
        """
        self.input_complexes = complexes
        self.target_vectors = target_vectors
        self.ligand_origins = ligand_origins
        self.output_complexes = []

    def tag_mono_coordinating_ligands(self):
        """
        Tag the mono-coordinating ligands in the complex
        """
        for tmc in self.input_complexes:
            pass

    def opt_mono_rotation(self, opt_target: str = "max_atomic_distance"):
        """
        Optimize the rotation of each mono-coordinating ligand around their respective coordination axis simultaneously
        """
        np.random.seed(42)  # todo don't remove this as it is important to make sure the optimizer is deterministic
        # Specify bounds for each parameter: (lower, upper) for x and y respectively.
        bounds = [[0, 360] for _ in self.target_vectors]
        # Loop through each of the inputted complexes
        for tmc in self.input_complexes:

            # Run the global optimizer.
            result = differential_evolution(self.objective_function, bounds=bounds, args=(self.target_vectors, self.ligand_origins, tmc))

            # Retrieve the multi_tags array
            multi_tags = tmc.get_array("multi_tags")

            # Get the unique set of all second tags
            unique_ligand_idc_set = np.unique(multi_tags[:, 1])
            unique_ligand_idc_set = unique_ligand_idc_set[unique_ligand_idc_set != 0]

            # Loop through the unique second tags
            for angle, axis, origin, tag in zip(list(result.x), self.target_vectors, self.ligand_origins, unique_ligand_idc_set):
                # Get indices where the second tag is the current tag
                indices = np.where(multi_tags[:, 1] == tag)[0]

                # Check if any of these indices have a first tag not equal to 1, if so, skip
                if np.any(multi_tags[indices, 0] != 1):
                    continue

                tmc = self.rotate(atoms=tmc, vector=np.array(list(axis.values())[0]), origin=np.array(origin), idc=indices, angle=angle).copy()

            self.output_complexes.append(tmc)
        return self.output_complexes

    def objective_function(self, x: np.ndarray, vectors_in: List[np.array], origins_in: List[np.array], TMC_in: Atoms, ):
        """
        Objective function to optimize the position of the ligands in the TMC complex.
        :param x:
        :param vectors_in:
        :param origins_in:
        :param TMC_in:
        :return:
        """
        # Generate a copy of the input complex
        TMC_worker = TMC_in.copy()

        # Retrieve the multi_tags array
        multi_tags = TMC_worker.get_array("multi_tags")

        # Get the unique set of all nonzero second tags (each tag represents a ligand, the zero tag represents the metals).
        unique_tags = [tag for tag in np.unique(multi_tags[:, 1]) if tag != 0]

        for tag, angle, axis, origin in zip(unique_tags, list(x), vectors_in, origins_in):

            # Get indices where the second tag is the current tag (essentially the indices of the atoms in this particular ligand)
            indices = np.where(multi_tags[:, 1] == tag)[0]

            # Check if any of these indices have a first tag not equal to 1, if so, skip
            if np.any(multi_tags[indices, 0] != 1):
                continue  # Skip this ligand group

            TMC_worker = self.rotate(atoms=TMC_worker, vector=np.array(list(axis.values())[0]), origin=np.array(origin), idc=indices, angle=angle).copy()

        # Get the interatomic distance matrix
        distance_matrix = TMC_worker.get_all_distances()

        # Set the diagonal to a large number to avoid self-interaction (or np.inf)
        np.fill_diagonal(distance_matrix, np.inf)

        # Calculate the penalty: for each pair with d <= 4, add 1/d^2 to the penalty
        penalty = np.sum(1.0 / (distance_matrix ** 2))

        return penalty

    @staticmethod
    def rotate(atoms: Atoms, vector: np.array, origin: np.array, idc: List[int], angle: int):
        """
        Rotate the atoms in the Atoms object (only atoms with indices=idc) around the vector by the specified angle.
        :param atoms: Atoms object to rotate.
        :param vector: vector to rotate around.
        :param origin: origin of the rotation.
        :param idc: indices of the atoms to rotate.
        :param angle: the angle to rotate the atoms by in degrees.
        :return: an ase.Atoms object with the rotated atoms.
        """

        # Normalize rotation vector
        vector = np.asarray(vector, dtype=float)
        vector /= np.linalg.norm(vector)

        # Create rotation object
        rotation = R.from_rotvec(np.radians(angle) * vector)

        # Copy the atoms object to avoid modifying the original
        rotated_atoms = atoms.copy()

        # Apply rotation to selected atoms
        for i in idc:
            pos = atoms.positions[i] - origin  # Translate to origin
            rotated_pos = rotation.apply(pos) + origin  # Rotate and translate back
            rotated_atoms.positions[i] = rotated_pos

        return rotated_atoms


class ReduceIsomers:
    """
    This class will take a list of ASE Atoms objects and reduce the number of isomers to the minimum number of unique isomers
    """

    def __init__(self, isomers: List[Atoms], metal_centres, rssd_threshold: float = 0.01, debug: bool = False):
        """
        Initializes the ReduceIsomers object
        :param: isomers: List of ASE Atoms objects
        :param: rssd_threshold: The threshold for the RSSD value
        :param: metal_centres: The metal centres of the complex
        """
        self.isomers = isomers
        self.rssd_threshold = rssd_threshold
        self.metal_centres = metal_centres

        # debug attributes
        self.debug = debug
        self.brute_force_grid_density = 9

    def get_unique_isomers(self) -> List[Atoms]:
        """
        Get the unique isomers
        :return: List[Atoms]
        """
        return self._reduce_isomers()

    def _reduce_isomers_debug(self, write_svg: bool = True, plot_plotly: bool = True) -> List[Atoms]:
        """
        Reduce the number of isomers and display a similarity matrix.

        :param write_svg: If True, exports a static Matplotlib SVG version of the similarity matrix.
        :param plot_plotly: If True, displays the Plotly heatmap.
        :return: List[Atoms]
        """
        n = len(self.isomers)
        sim_matrix = np.zeros((n, n))

        # Show fingerprint difference matrix first (unchanged)
        self.plot_fingerprint_difference_matrix(self.isomers)

        # Compute similarity scores
        for i, isomer1 in enumerate(self.isomers):
            for j, isomer2 in enumerate(self.isomers):
                print(f"Comparing isomers [{i}] and [{j}] ... [{i * n + j + 1}/{n * n}]")
                score = self.align_isomers(isomer1, isomer2)
                sim_matrix[i, j] = score

        # Create labeled DataFrame
        labels = [f"{i}" for i in range(n)]
        sim_df = pd.DataFrame(sim_matrix, index=labels, columns=labels)

        # Custom color scale: red-white-blue
        custom_colorscale = [
            [0.0, "rgb(238,100,97)"],
            [0.5, "rgb(255,255,255)"],
            [1.0, "rgb(12,171,185)"]
        ]

        if write_svg:
            # Define Matplotlib colormap
            cmap = LinearSegmentedColormap.from_list(
                "custom_rwb",
                [(0.0, (238 / 255, 100 / 255, 97 / 255)),
                 (0.5, (1.0, 1.0, 1.0)),
                 (1.0, (12 / 255, 171 / 255, 185 / 255))]
            )

            fig, ax = plt.subplots(figsize=(8, 8))
            cax = ax.imshow(sim_df.values, cmap=cmap, vmin=0, vmax=1)

            # Annotate with values
            for i in range(n):
                for j in range(n):
                    val = sim_df.values[i, j]
                    r, g, b = cmap(val)[:3]
                    luminance = 0.299 * r + 0.587 * g + 0.114 * b
                    text_color = 'black' if luminance > 0.5 else 'white'
                    ax.text(j, i, f"{val:.2f}", va='center', ha='center', color=text_color, fontsize=8)

            ax.set_xticks(np.arange(n))
            ax.set_yticks(np.arange(n))
            ax.set_xticklabels(labels, rotation=90)
            ax.set_yticklabels(labels)
            ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
            fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04, label="Similarity Score")
            ax.set_title("Similarity Matrix of Isomers based on rotation", pad=20)
            plt.tight_layout()
            plt.savefig("similarity_matrix.svg", format="svg")
            plt.close()

        if plot_plotly:
            fig_sim = px.imshow(
                sim_df,
                text_auto=".2f",
                color_continuous_scale=custom_colorscale,
                title="Similarity Matrix of Isomers based on rotation"
            )
            fig_sim.show()

        return self.isomers

    def _reduce_isomers(self) -> List[Atoms]:
        """
        Takes a list of 'isomers' and identifies which are super-imposable and removes duplicates.
        In debug mode, it generates and displays a similarity matrix and a corresponding runtime matrix.
        :return: List[Atoms]
        """
        unique_isomers = []
        if self.debug:
            return self._reduce_isomers_debug()

        else:
            for isomer in self.isomers:
                # Check if this isomer is already in the list within the threshold.
                if not any(self.align_isomers(isomer, u) < self.rssd_threshold for u in unique_isomers):
                    unique_isomers.append(isomer)

        return unique_isomers

    def align_isomers(self, stationary_atoms: ase.Atoms, rotated_atoms: ase.Atoms) -> float:
        """
        Using Scipy's global optimizer to find the optimal rotation that minimizes the distance between like atoms in the two isomers
        when the two isomers are optimally aligned the rmsd or another measure of inter-atomic distance can be calculated.
        This method with a large enough grid density is considered to be the ground truth for the optimal rotation of the two isomers.
        other Scipy global optimizers can be easily integrated here for testing purposes.
        :param: stationary_atoms: ase.Atoms object
        :param: rotated_atoms: ase.Atoms object
        :return: float
        """
        metal_center = np.array(self.metal_centres[0])  # Ensure it's a NumPy array
        # Here the number of metal centres and the fact that metal centres of different isomers must always be aligned
        # is taken into account to determine if the isomers are similar or not
        if len(self.metal_centres) == 1:
            # There are 3 axis which an isomer can be rotated around to align with another isomer
            bounds = [[0, 360] for _ in range(3)]
            # The 3 cardinal axes, properly centered on the metal centre
            axes = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # Standard Cartesian axes
            # Center the axes around the metal center
            axes = metal_center + axes - metal_center  # This is effectively just axes
        elif len(self.metal_centres) == 2:
            # The isomer can only be rotated around the metal-metal axis to determine if the isomers are similar or not
            bounds = [[0, 360] for _ in range(1)]
            axes = np.linalg.norm(np.array(self.metal_centres[0]) - np.array(self.metal_centres[1]))

        elif len(self.metal_centres) >= 3:
            # 3 metal centres means each isomer is fixed in space and their geometries can be directly compared
            # We simply return the energy heuristic
            return self.energy_heuristic(stat_atoms=stationary_atoms, rot_atoms=rotated_atoms, mode="max")
        else:
            raise ValueError("The number of metal centres must be greater than or equal to 1")

        result2 = brute(self.objective_function, ranges=bounds, args=(stationary_atoms, rotated_atoms, axes), Ns=self.brute_force_grid_density)
        return self.objective_function(np.array(result2), stationary_atoms, rotated_atoms, axes)

    @staticmethod
    def combined_rotation_matrix(angles, axes):
        """
        Compute the combined rotation matrix from a list of angles and corresponding axes.
        :param angles: List of angles in degrees.
        :param axes: List of rotation axes.
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
        :param atoms: ASE Atoms object to rotate.
        :param R_total: Combined rotation matrix.
        :param center: Center of rotation
        """
        # Shift positions relative to center, apply rotation, then shift back.
        shifted = atoms.positions - center
        atoms.positions = center + (shifted @ R_total.T)

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
                                    rot_atoms=rotated_isomer,
                                    mode="max")
        return val

    def energy_heuristic(self, stat_atoms: ase.Atoms, rot_atoms: ase.Atoms, mode: str = "sum"):
        """
        Heuristic to determine the energy of the isomer to be minimized.

        For each element type, pair each atom in stat_atoms with an atom in rot_atoms
        (both having the same number of atoms for that element) such that:
          - If mode="sum": the sum of distances is minimized.
          - If mode="max": the maximum paired distance is returned.

        :param stat_atoms: ase.Atoms object (stationary isomer)
        :param rot_atoms: ase.Atoms object (rotated isomer)
        :param mode: 'sum' to return the sum of paired distances, or 'max' to return the maximum distance.
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

            if mode == "sum":
                total_energy += distance_matrix[row_ind, col_ind].sum()
            elif mode == "max":
                element_max = distance_matrix[row_ind, col_ind].max()
                overall_max_distance = max(overall_max_distance, element_max)
            else:
                raise ValueError("Invalid mode. Use 'sum' or 'max'.")

        return total_energy if mode == "sum" else overall_max_distance

    @staticmethod
    def compute_sorted_distance_fingerprint(isomer) -> np.ndarray:
        """
        Compute an inter-atomic distance matrix for an isomer and sort the distances under to two conditions:
        1. entries in the matrix (atom-atom distances) are sorted in ascending order.
        2. the rows of the matrix are sorted in lexicographical order of the element symbols.
            i.e. C-C, C-H, H-H, etc. This ensures like interatomic distances are compared between isomers.
        The resulting 1D array acts as a fingerprint for the isomer which can be compared to other isomers.

        :param isomer: ASE Atoms object.
        :return: sorted inter-atomic distance fingerprint as a 1D numpy array.
        """

        # Get the positions and symbols of the atoms in the isomer.
        positions = isomer.get_positions()  # shape (N, 3)
        symbols = isomer.get_chemical_symbols()  # list of element symbols

        # Dictionary to collect distances keyed by element pair, e.g. ('C', 'H')
        pair_distances = {}

        # Compute distances only for the upper triangle (i < j)
        for i in range(len(symbols)):
            for j in range(i + 1, len(symbols)):
                # Create a key that is independent of the order of the atoms i.e ('C', 'H') == ('H', 'C') becomes ('C', 'H')
                key = tuple(sorted((symbols[i], symbols[j])))
                d = np.linalg.norm(positions[i] - positions[j])
                # Checks if the key exists if not the key is created and the distance is appended to the list
                pair_distances.setdefault(key, []).append(d)

        # For a consistent fingerprint, sort the distances for each element pair
        # and then concatenate them in a fixed order (e.g., lexicographical order of the keys)
        fingerprint_parts = []
        for key in sorted(pair_distances.keys()):
            sorted_group = np.sort(pair_distances[key])

            fingerprint_parts.extend(sorted_group.tolist())
        fingerprint = np.array(fingerprint_parts)

        # Return the fingerprint as a 1D numpy array. i.e. [d1, d2, d3, ...]
        return fingerprint

    @staticmethod
    def fingerprint_comparison(fp1: np.ndarray, fp2: np.ndarray, mode: str = "max_diff"):
        """
        Compute the maximum absolute difference between two sorted fingerprint vectors.

        :param fp1: 1D numpy array fingerprint for isomer 1.
        :param fp2: 1D numpy array fingerprint for isomer 2.
        :param mode: Comparison mode: 'max_diff', 'sum_diff', 'mean_diff', or 'rmsd'.
        :return: Maximum absolute difference.
        """
        if mode == "max_diff":
            return np.max(np.abs(fp1 - fp2))
        elif mode == "sum_diff":
            return np.sum(np.abs(fp1 - fp2))
        elif mode == "mean_diff":
            return np.mean(np.abs(fp1 - fp2))
        elif mode == "rmsd":
            return np.sqrt(np.mean((fp1 - fp2) ** 2))
        else:
            raise ValueError("Error: invalid comparison mode. Use 'max_diff', 'sum_diff', 'mean_diff', or 'rmsd'.")

    def compute_fingerprint_matrix(self, isomers: list) -> np.ndarray:
        """
        Compute the pairwise fingerprint difference matrix between isomers,
        where the difference between two isomers is defined as the maximum absolute
        difference between their sorted pairwise distance fingerprints.

        :param isomers: List of ASE Atoms objects.
        :return: 2D numpy array representing the pairwise differences.
        """
        # Precompute the fingerprint for each isomer.
        fingerprints = [self.compute_sorted_distance_fingerprint(isomer) for isomer in isomers]
        n = len(fingerprints)
        diff_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                diff_matrix[i, j] = self.fingerprint_comparison(fingerprints[i], fingerprints[j], mode="max_diff")
        return diff_matrix

    @staticmethod
    def _cluster_fingerprint_matrix(matrix: np.ndarray, quantile: float = 0.2) -> np.ndarray:
        """
        Cluster the fingerprint difference matrix values using MeanShift and return a matrix of group labels.
        The cluster whose center is closest to zero is labeled as "Close", all others are labeled as "Far".

        :param: matrix: 2D numpy array representing the fingerprint difference matrix.
        :param: quantile (float): Quantile for estimating the bandwidth (default is 0.2).
        :return: np.ndarray: An array of the same shape as 'matrix' containing group labels ("Close" or "Far").
        """
        # Flatten the matrix values for clustering.
        values = matrix.flatten().reshape(-1, 1)
        bandwidth = estimate_bandwidth(values, quantile=quantile, n_samples=len(values))
        ms = MeanShift(bandwidth=bandwidth)
        ms.fit(values)
        cluster_labels = ms.labels_
        cluster_centers = ms.cluster_centers_

        # Identify the cluster whose center is closest to zero.
        close_cluster_idx = np.argmin(np.abs(cluster_centers))
        group_labels = np.where(cluster_labels == close_cluster_idx, "Close", "Far")

        # Reshape back to the original matrix shape.
        return group_labels.reshape(matrix.shape)

    def plot_fingerprint_difference_matrix(self, isomers: List[Atoms], write_svg: bool = True, plot_plotly: bool = True) -> None:
        """
        Computes the fingerprint difference matrix for a list of isomers and plots it as a heatmap and histogram.
        The heatmap is annotated with both the numeric difference (centered) and a group tag (offset).
        Clustering is performed using MeanShift—via a helper method—and the cluster whose center is closest
        to zero is labeled as "Close", with all others labeled as "Far".

        :param: isomers (List[ase.Atoms]): A list of ASE Atoms objects.
        :param: write_svg (bool): If True, exports a static Matplotlib SVG version of the matrix instead of showing the Plotly heatmap.
        :param: plot_plotly: If True, shows the Plotly heatmap.
        """
        # Compute the fingerprint difference matrix.
        diff_matrix = self.compute_fingerprint_matrix(isomers)
        labels_list = [f"{i}" for i in range(len(isomers))]
        df = pd.DataFrame(diff_matrix, index=labels_list, columns=labels_list)

        # Cluster the fingerprint differences.
        group_labels_matrix = self._cluster_fingerprint_matrix(df.values, quantile=0.2)

        # Red-white-blue custom colour scale
        custom_colorscale = [
            [0.0, "rgb(238,100,97)"],  # Red
            [0.5, "rgb(255,255,255)"],  # White
            [1.0, "rgb(12,171,185)"]  # Blue
        ]

        if write_svg:
            # Define the custom matplotlib colormap
            cmap = LinearSegmentedColormap.from_list(
                "custom_rwb",
                [(0.0, (190 / 255, 140 / 255, 215 / 255)),
                 (0.5, (1.0, 1.0, 1.0)),
                 (1.0, (0 / 255, 186 / 255, 124 / 255))]
            )

            fig, ax = plt.subplots(figsize=(8, 8))
            cax = ax.imshow(df.values, cmap=cmap, vmin=0, vmax=1.0)

            def get_contrasting_text_color(value, vmin=0.0, vmax=0.5, cmap=None):
                """Return 'white' or 'black' based on background luminance."""
                norm_val = (value - vmin) / (vmax - vmin)
                r, g, b = cmap(norm_val)[:3]
                luminance = 0.299 * r + 0.587 * g + 0.114 * b  # perceived brightness
                return 'black' if luminance > 0.5 else 'white'

            # Annotate cells with numeric values and group tags
            nrows, ncols = df.shape
            for i in range(nrows):
                for j in range(ncols):
                    # Get the contrasting text color
                    text_color = get_contrasting_text_color(df.values[i, j], cmap=cmap)
                    ax.text(j, i, f"{df.values[i, j]:.2f}", va='center', ha='center', color=text_color, fontsize=8)
                    # ax.text(j - 0.3, i - 0.3, group_labels_matrix[i, j], va='top', ha='left', color='red', fontsize=10) # Temporarily commented out for aesthetic reasons

            ax.set_xticks(np.arange(ncols))
            ax.set_yticks(np.arange(nrows))
            ax.set_xticklabels(labels_list, rotation=90)
            ax.set_yticklabels(labels_list)
            ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

            fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04, label="Difference")
            ax.set_title("Fingerprint Difference Matrix with Annotations", pad=20)
            plt.tight_layout()
            plt.savefig("fingerprint_difference_matrix.svg", format="svg")
            plt.close()

        if plot_plotly:
            # Plotly version for display
            fig = px.imshow(
                df,
                color_continuous_scale=custom_colorscale,
                zmax=0.5,
                title=("Fingerprint Difference Matrix with Annotations\n"
                       "Numeric value (center) and group tag (corner)")
            )

            # Numeric annotations
            num_annotations = []
            nrows, ncols = df.shape
            for i in range(nrows):
                for j in range(ncols):
                    num_annotations.append(dict(
                        x=j,
                        y=i,
                        text=str(round(df.values[i, j], 3)),
                        showarrow=False,
                        font=dict(color="white", size=10),
                        xanchor="center",
                        yanchor="middle"
                    ))

            # Group tag annotations
            group_annotations = []
            for i in range(nrows):
                for j in range(ncols):
                    group_annotations.append(dict(
                        x=j - 0.35,
                        y=i - 0.35,
                        text=group_labels_matrix[i, j],
                        showarrow=False,
                        font=dict(color="red", size=12),
                        xanchor="left",
                        yanchor="top"
                    ))

            fig.update_layout(annotations=num_annotations + group_annotations)
            fig.show()

        # Histogram (always shown)
        fig_hist = px.histogram(
            x=df.values.flatten(),
            nbins=400,
            title="Histogram of Fingerprint Difference Values",
            labels={'x': 'Difference Value', 'y': 'Count'}
        )
        fig_hist.show()

# todo: Notes for Cian (need to review as some of these are now outdated)
# todo: How the following method needs to be modified to account for swapping ligands of a same coordination mode (e.g. 1-1-1-1-1-1 or 2-2-2 or 3-3 (but the 3-3 case kind of doesn't apply))
# todo: 1. we need a measure of ligand coordination similarity. should probs use string: "monodentate", "trigonal", "tetragonal offset", etc.
# todo: 2. For example, lets say we have Ligands A, B, C and D
# todo: A: monodentate
# todo: B: monodentate
# todo: C: monodentate
# todo: D: tridentate
# todo: 3. We need to generate all possible isomers, enantiomers of the complex
# todo: 4. not only can the tridentate be flipped but the monodentate ligands can be exchanged as well
# todo: 5. to swap ligands we need to transform (rotate) them based on their origin and target vectors (translate based on vector between origins and rotate based on the angle between the target vectors)
# todo: 6. this needs to be done for every possible combination of "swappable" ligands




#########################################################################################
# This file contains the classes and methods that are used to process the input data    #
# and generate the assembled transition metal complexes                                 #
#########################################################################################
import warnings

from DARTassembler.src.assembly.ligand_geometries import try_all_geometrical_isomer_possibilities
from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.DataBase import RCA_Ligand, LigandDB
from DARTassembler.src.constants import Periodic_Table as PerTab
from typing import Dict, Any, List, Optional, Tuple, Union
from ase.visualize import view
from ase import Atoms
import numpy as np
import itertools
import ase

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
        self.ligand_db = LigandDB.load_from_json(path=self.ligand_db_path)
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
    REQUIRED_KEYS = {"name", "random_seed", "max_num_complexes", "total_charge", "geometry", "isomer_instruction"}

    def __init__(self, batch: Dict[str, Any]):
        """
        This class will parse the input of the assembly YAML file, store, check and return the information
        """

        self.batch = batch
        self._validate_batch()

        self.batch_name = self._get_batch_name()
        self.random_seed = self._get_random_seed()
        self.max_num_complexes = self._get_max_num_complexes()
        self.total_charge = self._get_total_charge()
        self.geometry = self._get_geometry()

        # Metals and ligands
        self.ligands = []
        self.metals = []
        self._process_geometry()

        # calculate the total metal oxidation state
        self.total_metal_oxidation_state = int(sum([metal.metal_oxidation_state for metal in self.metals]))

    def _validate_batch(self) -> None:
        """
        Validates the provided batch dictionary, ensuring all required keys exist
        :raises ValueError: If a required key is missing
        """
        missing_keys = self.REQUIRED_KEYS - self.batch.keys()
        extra_keys = self.batch.keys() - self.REQUIRED_KEYS  # Helpful for debugging unexpected keys

        if missing_keys:
            raise ValueError(f"Fatal Error: Missing required keys in input file: {', '.join(missing_keys)}")

        if extra_keys:
            print(f"Warning: Unrecognized keys found in batch: [{', '.join(extra_keys)}] These will be ignored")

    def _get_batch_name(self) -> str:
        """
        Retrieves and returns the batch name, ensuring it exists
        :return: The batch name as a string
        """
        return str(self.batch["name"])

    def _get_random_seed(self) -> float:
        """
        Retrieves and returns the random seed, ensuring it exists
        :return: The random seed as a float
        """
        try:
            return float(self.batch["random_seed"])
        except (KeyError, ValueError, TypeError):
            raise ValueError("Fatal Error: Random seed must be float or int")

    def _get_max_num_complexes(self) -> int:
        """
        Retrieves and returns the maximum number of complexes, ensuring it exists
        :return: The maximum number of complexes as an int
        """
        try:
            return int(self.batch["max_num_complexes"])
        except (KeyError, ValueError, TypeError):
            raise ValueError("Fatal Error: Maximum number of complexes must be an integer")

    def _get_total_charge(self) -> int:
        """
        Retrieves and returns the total charge, ensuring it exists
        :return: The total charge as an int
        """
        try:
            return int(self.batch["total_charge"])
        except (KeyError, ValueError, TypeError):
            raise ValueError("Fatal Error: Total charge must be an integer")

    def _get_geometry(self) -> list:
        """
        Retrieves and returns the geometry, ensuring it exists
        :return: The geometry as a string
        """
        geometry = self.batch.get("geometry")
        if not isinstance(geometry, list):
            raise ValueError("Fatal Error: 'geometry' must be a list of dictionaries")
        return geometry

    @staticmethod
    def _validate_metal(metal_data: Dict[str, Any]) -> List[float]:
        """
        Validates the metal's coordinate entry
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

    def get(self, key: str, default: any = None) -> any:
        """
        Safe retrieval of optional values from the batch dictionary
        :param key: Key to retrieve
        :param default: Default value if the key is missing
        :return: Value associated with the key or the default
        """
        return self.batch.get(key, default)


class AssembledIsomer(RCA_Molecule):

    def __init__(self,
                    atomic_props: Union[ase.Atoms, Dict[str, Any]],
                    graph: nx.Graph,
                    metal_idc: List[int],
                    donor_idc: List[List[int]],
                    ligand_idc: List[List[int]],
                    ligand_info: Dict[str, Any] = None,
                    global_props: Dict[str, Any] = None,
                    validity_check: bool = False,
                    ):
        if global_props is None:
            global_props = {}
        if ligand_info is None:
            ligand_info = {}

        super().__init__(
                         atomic_props=atomic_props,
                         global_props=global_props,
                         graph=graph,
                         validity_check=validity_check
                         )

        self.metal_idc = metal_idc
        self.donor_idc = donor_idc
        self.ligand_idc = ligand_idc
        self.ligand_info = ligand_info

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

        # Generate all combinations; each combination is a tuple with one isomer per ligand.
        combinations = list(itertools.product(*rotated_ligands))
        unique_metal_centers = get_all_unique_metal_centers(metal_centers)
        ase_isomers = []
        for combo in combinations:
            combined = Atoms()  # Start with an empty Atoms object.
            # todo @Cian: I always found it extremely handy in DART that the metal center came first, so I have changed the order here (in your code now, the ligand came first). Pls remove comment if ok. Otherwise, the `metal_idc` needs to be updated as well.
            for atom in unique_metal_centers:
                combined += atom
            for ligand in combo:  # Iterate over the ligands in the combination.
                combined += ligand  # combining Atoms objects.
            ase_isomers.append(combined)  # Store all the new isomers
        metal_idc = [idx for idx in range(len(unique_metal_centers))]

        isomers = []
        for ase_isomer in ase_isomers:
            # Merge the graphs of the ligands and the metal centers to get the full graph of the complex.
            graph, ligand_indices, donor_indices = get_merged_graph_from_ligands_and_metal_centers(ligands=ligands, metal_centers=metal_centers)
            global_props = {}   # Will be populated during the DART workflow.
            # To save memory, each complex is saved with only the most important information about its ligands.
            ligand_info = {
                'unique_names': [lig.unique_name for lig in ligands],
                'geometries': [lig.geometry for lig in ligands],
                'donors': ['-'.join(sorted(lig.local_elements)) for lig in ligands],
                'charges': [lig.pred_charge for lig in ligands],
                'stoichiometries': [lig.stoichiometry for lig in ligands],
            }
            isomer = cls(
                            atomic_props=ase_isomer,
                            graph=graph,
                            metal_idc=metal_idc,
                            ligand_idc=ligand_indices,
                            donor_idc=donor_indices,
                            global_props=global_props,
                            ligand_info=ligand_info,
                            )
            isomers.append(isomer)

        # Warnings for each isomer. If an isomer has no issues, the note should be ''. If an isomer is excluded because of clashing ligands or symmetrical ligands, the note should be `clashing_ligands' or `symmetrical_ligands`. todo @Cian: Please add corresponding warnings here after implementing the checks.
        warnings = ['' for _ in isomers]

        return isomers, warnings


class ReduceIsomers:
    """
    This class will take a list of ASE Atoms objects and reduce the number of isomers to the minimum number of unique isomers
    """

    def __init__(self, isomers: List[Atoms], rssd_threshold: float = 0.01):
        """
        Initializes the ReduceIsomers object
        :param isomers: List of ASE Atoms objects
        """
        raise NotImplementedError("Fatal Error !!! This class is not yet implemented and is currently under construction. Please do not use !!!")
        self.isomers = isomers
        self.rssd_threshold = rssd_threshold
        self.unique_isomers = self._reduce_isomers()

    def get_unique_isomers(self) -> List[Atoms]:
        """
        Get the unique isomers
        :return: List[Atoms]
        """
        return self.unique_isomers

    def _reduce_isomers(self) -> List[Atoms]:
        """
        takes a list of 'isomers' and identifies which are super-imposable and removes duplicates
        :return: List[Atoms]
        """
        unique_isomers = []
        """
        for isomer in self.isomers:
            # Check if this isomer is already in the list within the threshold
            if not any(self.align_molecules(isomer, u) < self.rssd_threshold for u in unique_isomers):
                unique_isomers.append(isomer)"""

        for i, isomer1 in enumerate(self.isomers):
            for j, isomer2 in enumerate(self.isomers):
                view(isomer1)
                view(isomer2)
                print(i, j, self.align_molecules(isomer1, isomer2))
                print("stop")
        return unique_isomers

    @staticmethod
    def align_molecules(
            atoms1: ase.Atoms,
            atoms2: ase.Atoms,
    ):
        """
        Computes the Root Mean Square Deviation (RMSD) between two ASE Atoms objects
        after optimally aligning A onto B using the Kabsch algorithm.

        Assumes both Atoms objects have the same number of atoms and are correspondingly ordered.

        Parameters:
        ase_atoms_A (Atoms): The Atoms object to be rotated (aligned).
        ase_atoms_B (Atoms): The reference Atoms object.

        Returns:
        float: The RMSD value after optimal alignment.
        """
        assert len(atoms1) == len(atoms2), "Both Atoms objects must have the same number of atoms."

        # todo: I am pretty sure the issue is that the OH is not perfectly the same across isomers

        # Extract coordinates
        A = atoms1.get_positions()
        B = atoms2.get_positions()

        # Compute centroids
        centroid_A = np.mean(A, axis=0)
        centroid_B = np.mean(B, axis=0)

        # Center the coordinates
        A_centered = A - centroid_A
        B_centered = B - centroid_B

        # Compute covariance matrix
        H = A_centered.T @ B_centered

        # Perform Singular Value Decomposition
        U, S, Vt = np.linalg.svd(H)

        # Compute the optimal rotation matrix
        R = Vt.T @ U.T

        # Ensure a proper rotation (avoiding reflection)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        # Apply rotation to A
        A_aligned = (A - centroid_A) @ R + centroid_B

        # Compute RMSD
        rmsd = np.sqrt(np.mean(np.sum((A_aligned - B) ** 2, axis=1)))

        return rmsd

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










# # Notes for talk with Cian
# # algorithm above? Is this the align_vectors method?
# # New refactoring:
# # old way:
# ChemBuild = AssembledIsomer(
#                             ligands=...,
#                             target_vectors=...,
#                             ligand_origins=...,
#                             metal_centers=...,
#                             )
# isomers = ChemBuild.get_isomers()   # returns a list of ase.Atoms() objects
# # To save these, one needs to somehow combine the ase.Atoms() object with the graph and some other data like the metal OS:
# for isomer in isomers:
#     graph = ...
#     metal_os = ...
#     save_dict = {"graph": graph, "metal_os": metal_os, "atoms": isomer}
#
# # new way:
# isomers = AssembledIsomer.from_ligands_and_metal_centers(   # isomers is a list of AssembledIsomer() objects
#                                                             ligands=...,
#                                                             target_vectors=...,
#                                                             ligand_origins=...,
#                                                             metal_centers=...,
#                                                             )
# # To save this, we can simply do the following:
# for isomer in isomers:
#     isomer.save_to_file("path/to/save/file.json")
# # Is also be much simpler to generalize since we can use the RCA_Molecule class as parent class for the AssembledIsomer class



#########################################################################################
# This file contains the classes and methods that are used to process the input data    #
# and generate the assembled transition metal complexes                                 #
#########################################################################################
from DARTassembler.src.assembly.ligand_geometries import try_all_geometrical_isomer_possibilities
from DARTassembler.src.ligand_extraction.DataBase import RCA_Ligand, LigandDB
from DARTassembler.src.constants import Periodic_Table as PerTab
from typing import Dict, Any, List, Optional, Tuple
from ase.visualize import view
from ase import Atoms
import numpy as np
import itertools
import ase



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


class AssemblyComplex(object):

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
        self.ligands = ligands.values()         # for each ligand, get the ligand object
        self.target_vectors = target_vectors    # List of target vectors i.e [[0, 0, 0], [0, 0, 1], ...]
        self.ligand_origins = ligand_origins    # List of ligand_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_origins = metal_origins      # List of metal_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_types = metal_types          # List of metal types i.e ['Ru', 'Mn', ...]
        self.mono_metallic = monometallic       # Boolean flag for monometallic complexes

        # Validate the input
        self._validate_input()

        # Assemble the complex
        self.assembled_complexes = self._assemble_complex()

    def _validate_input(self) -> None:
        """
        Validates the input to the AssemblyComplex class
        :raises ValueError: If the input is invalid
        """
        if len(self.ligands) != len(self.target_vectors) or len(self.ligands) != len(self.ligand_origins):
            raise ValueError(
                f"Fatal Error: Ligand objects [{len(self.ligands)}], target vectors [{len(self.target_vectors)}], and Ligand_Origins [{len(self.ligand_origins)}] must have the same length")

        if self.mono_metallic:
            if (len(self.metal_origins) != 1) or (len(self.metal_types) != 1):
                raise ValueError("Fatal Error: Monometallic complexes must have exactly one metal origin and type.")
        else:
            if len(self.metal_origins) != len(self.metal_types):
                raise ValueError(f"Fatal Error: Metal origins [length: {self.metal_origins}] and metal types [length: {self.metal_types}] must have the same length for multi-metallic systems.")

    def _assemble_complex(self) -> List[Atoms]:
        """
        Assembles the complex using the provided geometry and ligands
        :return: ASE Atoms object representing the complex
        """
        rotated_ligands = []
        for ligand, target_vectors, origin in zip(self.ligands, self.target_vectors, self.ligand_origins):
            # Extract the geometry and donor atoms of the ligand
            geometry, donor_atoms = ligand.get_isomers_effective_ligand_atoms_with_effective_donor_indices()
            # cast the target vectors to numpy arrays
            target_vectors = [np.array(v) for v in target_vectors.values()]
            # Align the donor atoms of the ligand to the target vectors
            ligand_isomers, donor_atoms_ordered, rssd = try_all_geometrical_isomer_possibilities(atoms=geometry,
                                                                                                 donor_idc=donor_atoms[0],
                                                                                                 target_vectors=target_vectors)
            # Remove the dummy atom from the haptic ligands
            ligand_isomers = self._remove_haptic_dummy_atom(atoms_list=ligand_isomers, dummy_atom="Cu", donor_atoms_idc=ligand.hapdent_idc)

            # Append the rotated ligands to the list
            rotated_ligands.append(ligand_isomers)

        # Use the rotated ligands to generate all (for now it is most) possible isomers
        all_isomers = self._gen_all_isomers(rotated_ligands)

        return all_isomers

    def _add_metals(self, ligand_structure: Atoms):
        """
        Adds the metals to the complex
        """
        for metal_type, metal_origin in zip(self.metal_types, self.metal_origins):
            metal = Atoms(symbols=metal_type, positions=[metal_origin])
            ligand_structure += metal
        return ligand_structure

    def _gen_all_isomers(self, ligands: List[List[Atoms]]):
        """
        Generate all possible isomers from a list of ligands which have multiple isomers
        :param ligands: list: [[ligand1_isomer1, Ligand1_isomer2], [ligand2_isomer1, ligand2_isomer2], [ligand3_isomer1, ligand3_isomer2, ligand3_isomer3], ...]
        :return: list of ase objects
        """
        # Generate all combinations; each combination is a tuple with one isomer per ligand.
        combinations = list(itertools.product(*ligands))
        isomers = []

        for combo in combinations:
            combined = Atoms()                                          # Start with an empty Atoms object.
            for ligand in combo:                                        # Iterate over the ligands in the combination.
                combined += ligand                                      # combining Atoms objects.
            combined = self._add_metals(ligand_structure=combined)      # Add the metals to the complex
            isomers.append(combined)                                    # Store all the new isomers

        return isomers

    def get_isomers(self) -> List[Atoms]:
        """
        Get the assembled complexes
        :return:
        """
        return self.assembled_complexes

    @staticmethod
    def _remove_haptic_dummy_atom(atoms_list: List[Atoms], dummy_atom: str, donor_atoms_idc: Tuple[Tuple[int]]):
        """
        Removes the dummy atom from the generated isomers
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
                dummy_idc.sort(reverse=True)  # This is important so that the larger index is removed first so as not to change the index of the other atoms
                for dummy_idx in dummy_idc:
                    atoms.pop(dummy_idx)
            return atoms_list


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

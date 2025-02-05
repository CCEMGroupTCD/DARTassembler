import ase
from scipy.spatial.transform import Rotation as R
from DARTassembler.src.ligand_extraction.DataBase import MoleculeDB, RCA_Ligand, LigandDB
from DARTassembler.src.constants import Periodic_Table as PerTab
from DARTassembler.src.assembly.ligand_geometries import try_all_geometrical_isomer_possibilities
from typing import Dict, Any, List, Optional, Union
from ase import Atoms
import itertools
import numpy as np


def align_donor_atoms(
        atoms: ase.Atoms,
        donor_idc: list[int],
        target_vectors: list[list[float]],
        origin: list[float],
        return_rssd: bool = False
):
    """
    Align the donor atoms of a ligand to the target vectors.
    :param return_rssd: #todo add description
    :param origin: The origin of the ligand.
    :param atoms: ASE Atoms object of all the atoms of the ligand.
    :param donor_idc: Indices of the donor atoms of the ligand.
    :param target_vectors: A list of 3D target vectors to which the donor atoms should be aligned.
    :return: ASE Atoms object of the ligand with the donor atoms aligned to the target vectors.
    """
    atoms = atoms.copy()  # Make a copy of the original atoms object so that the original is not changed
    target_vectors = np.array(target_vectors)
    assert len(target_vectors) == len(donor_idc), 'The number of target vectors must match the number of donor atoms.'
    assert target_vectors.shape[1] == 3, 'The target vectors must be 3D vectors.'
    donor_idc = list(donor_idc)  # A tuple wouldn't work for indexing

    donor_vectors = atoms.positions[donor_idc]
    # Normalize the donor vectors and target vectors to unit vectors so that only the direction of the vectors counts, not the magnitude.
    donor_vectors = donor_vectors / np.linalg.norm(donor_vectors, axis=1)[:, None]
    target_vectors = target_vectors / np.linalg.norm(target_vectors, axis=1)[:, None]
    # Find the correct rotation to align the donor vectors with the target vectors
    rot, rssd = R.align_vectors(a=target_vectors, b=donor_vectors)  # the a and b are unintuitive but correct
    # Apply the rotation to all the atoms of the ligand
    rotated_coords = rot.apply(atoms.positions - origin) + origin
    atoms.set_positions(rotated_coords)

    if return_rssd:
        return atoms, rssd
    else:
        return atoms


class LigandSpec:
    """
    Represents a ligand entry in the batch file's geometry section
    Handles different numbers of vectors dynamically
    """

    def __init__(self, name: str, data: Dict[str, Any]):
        """
        Initializes a Ligand object
        :param name: The ligand's identifier (e.g., 'ligand_1')
        :param data: Dictionary containing ligand properties
        :raises ValueError: If required keys are missing
        """
        self.name = name
        self.data = data
        self.origin = self._get_vector(key="origin", required=True)
        self.vectors = self._extract_vectors()
        self.ligand_db_path = data.get("ligand_db", None)
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
        Initializes a Metal object
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
        This class will parse the input of the assembly YAML file and store, check and return the information
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

    def __init__(self, ligands: Dict[int, MoleculeDB], target_vectors: List[Dict[Any, List[float]]], ligand_origins: List[List[float]], metal_origins: List[List[float]],
                 metal_types: List[str], monometallic: bool = False):
        """
        Generates novel transition metal complexes from ligands and metals
        :param ligands:
        :param target_vectors:
        :param ligand_Origins:
        :param metal_origins:
        :param metal_types:
        :param monometallic:
        """
        # Define the class variables
        self.ligands = ligands.values()         # for each
        self.target_vectors = target_vectors    # List of target vectors i.e [[0, 0, 0], [0, 0, 1], ...]
        self.ligand_origins = ligand_origins    # List of ligand_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_origins = metal_origins      # List of metal_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_types = metal_types          # List of metal types i.e ['Ru', 'Mn', ...]
        self.mono_metallic = monometallic

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
            ligand: RCA_Ligand
            ase_ligand = ligand.get_ase_molecule()
            _, donor_atom_idx = ligand.get_all_effective_atoms_with_effective_donor_indices()
            target_vectors_2 = [np.array(v) for v in target_vectors.values()]
            ligand_isomers, donor_atoms_ordered, rssd = try_all_geometrical_isomer_possibilities(atoms=ase_ligand,
                                                                                                 donor_idc=donor_atom_idx,
                                                                                                 target_vectors=target_vectors_2)
            rotated_ligands.append(ligand_isomers)

        # Generate all possible isomers
        # todo: How the following method needs to be modified to account for swapping ligands of a same coordination mode (e.g. 1-1-1-1-1-1 or 2-2-2 or 3-3 (but the 3-3 case kind of doesn't apply))
        # todo: 1. we need a measure of ligand coordination similarity. should probs use string: "monodentate", "trigonal", "tetragonal offset", etc.
        # todo: 2. For example, lets say we have Ligands A, B, C and D
            # todo: A: monodentate
            # todo: B: monodentate
            # todo: C: monodenate
            # todo: D: tridenate
        # todo: 3. We need to generate all possible isomers, enantiomers of the complex
        # todo: 4. not only can the tridentate be flipped but the monodentates can be exchaged as well
        # todo: 5. to swap ligands we need to transform (rotate) them based on their origin and target vectors (translate based on vector between origins and rotate based on the angle between the target vectors)
        # todo: 6. this needs to be done for every possible combination of "swappable" ligands
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
            combined = Atoms()  # Start with an empty Atoms object.
            for ligand in combo:
                combined += ligand  # Assuming the '+' operator works for combining Atoms objects.
            combined = self._add_metals(ligand_structure=combined)
            isomers.append(combined)

        return isomers

    def get_isomers(self) -> List[Atoms]:
        """
        Get the assembled complexes
        :return:
        """
        return self.assembled_complexes

#########################################################################################
# This file contains the classes and methods that are used to process the input data    #
# and generate the assembled transition metal complexes                                 #
#########################################################################################


# Dart Assembler imports
from DARTassembler.src.assembly.ligand_geometries import try_all_geometrical_isomer_possibilities
from DARTassembler.src.ligand_extraction.DataBase import RCA_Ligand, LigandDB, RCA_Molecule
from DARTassembler.src.constants import Periodic_Table as PerTab

# Scientific package imports
from scipy.optimize import differential_evolution, linear_sum_assignment, brute
from sklearn.cluster import MeanShift, estimate_bandwidth
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import cdist
import numpy as np

# Standard library imports
from typing import Dict, Any, List, Optional, Tuple, Union
from dev.Assembler_revision_June_2025.assembler.logger import default_logger as logger
from collections import defaultdict
from dev.Assembler_revision_June_2025.assembler.logger import LoggedValueError
import itertools
import logging

# Data visualization imports
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

# Plotly imports
import dash
from dash import dcc, html, Input, Output
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go

# Data manipulation imports
import pathlib as pl
from ase.visualize import view
from ase import Atoms
from networkx.readwrite import json_graph
import networkx as nx
import pandas as pd
import ase

pio.renderers.default = "browser"
# Todo: we need to cite the original source of the covalent radii data (I am not sure where it was obtained)
elem_cov_radii = {'H': 0.32, 'He': 0.46, 'Li': 1.33, 'Be': 1.02, 'B': 0.85, 'C': 0.75, 'N': 0.71, 'O': 0.63, 'F': 0.64, 'Ne': 0.67, 'Na': 1.55, 'Mg': 1.39, 'Al': 1.26, 'Si': 1.16, 'P': 1.11,
                  'S': 1.03, 'Cl': 0.99, 'Ar': 0.96, 'K': 1.96, 'Ca': 1.71, 'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22, 'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11, 'Ni': 1.1, 'Cu': 1.12, 'Zn': 1.18,
                  'Ga': 1.24, 'Ge': 1.21, 'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17, 'Rb': 2.1, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25,
                  'Pd': 1.2, 'Ag': 1.28, 'Cd': 1.36, 'In': 1.42, 'Sn': 1.4, 'Sb': 1.4, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31, 'Cs': 2.32, 'Ba': 1.96, 'La': 1.8, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74,
                  'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68, 'Gd': 1.69, 'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.7, 'Lu': 1.62, 'Hf': 1.52, 'Ta': 1.46, 'W': 1.37, 'Re': 1.31,
                  'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33, 'Tl': 1.44, 'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01, 'Ac': 1.86, 'Th': 1.75,
                  'Pa': 1.69, 'U': 1.7, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66, 'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65, 'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57, 'Db': 1.49,
                  'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34, 'Mt': 1.29, 'Ds': 1.28, 'Rg': 1.21, 'Cn': 1.22, 'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75, 'Ts': 1.65, 'Og': 1.57}


class HelpFunc:
    def __init__(self):
        """
        This class is a collection of helper functions that are used throughout the assembly process
        """
        pass

    @staticmethod
    def get_and_cast(dictionary: dict, key: str, expected_type: type) -> Any:
        """
        Safely retrieve a value from the dictionary and convert it to the expected type.
        If the expected type is bool, interpret common string/int representations.
        """
        value = dictionary.get(key)

        # Special handling for booleans
        if expected_type is bool:
            if isinstance(value, bool):
                return value
            if isinstance(value, int):
                if value in (0, 1):
                    return bool(value)
            if isinstance(value, str):
                lowered = value.strip().lower()
                if lowered in {"true", "t", "1"}:
                    return True
                if lowered in {"false", "f", "0"}:
                    return False
            raise LoggedValueError(
                f"Fatal Error: '{key}' must be a boolean (true/false, 1/0, t/f)", logger
            )

        # Standard casting fallback
        if not isinstance(value, expected_type):
            try:
                value = expected_type(value)
            except Exception:
                raise LoggedValueError(
                    f"Fatal Error: '{key}' must be of type {expected_type.__name__}", logger
                )
        return value


class LigandSpec:
    REQUIRED_KEYS = {"origin", "ligand_db"}
    OPTIONAL_KEYS = {"n_max": np.inf, "temp_dent": None, "swap_group": None}

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
        self._validate_keys()
        self.origin = self._get_vector(key="origin", required=True)
        self.vectors = self._extract_vectors()
        self.ligand_db_path = data.get("ligand_db", None)
        self.temp_dent = data.get("temp_dent", None)  # TODO: WARNING this is a temporary fix and will need to be removed in the future
        self.n_max = HelpFunc.get_and_cast(self.data, "n_max", int) if "n_max" in self.data else self.OPTIONAL_KEYS["n_max"]
        self.ligand_db = LigandDB.load_from_json(path=self.ligand_db_path, n_max=self.n_max)
        self.swap_group = HelpFunc.get_and_cast(self.data, "swap_group", int) if "swap_group" in self.data else self.OPTIONAL_KEYS["swap_group"]

        # Log the ligand data
        logger.info(f"Ligand '{self.name}' initialized with {len(self.vectors)} vectors and origin {self.origin}.")

    def _validate_keys(self) -> None:
        """
        Checks for missing required keys, duplicate keys, unrecognized keys,
        and ensures at least one 'vector_*' key is present.
        """
        # 1. Check required keys
        all_valid_keys = self.REQUIRED_KEYS.union(self.OPTIONAL_KEYS)
        missing_keys = self.REQUIRED_KEYS - self.data.keys()
        if missing_keys:
            raise LoggedValueError(f"Fatal Error: Missing required keys for ligand '{self.name}': {', '.join(missing_keys)}", logger)

        # 2. Ensure at least one vector_* key exists
        vector_keys = [k for k in self.data if k.startswith("vector")]
        if not vector_keys:
            raise LoggedValueError(f"Fatal Error: Ligand '{self.name}' must have at least one 'vector_*' key (e.g., 'vector_1').", logger)

        # 3. Check for unrecognized keys (excluding vectors)
        non_vector_keys = set(self.data.keys()) - set(vector_keys)
        extra_keys = non_vector_keys - all_valid_keys
        if extra_keys:
            logger.warning(f"Unrecognized keys in ligand '{self.name}': {', '.join(extra_keys)}")

        # 4. Check for duplicate keys (except vector_* keys)
        for key in non_vector_keys:
            count = sum(1 for k in self.data if k == key)
            if count > 1:
                raise LoggedValueError(f"Fatal Error: Ligand '{self.name}' has duplicate key '{key}'. Each non-vector key must occur exactly once.", logger)

    def _extract_vectors(self) -> Dict[str, List[float]]:
        return {key: self._get_vector(key=key, required=False) for key in self.data if key.startswith("vector")}

    def _get_vector(self, key: str, required: bool = False) -> Optional[List[float]]:
        """
        Gets a vector from the ligand data, ensuring it is either:
        - a 3-element list of floats
        - a valid symbolic string ('x', 'y', 'z', '-x', '-y', '-z')

        :param key: The key to retrieve from the ligand data
        :param required: If True, raises an error if the key is not found or is None
        :return: The vector as a list of floats or None if not found
        """
        key_count = sum(1 for k in self.data if k == key)
        if key_count > 1:
            raise LoggedValueError(
                f"Fatal Error: Ligand '{self.name}' must have exactly one '{key}' key, found {key_count}.",
                logger
            )

        value = self.data.get(key)

        if required and value is None:
            raise LoggedValueError(
                f"Fatal Error: '{key}' is required for ligand '{self.name}'.",
                logger
            )

        if value is not None:
            axis_map = {
                "x": [1.0, 0.0, 0.0],
                "y": [0.0, 1.0, 0.0],
                "z": [0.0, 0.0, 1.0],
                "-x": [-1.0, 0.0, 0.0],
                "-y": [0.0, -1.0, 0.0],
                "-z": [0.0, 0.0, -1.0],
            }

            if isinstance(value, str):
                vec = axis_map.get(value.strip().lower())
                if vec is None:
                    raise LoggedValueError(
                        f"Fatal Error: Invalid string '{value}' in '{key}' for ligand '{self.name}'. "
                        f"Expected one of: {', '.join(axis_map.keys())}.",
                        logger
                    )
                return vec

            elif isinstance(value, list) and len(value) == 3:
                if not all(isinstance(v, (int, float)) for v in value):
                    raise LoggedValueError(
                        f"Fatal Error: '{key}' in ligand '{self.name}' must be a 3-element list of numbers.",
                        logger
                    )
                return [float(v) for v in value]

            else:
                raise LoggedValueError(
                    f"Fatal Error: '{key}' in ligand '{self.name}' must be either a symbolic axis string "
                    f"or a 3-element list of numbers.",
                    logger
                )

        return None


class MetalSpec:
    REQUIRED_KEYS = {"metal_type", "metal_oxidation_state", "origin"}
    OPTIONAL_KEYS = {"connectivity": None}

    def __init__(self, name: str, data: Dict[str, Any]):
        """
        Initializes a MetalSpec object. This object is used to store input instructions concerning a metal atom
        :param name: The metal's identifier (e.g., 'metal_1')
        :param data: Dictionary containing metal properties
        :raises ValueError: If required keys are missing
        """

        # validate that the metal data contains the required keys
        self.name = name
        self.data = data
        self.metal_type = self._get_metal_type()
        self.metal_oxidation_state = HelpFunc.get_and_cast(dictionary=self.data,
                                                           key="metal_oxidation_state",
                                                           expected_type=int)
        self.coord = self._get_origin()
        self.connectivity = self.data.get("connectivity", self.OPTIONAL_KEYS["connectivity"])

        # log the metal data
        logger.info(f"Metal '{self.name}' initialized with type '{self.metal_type}' and coordinates {self.coord}.")

    def _validate_keys(self) -> None:
        """
        Ensures all required keys are present and flags unrecognized keys.
        """
        missing = self.REQUIRED_KEYS - self.data.keys()
        if missing:
            raise LoggedValueError(
                f"Fatal Error: Missing required keys for metal '{self.name}': {', '.join(missing)}", logger
            )

        extra = set(self.data.keys()) - self.REQUIRED_KEYS - self.OPTIONAL_KEYS
        if extra:
            logger.warning(f"Unrecognized keys in metal '{self.name}': {', '.join(extra)}")

    def _get_origin(self) -> List[float]:
        coords = HelpFunc.get_and_cast(self.data, "origin", list)
        if len(coords) != 3:
            raise LoggedValueError(f"Fatal Error: 'origin' for metal '{self.name}' must be a 3-element list.", logger)
        return [float(x) for x in coords]

    def _get_metal_type(self) -> str:
        metal_val = self.data.get("metal_type")
        if metal_val is None:
            raise LoggedValueError(f"Fatal Error: 'metal_type' is missing for metal '{self.name}'.", logger)
        metal_str = str(metal_val)
        if metal_str not in PerTab.all_atomic_symbols:
            raise LoggedValueError(f"Fatal Error: Metal '{metal_str}' not found in the periodic table.", logger)
        return metal_str


class BatchInput:
    # Define the required and optional keys for the batch input
    REQUIRED_KEYS = {"name", "random_seed", "max_num_complexes", "total_charge", "geometry", "complex_name_appendix"}
    OPTIONAL_KEYS = {"debug": False,
                     "opt_mono_rot": True,
                     "filter_duplicate_isomers": True,
                     "filter_clashing_structures": True,
                     "filter_clashing_structures_cov_radii_buffer": -0.3,
                     "check_metal_clashes": True,
                     "filter_duplicate_isomers_method": "fingerprint",
                     "filter_duplicate_isomers_grid_size": 9,
                     "isomer_comparison_mode": "max_diff",
                     "isomer_comparison_grouping_mode": "cluster",
                     "isomer_comparison_grouping_cutoff": 1.0,
                     "auxiliary_structure_path": None
                     }

    def __init__(self, batch: Dict[str, Any]):
        """
        Parses and validates a single batch entry from the assembly YAML input.
        Initializes ligand and metal specifications and computes total oxidation state.
        """

        self.batch = batch
        self._validate_batch()

        # Extract and validate the required keys
        self.batch_name = HelpFunc.get_and_cast(dictionary=self.batch,
                                                key="name",
                                                expected_type=str)

        self.random_seed = HelpFunc.get_and_cast(dictionary=self.batch,
                                                 key="random_seed",
                                                 expected_type=int)

        self.max_num_complexes = HelpFunc.get_and_cast(dictionary=self.batch,
                                                       key="max_num_complexes",
                                                       expected_type=int)

        self.total_charge = HelpFunc.get_and_cast(dictionary=self.batch,
                                                  key="total_charge",
                                                  expected_type=int)

        self.geometry = HelpFunc.get_and_cast(dictionary=self.batch,
                                              key="geometry",
                                              expected_type=list)

        self.debug = HelpFunc.get_and_cast(dictionary=self.batch,
                                           key="debug",
                                           expected_type=bool) if "debug" in self.batch else self.OPTIONAL_KEYS["debug"]

        self.opt_mono_rot = HelpFunc.get_and_cast(dictionary=self.batch,
                                                  key="opt_mono_rot",
                                                  expected_type=bool) if "opt_mono_rot" in self.batch else self.OPTIONAL_KEYS["opt_mono_rot"]

        self.filter_duplicate_isomers = HelpFunc.get_and_cast(dictionary=self.batch,
                                                              key="filter_duplicate_isomers",
                                                              expected_type=bool) if "filter_duplicate_isomers" in self.batch else self.OPTIONAL_KEYS["filter_duplicate_isomers"]

        self.filter_clashing_structures = HelpFunc.get_and_cast(dictionary=self.batch,
                                                                key="filter_clashing_structures",
                                                                expected_type=bool) if "filter_clashing_structures" in self.batch else self.OPTIONAL_KEYS["filter_clashing_structures"]

        self.filter_clashing_structures_cov_radii_buffer = HelpFunc.get_and_cast(dictionary=self.batch,
                                                                                 key="filter_clashing_structures_cov_radii_buffer",
                                                                                 expected_type=float) \
            if "filter_clashing_structures_cov_radii_buffer" in self.batch else self.OPTIONAL_KEYS["filter_clashing_structures_cov_radii_buffer"]

        self.check_metal_clashes = HelpFunc.get_and_cast(dictionary=self.batch,
                                                         key="check_metal_clashes",
                                                         expected_type=bool) if "check_metal_clashes" in self.batch else self.OPTIONAL_KEYS["check_metal_clashes"]

        self.filter_duplicate_isomers_method = HelpFunc.get_and_cast(dictionary=self.batch,
                                                                     key="filter_duplicate_isomers_method",
                                                                     expected_type=str) \
            if "filter_duplicate_isomers_method" in self.batch else self.OPTIONAL_KEYS["filter_duplicate_isomers_method"]

        self.filter_duplicate_isomers_grid_size = HelpFunc.get_and_cast(dictionary=self.batch,
                                                                        key="filter_duplicate_isomers_grid_size",
                                                                        expected_type=int) \
            if "filter_duplicate_isomers_grid_size" in self.batch else self.OPTIONAL_KEYS["filter_duplicate_isomers_grid_size"]

        self.isomer_comparison_mode = HelpFunc.get_and_cast(dictionary=self.batch,
                                                            key="isomer_comparison_mode",
                                                            expected_type=str) \
            if "isomer_comparison_mode" in self.batch else self.OPTIONAL_KEYS["isomer_comparison_mode"]
        self.isomer_comparison_grouping_mode = HelpFunc.get_and_cast(dictionary=self.batch,
                                                                     key="isomer_comparison_grouping_mode",
                                                                     expected_type=str) \
            if "isomer_comparison_grouping_mode" in self.batch else self.OPTIONAL_KEYS["isomer_comparison_grouping_mode"]

        self.isomer_comparison_grouping_cutoff = HelpFunc.get_and_cast(dictionary=self.batch,
                                                                       key="isomer_comparison_grouping_cutoff",
                                                                       expected_type=float) \
            if "isomer_comparison_grouping_cutoff" in self.batch else self.OPTIONAL_KEYS["isomer_comparison_grouping_cutoff"]

        self.auxiliary_structure_path = HelpFunc.get_and_cast(dictionary=self.batch,
                                                              key="auxiliary_structure_path",
                                                              expected_type=str) if "auxiliary_structure_path" in self.batch else self.OPTIONAL_KEYS["auxiliary_structure_path"]

        self.complex_name_appendix = HelpFunc.get_and_cast(dictionary=self.batch, key="complex_name_appendix", expected_type=str)

        # Metals and ligands
        self.ligands = []
        self.metals = []
        self._process_geometry()

        # Calculate the total metal oxidation state by summing the oxidation states of all metals
        self.total_metal_oxidation_state = int(sum([metal.metal_oxidation_state for metal in self.metals]))

        # Validate inter-ligand inputs
        self._compare_ligand_inputs()

        # Validate metal connectivity specifications
        self._validate_metal_connectivity_specifications()

        # set logger level based on debug flag
        logger.setLevel(logging.DEBUG if self.debug else logging.INFO)

        # log data
        logger.info(f"Batch '{self.batch_name}' initialized")
        logger.info(f"Provided keys: {', '.join(self.batch.keys())}")
        logger.info(f"Random seed: {self.random_seed}")
        logger.info(f"Max number of complexes: {self.max_num_complexes}")
        logger.info(f"Total charge: {self.total_charge}")
        logger.info(f"Total metal oxidation state: {self.total_metal_oxidation_state}")
        logger.info(f"debug: {self.debug}")
        logger.info(f"opt_mono_rot: {self.opt_mono_rot}")
        logger.info(f"filter_duplicate_isomers: {self.filter_duplicate_isomers}")
        logger.info(f"filter_clashing_structures: {self.filter_clashing_structures}")
        logger.info(f"filter_clashing_structures_cov_radii_buffer: {self.filter_clashing_structures_cov_radii_buffer}")

    def _validate_batch(self) -> None:
        """
        Validates the provided batch dictionary, ensuring all required keys exist
        :raises ValueError: If a required key is missing
        """
        # Ensure the batch is a dictionary
        if not isinstance(self.batch, dict):
            raise LoggedValueError("Fatal Error: Batch input must be a dictionary.", logger)

        # Check for missing and extra keys
        missing_keys = self.REQUIRED_KEYS - set(self.batch.keys())
        extra_keys = set(self.batch.keys()) - self.REQUIRED_KEYS - set(self.OPTIONAL_KEYS)

        if missing_keys:
            raise LoggedValueError(f"Fatal Error: Missing required keys --> input: [{', '.join(missing_keys)}]", logger)

        if extra_keys:
            logger.warning(f"Unrecognized keys --> input: [{', '.join(extra_keys)}]")

    def _validate_metal_connectivity_specifications(self) -> None:
        """
        Validates that the connectivity specifications for each metal are structurally and semantically correct.

        Rules enforced:
        - Connectivity must be a list of dicts.
        - Each referenced species must exist and be unique among ligands and metals.
        - Metals can only connect to ligands via named vectors (e.g., 'vector_1', not 'origin').
        - Metals can only connect to other metals via 'origin'.
        - Metals cannot connect to themselves.

        Raises:
            LoggedValueError if any condition is violated.
        """
        all_species = self.ligands + self.metals

        for metal in self.metals:
            if metal.connectivity is not None:
                if not isinstance(metal.connectivity, list):
                    raise LoggedValueError(
                        f"Connectivity specification for metal '{metal.name}' must be a list of dictionaries.",
                        logger
                    )

                for connection in metal.connectivity:
                    if not isinstance(connection, dict):
                        raise LoggedValueError(
                            f"Each entry in the connectivity list for metal '{metal.name}' must be a dictionary.",
                            logger
                        )

                    species_name = list(connection.keys())[0]
                    targets = list(connection.values())[0]

                    # Prevent self-connection
                    if species_name == metal.name:
                        raise LoggedValueError(
                            f"Fatal Error: Metal '{metal.name}' cannot be connected to itself in the connectivity specification.",
                            logger
                        )

                    # Ensure exactly one matching species exists
                    matching_species = [s for s in all_species if s.name == species_name]
                    if not matching_species:
                        raise LoggedValueError(
                            f"Fatal Error: Unrecognized species '{species_name}' in connectivity for metal '{metal.name}'. "
                            f"Available species: {', '.join([s.name for s in all_species])}",
                            logger
                        )
                    if len(matching_species) > 1:
                        raise LoggedValueError(
                            f"Fatal Error: Duplicate species named '{species_name}' found in connectivity for metal '{metal.name}'.",
                            logger
                        )

                    corresponding = matching_species[0]

                    # Metal–Ligand connections: must use named vectors
                    if corresponding in self.ligands:
                        if not hasattr(corresponding, "vectors") or not isinstance(corresponding.vectors, dict):
                            raise LoggedValueError(
                                f"Ligand '{species_name}' has no valid 'vectors' defined.",
                                logger
                            )
                        for t in targets:
                            if t not in corresponding.vectors:
                                raise LoggedValueError(
                                    f"Fatal Error: Metal '{metal.name}' is attempting to connect to ligand '{species_name}' "
                                    f"using an invalid vector '{t}'. Valid vectors: {list(corresponding.vectors.keys())}",
                                    logger
                                )

                    # Metal–Metal connections: must use origin only
                    elif corresponding in self.metals:
                        for t in targets:
                            if t != "origin":
                                raise LoggedValueError(
                                    f"Fatal Error: Metal '{metal.name}' is attempting to connect to metal '{species_name}' using "
                                    f"'{t}'. Only 'origin' is allowed for metal–metal connections.",
                                    logger
                                )

        logger.info(f"Metal connectivity specifications for batch '{self.batch_name}' validated successfully.")

    def _validate_filter_duplicate_isomers_method(self) -> None:
        """
        Validates the method for filtering duplicate isomers.
        :raises ValueError: If the method is not recognized
        """
        valid_methods = ["fingerprint", "alignment"]
        if self.filter_duplicate_isomers_method not in valid_methods:
            raise LoggedValueError(f"Fatal Error: Unrecognized method for filtering duplicate isomers: {self.filter_duplicate_isomers_method}. "
                                   f"Valid methods are: {', '.join(valid_methods)}", logger)

    def _compare_ligand_inputs(self):
        """
        Compares the ligand inputs to ensure they are consistent across ligands
        :return: None
        """
        if any(ligand.swap_group is not None for ligand in self.ligands):
            if not all(ligand.swap_group is not None for ligand in self.ligands):
                ligands_missing_swap_group = [ligand.name for ligand in self.ligands if ligand.swap_group is None]
                raise LoggedValueError(f"Fatal Error: If a swap_group is specified for any "
                                       f"ligand it MUST be specified for all ligands. swap_group not specified for {ligands_missing_swap_group}", logger)

        return None

    @staticmethod
    def _validate_metal(metal_data: Dict[str, Any]) -> List[float]:
        """
        Validate the metal's coordinate entry
        :param metal_data: Dictionary containing metal position data
        :return: The metal coordinates as a list of floats
        :raises ValueError: If the coordinate format is incorrect
        """
        if "origin" not in metal_data or not isinstance(metal_data["origin"], list) or len(metal_data["origin"]) != 3:
            raise LoggedValueError("Fatal Error: 'metal' must contain a 'origin' key with a 3-element list", logger)
        logger.debug(f"Metal coordinates for '{metal_data.get('metal_type', 'unknown')}' validated: {metal_data['origin']}")
        return [float(x) for x in metal_data["origin"]]

    def _process_geometry(self) -> None:
        """
        Processes the 'geometry' key to extract metal and ligand entries
        """

        for entry in self.geometry:
            if not isinstance(entry, dict) or len(entry) != 1:
                raise LoggedValueError("Fatal Error: 'geometry' must be a list of dictionaries with a single key-value pair", logger)

            key, value = next(iter(entry.items()))

            if key.startswith("metal"):
                self._validate_metal(value)
                self.metals.append(MetalSpec(key, value))
            elif key.startswith("ligand"):
                self.ligands.append(LigandSpec(key, value))
            else:
                raise LoggedValueError(f"Fatal Error: Unexpected key '{key}' in geometry", logger)


class AssemblyComplex(object):

    def __init__(self, ligands: Dict[int, RCA_Ligand], target_vectors: List[Dict[Any, List[float]]], ligand_origins: List[List[float]], metal_origins: List[List[float]],
                 metal_types: List[str], opt_mono_rot: bool = True, filter_duplicate_isomers: bool = True, filter_clashing_structures: bool = True,
                 filter_clashing_structures_cov_radii_buffer: float = -0.3,
                 check_metal_clashes: bool = True,
                 filter_duplicate_isomers_method: str = "fingerprint",
                 filter_duplicate_isomers_grid_size: int = 9,
                 isomer_comparison_mode: str = "max_diff",
                 isomer_comparison_grouping_mode: str = "cluster",
                 isomer_comparison_grouping_cutoff: float = 1.0,
                 swap_groups: Optional[List[int]] = None,
                 debug: bool = False):
        """
        Generates novel transition metal complexes from ligands and metals
        :param ligands:         Dictionary of ligand objects
        :param target_vectors:  List of dictionaries containing target vectors for each ligand
        :param ligand_origins:  List of ligand origins
        :param metal_origins:   List of metal origins
        :param metal_types:     List of metal types

        """
        # Define the class variables
        self.ligands = ligands  # for each ligand, get the ligand object
        self.target_vectors = target_vectors  # List of target vectors i.e [[0, 0, 0], [0, 0, 1], ...]
        self.ligand_origins = ligand_origins  # List of ligand_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_origins = metal_origins  # List of metal_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_types = metal_types  # List of metal types i.e ['Ru', 'Mn', ...]
        self.opt_mono_rot = opt_mono_rot  # Boolean to determine if mono-coordinating ligands should be optimized
        self.filter_duplicate_isomers = filter_duplicate_isomers  # Boolean to determine if duplicate isomers should be filtered
        self.filter_clashing_structures = filter_clashing_structures  # Boolean to determine if clashing structures should be filtered
        self.filter_clashing_structures_cov_radii_buffer = filter_clashing_structures_cov_radii_buffer  # Buffer for clashing structures filtering
        self.check_metal_clashes = check_metal_clashes  # Boolean to determine if metal clashes should be checked
        self.filter_duplicate_isomers_method = filter_duplicate_isomers_method  # Method for filtering duplicate isomers
        self.filter_duplicate_isomers_grid_size = filter_duplicate_isomers_grid_size
        self.isomer_comparison_mode = isomer_comparison_mode
        self.isomer_comparison_grouping_mode = isomer_comparison_grouping_mode
        self.isomer_comparison_grouping_cutoff = isomer_comparison_grouping_cutoff
        self.debug = debug  # Boolean to determine if debug mode is enabled
        self.swap_groups = swap_groups

        # set logger level based on debug flag
        logger.setLevel(logging.DEBUG if self.debug else logging.INFO)

        # Validate the input
        self._validate_input()

    @classmethod
    def from_batch_input(cls, batch_input: BatchInput, ligands):
        """
        Create an AssemblyComplex instance directly from a BatchInput object.
        :param: ligands: the ligands objects to be used in the assembly
        :param: batch_input: BatchInput instance containing parsed input data
        :return: AssemblyComplex instance
        """
        return cls(
            ligands=ligands,
            target_vectors=[ligand.vectors for ligand in batch_input.ligands],
            ligand_origins=[ligand.origin for ligand in batch_input.ligands],
            metal_types=[metal.metal_type for metal in batch_input.metals],
            metal_origins=[metal.coord for metal in batch_input.metals],
            opt_mono_rot=batch_input.opt_mono_rot,
            filter_duplicate_isomers=batch_input.filter_duplicate_isomers,
            filter_clashing_structures=batch_input.filter_clashing_structures,
            filter_clashing_structures_cov_radii_buffer=batch_input.filter_clashing_structures_cov_radii_buffer,
            check_metal_clashes=batch_input.check_metal_clashes,
            filter_duplicate_isomers_method=batch_input.filter_duplicate_isomers_method,
            filter_duplicate_isomers_grid_size=batch_input.filter_duplicate_isomers_grid_size,
            isomer_comparison_mode=batch_input.isomer_comparison_mode,
            isomer_comparison_grouping_mode=batch_input.isomer_comparison_grouping_mode,
            isomer_comparison_grouping_cutoff=batch_input.isomer_comparison_grouping_cutoff,
            swap_groups=[ligand.swap_group for ligand in batch_input.ligands],
            debug=batch_input.debug

        )

    def _validate_input(self) -> None:
        """
        Validates the input to the AssemblyComplex class
        :raises ValueError: If the input is invalid
        """
        if len(self.ligands) != len(self.target_vectors) or len(self.ligands) != len(self.ligand_origins):
            raise LoggedValueError(
                f"Fatal Error: Ligand objects [{len(self.ligands)}], target vectors [{len(self.target_vectors)}], and ligand origins [{len(self.ligand_origins)}] must have the same length",
                logger)

        if len(self.metal_origins) != len(self.metal_types):
            raise LoggedValueError(
                f"Fatal Error: Metal origins [length: {self.metal_origins}] and metal types [length: {self.metal_types}] must have the same length for multi-metallic systems.",
                logger)

    def _validate_swap_groups(self) -> None:
        """
        Validates the swap groups for ligands. Ligands of the same 'swap_group' must have the same effective ligand
        coordination number (elcn) to be considered for swapping. Will raise an error if this is not the case.
        """

        if any(swap_tag is None for swap_tag in self.swap_groups) and not all(swap_tag is None for swap_tag in self.swap_groups):
            raise LoggedValueError("Fatal Error: If a swap_group is specified for any ligand, it MUST be specified for all ligands.", logger)

        elif all(swap_tag is None for swap_tag in self.swap_groups):
            logger.warning("No swap groups specified. Using effective ligand coordination numbers (elcn) as swap groups.")
            self.swap_groups = [ligand.elcn for ligand in self.ligands]

        else:
            # All swap_groups are present; now check if ligands that are specified by the user to be swapped have the
            # same effective ligand coordination number (elcn) and make sense to be swapped
            list1 = self.swap_groups
            list2 = [ligand.elcn for ligand in self.ligands]
            groups = defaultdict(list)

            for i, v in enumerate(list1):
                groups[v].append(i)

            for group_id, indices in groups.items():
                elcns = {list2[i] for i in indices}
                if len(elcns) > 1:
                    details = ", ".join(f"{self.ligands[i].name} (elcn={list2[i]})" for i in indices)
                    raise LoggedValueError(
                        f"Fatal Error: Ligands in swap_group {group_id} must all have the same elcn, but found mismatches: {details}",
                        logger
                    )
        logger.debug("Swap groups validated successfully.")
        return None

    def generate(self) -> List[Atoms]:
        # Generate all possible isomers by exchanging compatible ligands
        ligand_lists = self._assign_ligands_to_vectors(
            ligands=list(self.ligands),
            vectors=self.target_vectors,
            swap_groups=self.swap_groups
        )

        # Loop through each ligand set permutation and generate the each isomer
        all_isomers = []
        all_DARTLigands = []  # List to store the DARTLigand objects used in each isomer
        for ligand_swapped_combo in ligand_lists:
            aligned_ligands = []  # list of a complexes ligands i.e. [l1, l2, [l3-i1, l3-i2], l4, ...] where i1, i2, etc. correspond to a unique coordination of a ligand within a binding site
            aligned_DARTLigands = []  # list of DARTLigand objects used in the complex
            for idx, (ligand, ligand_target_vectors, origin) in enumerate(zip(ligand_swapped_combo, self.target_vectors, self.ligand_origins)):
                # Retrieve the ligand's geometry (ASE atoms object) and donor atom indices (List[List[index]).
                geometry, donor_atoms = ligand.get_isomers_effective_ligand_atoms_with_effective_donor_indices()

                # Tag each atom in the ligand geometry with a 2D label:
                #   [0] = ligand index in the complex (elcn, e.g., for filtering/grouping)
                #   [1] = position of the ligand in the coordination sequence (1-based index)
                #   [2] = 1 if atom is a donor atom else zero
                # This enables downstream tracking of ligand identity and spatial assignment during isomer generation and filtering.
                donor_idc = list(self._flatten(ligand.get_denticities_and_hapticities_idc()))
                n_atoms = len(geometry)
                tags = np.full((n_atoms, 3), [ligand.n_eff_denticities, idx + 1, 0], dtype=int)
                tags[donor_idc, 2] = 1

                geometry.new_array('multi_tags', tags)

                # Convert target vector dictionary values to numpy arrays.
                vector_labels = list(ligand_target_vectors.keys())
                target_vectors = [np.array(v) for v in ligand_target_vectors.values()]

                # Align the donor atoms of the ligand to the target vectors
                # TODO this code may be replaced in the future, I think Timo might have a more robust version of this code
                ligand_isomers, donor_atoms_ordered, rssd = try_all_geometrical_isomer_possibilities(atoms=geometry,
                                                                                                     donor_idc=donor_atoms[0],
                                                                                                     vector_labels=vector_labels,
                                                                                                     target_vectors=target_vectors)

                # Translate the ligand to its correct location in the complex
                for isomer in ligand_isomers:
                    isomer.set_positions(isomer.get_positions() + np.array(origin))  # Note: I believe this method assumes that the ligand has been pre-translated to 0,0,0

                # Remove dummy atoms (e.g., "Cu") from haptic ligands.
                cleaned_isomers = self._remove_haptic_dummy_atom(
                    atoms_list=ligand_isomers,
                    dummy_atom="Cu",
                    donor_atoms_idc=ligand.hapdent_idc
                )
                # Append the rotated ligands to the list
                aligned_ligands.append(cleaned_isomers)
                aligned_DARTLigands.append([ligand for _ in range(len(cleaned_isomers))])

            all_isomers.append(aligned_ligands)
            all_DARTLigands.append(aligned_DARTLigands)

        # Assemble the complexes using the ligand permutations
        all_isomers, ligands_used = self._gen_all_isomers(ligands=all_isomers, DARTLigands=all_DARTLigands)

        # Optimize mono-coordinating ligands if requested
        if self.opt_mono_rot:
            logger.info("Optimizing rotation of mono-coordinating ligands...")
            opt_mono_rot = AxialOpt(complexes=all_isomers,
                                    target_vectors=self.target_vectors,
                                    ligand_origins=self.ligand_origins)
            all_isomers = opt_mono_rot.opt_mono_rotation()

        # Remove duplicate isomers if requested
        if self.filter_duplicate_isomers:
            logger.info("Filtering duplicate isomers...")
            reduce_isom = ReduceIsomersV2(isomers=all_isomers,
                                          method=self.filter_duplicate_isomers_method,
                                          grid_size=self.filter_duplicate_isomers_grid_size,
                                          metal_centres=self.metal_origins,
                                          isomer_comparison_mode=self.isomer_comparison_mode,
                                          isomer_comparison_grouping_mode=self.isomer_comparison_grouping_mode,
                                          fingerprint_grouping_cutoff=self.isomer_comparison_grouping_cutoff)
            all_isomers = reduce_isom.reduce_isomers()
            ligands_used = [ligands_used[i] for i in reduce_isom.output_indices]
            # Note the code can be paused here and you can inspect the isomers with reduce_isom.plot_fingerprint_difference_matrix()

        # Filter clashing structures if requested
        if self.filter_clashing_structures:
            logger.info("Filtering clashing structures...")
            clashing = ClashFilter(isomers=all_isomers,
                                   covalent_radii=elem_cov_radii,
                                   buffer=self.filter_clashing_structures_cov_radii_buffer,
                                   check_metal_clashes=self.check_metal_clashes)
            all_isomers, retained_indices = clashing.process_isomers()
            ligands_used = [ligands_used[i] for i in retained_indices]

        return all_isomers, ligands_used

    def _flatten(self, nested):
        """
        Recursively flattens a nested list or tuple into a flat generator of elements.
        :param: nested (list or tuple): Arbitrarily nested list/tuple of elements (e.g., indices)
        :return: Individual flattened elements
        """
        for item in nested:
            if isinstance(item, (list, tuple)):
                yield from self._flatten(item)
            else:
                yield item

    @staticmethod
    def _assign_ligands_to_vectors(ligands: List[RCA_Ligand],
                                   vectors: List[Dict[str, List[float]]],
                                   swap_groups: List[int]) -> List[List[RCA_Ligand]]:
        """
        Assigns ligands to vector entries based on swap group IDs.
        Ligands in the same swap group can be swapped among vector entries assigned the same swap group ID.

        :param ligands:       A list of ligand objects.
        :param vectors:       A list of dictionaries, each representing a coordination site with arbitrary keys.
        :param swap_groups:   A list of integers where the i-th value defines the swap group for the i-th ligand.
                              Ligands with the same group ID can be swapped; those with different IDs cannot.
                              This list must match the length and order of the `vectors` list.
        :return:              A list of ligand groupings for each vector entry.
                              Each entry contains the swappable ligands for that vector,
                              used later for permutation generation.
        :raise:               LoggedValueError if swap_groups are inconsistent with ligand or vector count.
        """

        # Validate swap_groups length
        if len(vectors) != len(swap_groups):
            raise LoggedValueError("Fatal Error: swap_groups must match the number of target vectors.", logger)
        if len(ligands) != len(swap_groups):
            raise LoggedValueError("Fatal Error: swap_groups must match the number of ligands.", logger)

        # Step 1: Build a mapping of swap_group ID → list of ligands that belong to that group.
        group_to_ligands = defaultdict(list)
        for ligand, group_id in zip(ligands, swap_groups):
            group_to_ligands[group_id].append(ligand)

        # Step 2: Build the vector-to-ligand group mapping.
        # For each vector entry, assign the appropriate ligand group based on the vector’s corresponding swap group ID.
        _list = [group_to_ligands[group_id] for group_id in swap_groups]

        # Step 3: Generate all valid permutations of ligand assignments.
        # For each vector position, choose one ligand from the allowed group — ensuring no ligand is reused.
        results = []
        current_permutation = []
        used = set()

        def backtrack(index: int):
            # Base case: if all vector positions have been filled, record the permutation.
            if index == len(_list):
                results.append(current_permutation.copy())
                return

            # Iterate over allowed ligands for the current position.
            for ligand in _list[index]:
                if ligand in used:
                    continue  # Skip ligands already assigned in this permutation.
                used.add(ligand)
                current_permutation.append(ligand)
                backtrack(index + 1)
                # Backtrack: remove the ligand from current permutation and reuse it elsewhere.
                current_permutation.pop()
                used.remove(ligand)

        backtrack(0)
        return results

    def _add_metals(self, ligand_structure: Atoms) -> Atoms:
        """
        Adds the metals to the complex, ensuring they are first in the Atoms object.
        :param ligand_structure: Atoms object representing the ligand structure
        :return: Atoms object with metals prepended
        """
        # Start with an empty Atoms object to accumulate metals
        metal_atoms = Atoms()

        for idx, (metal_type, metal_origin) in enumerate(zip(self.metal_types, self.metal_origins)):
            # Create a single-atom ASE structure for the metal
            metal = Atoms(symbols=metal_type, positions=[metal_origin])
            metal.new_array('multi_tags', np.full((1, 3), [0, 0, 0], dtype=int))  # Only one atom
            metal_atoms += metal
            logger.debug(f"Metal [{idx + 1}/{len(self.metal_types)}]: [{metal_type}] added at origin [{metal_origin}]")

        # Combine metals first, then ligand structure
        combined = metal_atoms + ligand_structure
        return combined

    def _gen_all_isomers(
            self,
            ligands: List[List[List[Atoms]]],
            DARTLigands: List[List[List[RCA_Ligand]]]
    ) -> Tuple[List[Atoms], List[List[RCA_Ligand]]]:
        """
        Generate all possible isomers from a list of ligands which have multiple isomers.

        :param ligands: Nested list of ASE ligand structures.
        :param DARTLigands: Nested list of corresponding RCA_Ligand objects.
        :return: Tuple of (assembled ASE Atoms objects, per-isomer list of RCA_Ligand objects used).
        """
        isomers = []
        ligand_compositions = []

        logger.info(f"Generating isomers from {len(ligands)} ligand slots.")

        for ligand_atoms_list, ligand_obj_list in zip(ligands, DARTLigands):
            combinations_atoms = list(itertools.product(*ligand_atoms_list))
            combinations_ligands = list(itertools.product(*ligand_obj_list))

            for idx, (combo_atoms, combo_ligands) in enumerate(zip(combinations_atoms, combinations_ligands)):
                combined = Atoms()
                for ligand_atoms in combo_atoms:
                    combined += ligand_atoms

                combined = self._add_metals(ligand_structure=combined)
                isomers.append(combined)
                ligand_compositions.append(list(combo_ligands))

                logger.debug(f"Isomer {idx + 1}/{len(combinations_atoms)} assembled.")

        logger.info(f"Total number of assembled isomers: {len(isomers)}.")
        return isomers, ligand_compositions

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


class ClashFilter:
    """
    Filters out isomers that have atomic clashes based on covalent radii.
    Optionally considers ligand-metal clashes separately.
    """

    def __init__(
            self,
            isomers: List[Atoms],
            covalent_radii: dict,
            buffer: float = -0.3,
            check_metal_clashes: bool = True
    ):
        """
        :param isomers: List of ASE Atoms objects.
        :param covalent_radii: Dictionary mapping atomic symbols to covalent radii.
        :param buffer: Float to tune clash threshold (default -0.3 Å).
        :param check_metal_clashes: Whether to consider metal-ligand clashes (default True).
        """
        self.isomers = isomers
        self.cov_radii = covalent_radii
        self.buffer = buffer
        self.check_metal_clashes = check_metal_clashes
        self.filtered_isomers = []
        self.rejected_isomers = []
        logger.debug(f"Initialized ClashFilter with {len(self.isomers)} isomers, "
                     f"covalent radii dict, "
                     f"buffer: {self.buffer}, "
                     f"check_metal_clashes: {self.check_metal_clashes}")

    def process_isomers(self) -> Tuple[List[Atoms], List[int]]:
        """
        Filters out isomers that have atomic clashes.
        :return: Tuple of (filtered isomers, indices of retained isomers in original list)
        """
        retained_indices = []
        for idx, isomer in enumerate(self.isomers):
            if not self.has_clash(isomer):
                self.filtered_isomers.append(isomer)
                retained_indices.append(idx)
            else:
                self.rejected_isomers.append(isomer)
                logger.debug(f"Clash detected in isomer [{idx}] of [{len(self.isomers)}], skipping.")

        logger.info(f"Filtered out {len(self.isomers) - len(self.filtered_isomers)} isomers due to atomic clashes.")
        return self.filtered_isomers, retained_indices

    def has_clash(self, atoms: Atoms) -> bool:
        """
        Determines whether any atomic clashes exist in the structure.
        Atom pairs are checked based on their multi_tags.

        :param atoms: ASE Atoms object with a 'multi_tags' array.
        :return: True if a clash is found, False otherwise.
        """
        if 'multi_tags' not in atoms.arrays:
            raise LoggedValueError(
                "Fatal Error: 'multi_tags' array not found in atoms object. "
                "Ensure that the atoms object has been properly tagged.",
                logger
            )

        tags = atoms.get_array('multi_tags')  # shape: (n_atoms, 2)
        n_atoms = len(atoms)
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                tag_i = tags[i]
                tag_j = tags[j]

                # Skip same-ligand checks
                if tag_i[1] == tag_j[1]:
                    continue

                # Skip ligand-metal checks if disabled
                if not self.check_metal_clashes and (tag_i[1] == 0 or tag_j[1] == 0):
                    continue

                dist = np.linalg.norm(positions[i] - positions[j])
                r1 = self.cov_radii.get(symbols[i], 1.0)
                r2 = self.cov_radii.get(symbols[j], 1.0)

                if dist < (r1 + r2 + self.buffer):
                    logger.debug(
                        f"Clash detected between {symbols[i]}{i} (tag {tag_i}) and "
                        f"{symbols[j]}{j} (tag {tag_j}) at distance {dist:.2f} Å."
                    )
                    return True

        return False

    def view_rejected_isomers(self):
        """
        View the rejected isomers in ASE GUI
        """
        if not self.rejected_isomers:
            logger.info("No rejected isomers to view.")
            return

        view(self.rejected_isomers)


class AxialOpt:
    def __init__(self, complexes: List[Atoms], target_vectors: List[Dict[str, List[float]]], ligand_origins: List[List[float]]):
        """
        This class will take a list of ASE Atoms objects and optimize mono-coordinating ligands around their coordination axis
        """
        self.input_complexes = complexes
        self.target_vectors = target_vectors
        self.ligand_origins = ligand_origins
        self.output_complexes = []

        # log the initialization
        logger.info(f"AxialOpt initialized with {len(self.input_complexes)} complexes, "
                    f"{len(self.target_vectors)} target vectors, "
                    f"and {len(self.ligand_origins)} ligand origins.")

    def opt_mono_rotation(self):
        """
        Optimize the rotation of each mono-coordinating ligand around their respective coordination axis simultaneously
        """
        # Set a random seed for reproducibility
        np.random.seed(42)

        # Specify bounds for the optimization: each angle can vary from 0 to 360 degrees.
        bounds = [[0, 360] for _ in self.target_vectors]  # Bounds is set for each ligand regardless of its coordination number

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
                logger.debug(f"Rotated ligand with tag {tag} by {angle:.2f} degrees around vector {axis} at origin {origin}.")

            self.output_complexes.append(tmc)

        logger.info(f"Optimized {len(self.output_complexes)} complexes with mono-coordinating ligand rotations.")
        return self.output_complexes

    def visualize_structures(self):
        """
        Visualize input and output complexes interleaved using ASE's GUI.
        Each input is followed by its corresponding optimized output.
        """
        if not self.output_complexes:
            logger.warning("No output complexes found. Run opt_mono_rotation() first.")
            return

        structures_to_view = []

        for i, (input_tmc, output_tmc) in enumerate(zip(self.input_complexes, self.output_complexes)):
            input_copy = input_tmc.copy()
            input_copy.info["label"] = f"Input {i}"

            output_copy = output_tmc.copy()
            output_copy.info["label"] = f"Optimized {i}"

            structures_to_view.extend([input_copy, output_copy])  # interleave input/output

        logger.info(f"Launching viewer for {len(structures_to_view)} structures...")
        view(structures_to_view)

    def objective_function(self, x: np.ndarray, vectors_in: List[np.array], origins_in: List[np.array], TMC_in: Atoms, ):
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

        # Retrieve the multi_tags array
        multi_tags = TMC_worker.get_array("multi_tags")

        # Get the unique set of all nonzero second tags (each unique tag represents a unique ligand, a zero tag represents a metal).
        unique_tags = [tag for tag in np.unique(multi_tags[:, 1]) if tag != 0]

        # Loop through each unique tag and apply the rotation
        for tag, angle, axis, origin in zip(unique_tags, list(x), vectors_in, origins_in):

            # Get indices where the second tag in the list is equal to the current "tag" (essentially the indices of the atoms in this particular ligand)
            indices = np.where(multi_tags[:, 1] == tag)[0]

            # Check if any of these indices have a first tag not equal to 1, if so, skip
            # This ensures only ligands that have a effective coordination number of 1 are rotated
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


class ReduceIsomersV2:
    """
    Class to reduce the number of isomers based on alignment or fingerprint similarity.
    """

    def __init__(self, isomers: List[Atoms],
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
        self.output_indices = []
        self.similarity_cutoff_used = None

    def reduce_isomers(self) -> List[Atoms]:
        """
        Reduce the number of isomers based on the specified method.
        :return: unique isomers as a list of ASE Atoms objects
        """

        if self.method == "alignment":
            self.output_isomers = self._reduce_by_alignment()
        elif self.method == "fingerprint":
            self.output_isomers = self._reduce_by_fingerprint()
        else:
            raise LoggedValueError(f"Fatal Error: Unsupported reduction method '{self.method}'. "
                                   f"Supported methods are 'alignment' and 'fingerprint'.", logger)
        logger.info(f"Reduced isomers from {len(self.isomers)} to {len(self.output_isomers)} using method '{self.method}'.")
        return self.output_isomers

    def _reduce_by_alignment(self) -> List[Atoms]:
        """
        Reduce isomers by aligning them and calculating RMSD or another distance metric.
        :return:
        """

        n = len(self.isomers)
        diff_matrix = np.zeros((n, n))
        counter = 0
        total = n * (n - 1) // 2  # Total comparisons in upper triangle
        for i in range(n):
            isomer1 = self.isomers[i]
            for j in range(i + 1, n):  # only compute upper triangle
                isomer2 = self.isomers[j]
                logger.debug(f"Comparing isomers [{i}] and [{j}] ... [{counter}/{total}]")
                score = self.align_isomers(isomer1, isomer2)
                diff_matrix[i, j] = score
                diff_matrix[j, i] = score  # mirror to lower triangle
                counter += 1

        self.diff_matrix = diff_matrix

        # Cluster labels: "Close" or "Far" for each pair
        group_labels_matrix = self._analyze_similarity(self.diff_matrix, quantile=0.2,
                                                       method=self.isomer_comparison_grouping_mode,
                                                       cutoff=self.fingerprint_grouping_cutoff)

        visited = set()
        unique_indices = []

        for i in range(n):
            if i in visited:
                continue
            unique_indices.append(i)
            # Mark all j similar to i as visited
            for j in range(n):
                if group_labels_matrix[i, j] == "Close":
                    visited.add(j)

        # Use selected indices to form output isomers
        self.output_indices = unique_indices
        self.output_isomers = [self.isomers[i] for i in unique_indices]
        logger.info(f"Fingerprint reduction retained [{len(self.output_isomers)}] unique isomers of a total [{n}] isomers.")

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
        assert hasattr(self, "metal_centres"), LoggedValueError("Fatal Error: metal_centres must be defined before calling align_isomers. ")
        logger.debug(f"Aligning isomers based on {len(self.metal_centres)} metal centre(s).")

        # Here the number of metal centres and the fact that metal centres of different isomers
        # must always be aligned is taken into account to determine if the isomers are similar or not
        if len(self.metal_centres) == 1:
            # There are 3 axes which an isomer can be rotated around to align it with another isomer
            bounds = [[0, 360] for _ in range(3)]
            # The 3 cardinal axes, properly centered on the metal centre
            axes = np.eye(3)  # Standard Cartesian axes (x, y, z)
            logger.debug("Performing 3D brute-force alignment over x, y, z axes.")

        elif len(self.metal_centres) == 2:
            # The isomer can only be rotated around the metal-metal axis to determine if the isomers are similar or not
            bounds = [[0, 360] for _ in range(1)]
            axis_vector = np.array(self.metal_centres[1]) - np.array(self.metal_centres[0])
            axis_vector /= np.linalg.norm(axis_vector)  # Normalize
            axes = [axis_vector]
            logger.debug(f"Performing 1D brute-force alignment around axis: {axis_vector.tolist()}")

        elif len(self.metal_centres) >= 3:
            # 3 metal centres means each isomer is fixed in space and their geometries can be directly compared
            # We simply return the energy heuristic
            logger.debug("Three or more metal centres detected — skipping alignment and using direct heuristic.")
            return self.energy_heuristic(stat_atoms=stationary_atoms, rot_atoms=rotated_atoms)

        else:
            raise LoggedValueError(f"Fatal Error: Unsupported number of metal centres ({len(self.metal_centres)}). ", logger)

        # Perform brute-force global optimization using the configured grid density
        result_angles = brute(
            func=self.objective_function,
            ranges=bounds,
            args=(stationary_atoms, rotated_atoms, axes),
            Ns=self.grid_size
        )

        min_score = self.objective_function(np.array(result_angles), stationary_atoms, rotated_atoms, axes)
        logger.debug(f"Best alignment score after brute-force search: {min_score:.4f}")
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

    def _reduce_by_fingerprint(self) -> List[Atoms]:
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

        n = len(self.isomers)
        visited = set()
        unique_indices = []

        for i in range(n):
            if i in visited:
                continue
            unique_indices.append(i)
            # Mark all j similar to i as visited
            for j in range(n):
                if group_labels_matrix[i, j] == "Close":
                    visited.add(j)

        # Use selected indices to form output isomers
        self.output_isomers = [self.isomers[i] for i in unique_indices]
        self.output_indices = unique_indices
        logger.info(f"Fingerprint reduction retained [{len(self.output_isomers)}] unique isomers of a total [{n}] isomers.")

        return self.output_isomers

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
            # Flatten the matrix values for clustering.
            # Upper triangle indices are used to avoid redundancy.
            triu_indices = np.triu_indices_from(matrix, k=1)
            values = matrix[triu_indices].reshape(-1, 1)
            bandwidth = estimate_bandwidth(values, quantile=quantile, n_samples=len(values))
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
                raise LoggedValueError("Fatal Error: Cutoff value must be provided when using 'cutoff' method.", logger)
            group_labels = np.where(matrix <= cutoff, "Close", "Far")
            self.similarity_cutoff_used = cutoff
            return group_labels

        else:
            raise LoggedValueError(
                f"Fatal Error: Unsupported similarity analysis method '{method}'. "
                f"Supported options: 'cluster', 'cutoff'.", logger)

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

        # For a consistent fingerprint, the distances for each element pair is sorted
        # and then concatenate them in a fixed order (e.g., lexicographical order of the keys)
        fingerprint_parts = [
            d for key in sorted(pair_distances)
            for d in sorted(pair_distances[key])
        ]
        return np.array(fingerprint_parts)

    @staticmethod
    def _fingerprint_comparison(fp1: np.ndarray, fp2: np.ndarray, mode: str = "max_diff"):
        """
        Compute the maximum absolute difference between two sorted fingerprint vectors.
        :param fp1: 1D numpy array fingerprint for isomer 1.
        :param fp2: 1D numpy array fingerprint for isomer 2.
        :param mode: Comparison mode: 'max_diff', 'sum_diff', 'mean_diff', or 'rmsd'.
        :return: Maximum absolute difference.
        """
        logger.debug(f"Comparing fingerprints with mode '{mode}'.")
        if mode == "max_diff":
            return np.max(np.abs(fp1 - fp2))
        elif mode == "sum_diff":
            return np.sum(np.abs(fp1 - fp2))
        elif mode == "mean_diff":
            return np.mean(np.abs(fp1 - fp2))
        elif mode == "rmsd":
            return np.sqrt(np.mean((fp1 - fp2) ** 2))
        else:
            raise LoggedValueError(
                f"Fatal Error: Unsupported fingerprint comparison mode '{mode}'. "
                f"Supported modes are 'max_diff', 'sum_diff', 'mean_diff', and 'rmsd'.", logger)

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

        assert cell_label_mode in {"value", "group", "none"}, LoggedValueError(
            f"Fatal Error: Unsupported cell label mode '{cell_label_mode}'. "
            "Supported modes are 'value', 'group', and 'none'.", logger)

        # Fallback defaults for colors
        min_color = {"r": 238, "g": 100, "b": 97} if min_color is None else min_color
        mid_color = {"r": 255, "g": 255, "b": 255} if mid_color is None else mid_color
        max_color = {"r": 12, "g": 171, "b": 185} if max_color is None else max_color

        if self.diff_matrix is None:
            self.reduce_isomers()

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

        app.layout = html.Div([
            dcc.Graph(id="heatmap", figure=fig),
            dcc.Graph(id="histogram", figure=fig_hist),

        ])

        @app.callback(
            Output("heatmap", "figure"),
            Input("heatmap", "clickData")
        )
        def display_click(clickData):
            if clickData:
                point = clickData["points"][0]
                i_idx = int(point["y"])
                j_idx = int(point["x"])
                logger.info(f"Clicked cell: ({i_idx}, {j_idx}) — triggering isomer viewer.")
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

        isomer1 = self.isomers[index1].copy()
        isomer2 = self.isomers[index2].copy()
        isomer2_aligned = self.isomers[index2].copy()

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
            logger.warning("Three or more metal centres — skipping rotation; showing structures unaligned.")
        else:
            raise LoggedValueError("Fatal Error: Invalid number of metal centres for alignment.", logger)

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


class AtomsCombiner:
    def __init__(self, base_atoms: Atoms, xyz_path: Union[str, pl.Path]):
        """
        Initialize with a base Atoms object and a path to an .xyz file.

        :param base_atoms: ASE Atoms object
        :param xyz_path: Path to .xyz file (str or pathlib.Path)
        """
        if not isinstance(base_atoms, Atoms):
            raise LoggedValueError("Fatal Error: base_atoms must be an ASE Atoms object.", logger)

        self.base_atoms = base_atoms
        self.xyz_path = pl.Path(xyz_path)
        self.xyz_atoms = self._load_xyz()

    @classmethod
    def from_batch_input(cls, batch_input: BatchInput, base_atoms: Atoms):
        """
        Create an instance from a BatchInput object and a base Atoms object.
        :param batch_input: BatchInput object containing auxiliary_structure_path
        :param base_atoms: ASE Atoms object to combine with the .xyz file
        """
        return cls(
            base_atoms=base_atoms,
            xyz_path=batch_input.auxiliary_structure_path
        )

    def _load_xyz(self) -> Atoms:
        """Read the .xyz file and return an ASE Atoms object."""
        if not self.xyz_path.exists():
            return None
        atoms = ase.io.read(str(self.xyz_path))  # ASE does not accept Path directly
        if not isinstance(atoms, Atoms):
            raise LoggedValueError(f"Fatal Error: The file {self.xyz_path} does not contain valid ASE Atoms data.", logger)
        return atoms

    def combine(self) -> Atoms:
        """Return a new Atoms object with the base and xyz Atoms combined."""
        if self.xyz_atoms is None:
            logger.warning("No auxiliary structure found in .xyz file; returning base atoms only.")
            return self.base_atoms
        return self.base_atoms + self.xyz_atoms


class DARTIsomer(RCA_Molecule):
    """
    generates a RCA_Molecule object from an ASE Atoms object and
    a list of ligands
    """

    def __init__(self,
                 atoms: Atoms,
                 ligands: List[RCA_Ligand],
                 default_graph: bool = True,
                 ligand_target_vectors: List[List[float]] = None,
                 ligand_origins: List[List[float]] = None,
                 connectivity: Dict = None,
                 metal_centers: Union[List[List[Union[str, List[float]]]], str] = None):

        self._ani = None
        self.atoms = atoms
        self.ligands = ligands
        self.default_graph = default_graph
        self.ligand_target_vectors = ligand_target_vectors
        self.ligand_origins = ligand_origins
        self.metal_centers = metal_centers
        self.ligand_info = self._get_ligand_info()
        self.connectivity = connectivity

        if not len(self.ligands) == len(self.ligand_target_vectors) == len(self.ligand_origins):
            raise LoggedValueError(f"Fatal Error: The input of ligands, target vectors, ligand origins and metal centers must have the same length."
                                   f"Respective lenghts given: [ligands: {len(self.ligands)}, target_vectors: {len(self.ligand_target_vectors)}, "
                                   f"ligand_origins: {len(self.ligand_origins)}.", logger)

        assert [len(vector_set) for vector_set in self.ligand_target_vectors] == [ligand.n_eff_denticities for ligand in self.ligands], \
            LoggedValueError("Fatal Error: The number of target vectors must match the number of ligand donor atoms.", logger)

        super().__init__(
            atomic_props=self.atoms,
            graph=self._get_graph(),
            global_props=None)

    @classmethod
    def from_batch_input(cls, atoms: Atoms, ligands: List[RCA_Ligand], batch_input: BatchInput, default_graph: bool = True) -> "DARTIsomer":
        """
        Create a DARTIsomer instance from a BatchInput object.
        :param atoms:
        :param ligands:
        :param batch_input: BatchInput object containing auxiliary_structure_path, ligands, target_vectors, and ligand_origins.
        :return: DARTIsomer instance
        """
        logger.debug(f"Creating DARTIsomer from BatchInput with {len(ligands)} ligands and {len(batch_input.metals)} metal centers.")
        target_vectors = [ligand.vectors for ligand in batch_input.ligands]
        ligand_origins = [ligand.origin for ligand in batch_input.ligands]
        metal_centers = [metal for metal in batch_input.metals]

        return cls(
            atoms=atoms,
            default_graph=default_graph,
            ligands=ligands,
            ligand_target_vectors=target_vectors,
            ligand_origins=ligand_origins,
            metal_centers=metal_centers,
            connectivity=[metal.connectivity for metal in batch_input.metals]
        )

    def _get_ligand_info(self) -> Dict[str, Any]:
        """
        Extract information about the ligands in this isomer.
        - unique_names: A list of unique ligand names that appears in the same
        order as they appear in the ASE Atoms object.
        :return: Dictionary containing:
                - unique_names: List of unique ligand names.
                - geometries: List of ligand geometries.
                - donors: List of donor atoms for each ligand, represented as strings.
                - charges: List of predicted charges for each ligand.
                - stoichiometries: List of ligand stoichiometries.
        """

        return {'unique_names': list(set(ligand.unique_name for ligand in self.ligands)),
                'geometries': [ligand.geometry for ligand in self.ligands],
                'donors': ['-'.join(sorted(ligand.local_elements)) for ligand in self.ligands],
                'charges': [ligand.pred_charge for ligand in self.ligands],
                'stoichiometries': [lig.stoichiometry for lig in self.ligands]}

    def _get_graph(self):
        """
        Generate a graph where (default):
            - Every metal is connected to every other metal.
            - Every metal is connected to every ligand donor atom.
            - Ligand donor atoms are not connected to each other.
        :return: networkx.Graph
        """

        graph = nx.Graph()
        tags = self.atoms.get_array("multi_tags")
        vector_keys = self.atoms.get_array("vector_keys")
        metal_idc = self._get_metal_idc()
        donor_idc = self._get_donor_idc()

        if self.default_graph and all(connectivity is None for connectivity in self.connectivity):
            # --- Step 1: Add all atoms as graph nodes ---
            for i, atom in enumerate(self.atoms):
                graph.add_node(i, node_label=atom.symbol)

            # --- Step 3: Connect all metal atoms to each other ---
            for i in metal_idc:
                for j in metal_idc:
                    if i < j:
                        graph.add_edge(i, j)

            # --- Step 4: Connect each metal to all donor atoms ---
            for metal_idx in metal_idc:
                for donor_idx in donor_idc:
                    graph.add_edge(metal_idx, donor_idx)

        elif all(connectivity is not None for connectivity in self.connectivity):
            # --- Guided connectivity via user-defined input in self.connectivity ---
            # This block builds the graph based on explicit user-provided metal–ligand and metal–metal connections.

            for metal_idx, metal_conn in zip(metal_idc, self.connectivity):
                # metal_idx: global atom index for the current metal atom
                # metal_conn: list of connection dicts for this metal, e.g., [{'ligand_1': ['vector_1']}, ...]

                for conn in metal_conn:
                    # Each conn is a dictionary defining connections to another metal or ligand
                    # e.g., {'ligand_3': ['vector_1', 'vector_2']}

                    for partner_key, vectors in conn.items():
                        # partner_key: string identifier of the connection partner ('ligand_X' or 'metal_Y')
                        # vectors: list of vector keys (e.g., ['vector_1']) to match against vector_keys array

                        if partner_key.startswith("ligand_"):
                            # --- Handle metal–ligand connections ---
                            # Extract ligand number from string key (1-based indexing from user)
                            ligand_num = int(partner_key.split("_")[1]) - 1

                            # ligand_tag corresponds to tags[:, 1] in multi_tags, which starts from 1
                            ligand_tag = ligand_num + 1

                            # Construct boolean mask:
                            # - (tags[:, 1] == ligand_tag): atoms belonging to the specific ligand
                            # - np.isin(vector_keys, vectors): atoms whose assigned vector_key matches one of the specified vectors
                            atom_mask = (tags[:, 1] == ligand_tag) & np.isin(vector_keys, vectors)

                            # Extract atom indices matching both ligand identity and vector key
                            atom_idcs = np.where(atom_mask)[0]

                            # Add edges from current metal to each valid donor atom
                            for atom_idx in atom_idcs:
                                graph.add_edge(metal_idx, atom_idx)

                        elif partner_key.startswith("metal_"):
                            # --- Handle metal–metal connections ---
                            # Extract the second metal's index from the string
                            other_metal_num = int(partner_key.split("_")[1]) - 1

                            # Get global atom index for this second metal
                            other_metal_idx = metal_idc[other_metal_num]

                            # Avoid adding duplicate edges by enforcing i < j
                            if metal_idx < other_metal_idx:
                                graph.add_edge(metal_idx, other_metal_idx)
        else:
            return LoggedValueError("Fatal Error: Alternative graph generation not implemented yet. ", logger)

        # --- Step 5: Add intra-ligand bonds (all atoms with same ligand tag > 0) ---
        for ligand_idx, ligand in enumerate(self.ligands):
            # Get atom indices in self.atoms that belong to this ligand
            atom_indices = np.where(tags[:, 1] == ligand_idx + 1)[0]  # ligand indices start at 1

            if len(atom_indices) != len(ligand.graph.nodes):
                raise LoggedValueError(f"Mismatch: ligand {ligand_idx} has {len(ligand.graph.nodes)} nodes but {len(atom_indices)} atoms in Atoms object.", logger)

            # Map ligand-local node indices to global atom indices
            index_map = dict(zip(ligand.graph.nodes, atom_indices))

            for u, v in ligand.graph.edges():
                graph.add_edge(index_map[u], index_map[v])

        graph = nx.relabel_nodes(graph, {n: int(n) for n in graph.nodes})

        # Add node labels for RCA_Molecule
        for i, atom in enumerate(self.atoms):
            graph.add_node(i, node_label=atom.symbol)

        return graph

    def _get_metal_idc(self) -> List[int]:
        """
        Get the indices of metal nodes from the ase object
        :return: List of indices of metal nodes in the graph
        """
        tags = self.atoms.get_array('multi_tags')
        metal_idc = np.where((tags == [0, 0, 0]).all(axis=1))[0]
        return metal_idc.tolist()

    def _get_ligand_idc(self):
        pass

    def _get_donor_idc(self):
        tags = self.atoms.get_array('multi_tags')
        Donor_idc = np.where((tags[:, 2] == [1]))[0]
        return Donor_idc.tolist()

    def to_dict(self):
        """
        Converts the DARTIsomer object to a fully JSON-serializable dictionary.
        """

        # Note: I seem to have encountered an issue with non serializable types in the graph data.
        # Note: This function sanitizes the graph data to ensure it can be serialized.
        def sanitize(obj):
            if isinstance(obj, dict):
                return {k: sanitize(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [sanitize(v) for v in obj]
            elif isinstance(obj, (np.integer, np.int64)):
                return int(obj)
            elif isinstance(obj, (np.floating, np.float64)):
                return float(obj)
            else:
                return obj

        raw_graph_data = json_graph.node_link_data(self.graph)
        clean_graph_data = sanitize(raw_graph_data)

        return {
            "atomic_props": {
                "x": [float(x) for x in self.atoms.positions[:, 0]],
                "y": [float(y) for y in self.atoms.positions[:, 1]],
                "z": [float(z) for z in self.atoms.positions[:, 2]],
                "atoms": [str(s) for s in self.atoms.get_chemical_symbols()],
            },
            "graph": clean_graph_data,
            "global_props": self.global_props,
            "metal_idc": [int(i) for i in self._get_metal_idc()],
            "donor_idc": [int(i) for i in self._get_donor_idc()],
            "ligand_idc": self._get_ligand_idc(),  # make sure this returns a serializable type too
            "ligand_info": sanitize(self.ligand_info),
        }

    def view_3d_graph_interactive(self) -> None:
        """
        Interactive 3D visualization of the ASE Atoms object and bonding graph using Plotly.
        """
        atoms = self.atoms
        graph = self._get_graph()
        pos = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        tags = atoms.get_array('multi_tags')

        # Color and size scheme using CPK colors
        colors = []
        sizes = []
        for symbol, tag in zip(symbols, tags):
            color = self.get_cpk_colors().get(symbol, 'grey')  # fallback to grey
            colors.append(color)

            # Size logic based on tags
            if np.all(tag == [0, 0, 0]):  # metal
                sizes.append(40)
            elif tag[2] == 1:  # donor
                sizes.append(30)
            else:
                sizes.append(20)
        # Atom coordinates
        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]

        # --- Main atom spheres (no text) ---
        node_trace = go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(size=sizes, color=colors, opacity=1.0),
            name='Atoms',
            hoverinfo='text',
            text=[f"{i}: {s}" for i, s in enumerate(symbols)]
        )

        # --- Optional labels (separate trace) ---
        label_trace = go.Scatter3d(
            x=x, y=y, z=z,
            mode='text+markers',  # Add markers to trigger legend
            text=[f"{i}: {s}" for i, s in enumerate(symbols)],
            textposition="top center",
            textfont=dict(size=12, color='black'),
            marker=dict(size=1, color='rgba(0,0,0,0)'),  # Fully transparent marker
            name='Atom Labels',
            showlegend=True,
            hoverinfo='none',
            visible='legendonly'  # Hidden by default, toggleable
        )

        edge_x_solid, edge_y_solid, edge_z_solid = [], [], []
        edge_x_dashed, edge_y_dashed, edge_z_dashed = [], [], []

        metal_idc = set(self._get_metal_idc())
        donor_idc = set(self._get_donor_idc())

        for u, v in graph.edges():
            if u in graph.nodes and v in graph.nodes:
                x0, y0, z0 = pos[u]
                x1, y1, z1 = pos[v]
                is_hashed = (
                        (u in metal_idc and v in donor_idc) or
                        (v in metal_idc and u in donor_idc) or
                        (u in metal_idc and v in metal_idc and u != v)
                )

                if is_hashed:
                    edge_x_dashed += [x0, x1, None]
                    edge_y_dashed += [y0, y1, None]
                    edge_z_dashed += [z0, z1, None]
                else:
                    edge_x_solid += [x0, x1, None]
                    edge_y_solid += [y0, y1, None]
                    edge_z_solid += [z0, z1, None]

        # Solid bonds (non donor-metal)
        edge_trace_solid = go.Scatter3d(
            x=edge_x_solid, y=edge_y_solid, z=edge_z_solid,
            mode='lines',
            line=dict(width=5, color='gray'),
            hoverinfo='none',
            name='Normal Bonds'
        )

        # Dashed donor-metal bonds
        edge_trace_dashed = go.Scatter3d(
            x=edge_x_dashed, y=edge_y_dashed, z=edge_z_dashed,
            mode='lines',
            line=dict(width=10, color='cornflowerblue', dash='solid'),
            hoverinfo='none',
            name='Donor–Metal Bonds'
        )

        # Assemble and render
        fig = go.Figure(data=[edge_trace_dashed, edge_trace_solid, node_trace, label_trace])
        fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
            ),
            margin=dict(l=0, r=0, t=0, b=0),
            showlegend=True,  # <-- Enable legend
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.2,
                xanchor="center",
                x=0.5
            ),
        )
        fig.show()

    def get_cpk_colors(self):
        """
        Returns a dictionary mapping element symbols to their CPK colors.
        :return:
        """
        return {
            'H': '#D3D3D3', 'D': '#FFFFC0', 'T': '#FFFFA0', 'He': '#D9FFFF', 'Li': '#CC80FF',
            'Be': '#C2FF00', 'B': '#FFB5B5', 'C': '#909090', 'C-13': '#505050', 'C-14': '#404040',
            'N': '#3050F8', 'N-15': '#105050', 'O': '#FF0D0D', 'F': '#90E050', 'Ne': '#B3E3F5',
            'Na': '#AB5CF2', 'Mg': '#8AFF00', 'Al': '#BFA6A6', 'Si': '#F0C8A0', 'P': '#FF8000',
            'S': '#FFFF30', 'Cl': '#1FF01F', 'Ar': '#80D1E3', 'K': '#8F40D4', 'Ca': '#3DFF00',
            'Sc': '#E6E6E6', 'Ti': '#BFC2C7', 'V': '#A6A6AB', 'Cr': '#8A99C7', 'Mn': '#9C7AC7',
            'Fe': '#E06633', 'Co': '#F090A0', 'Ni': '#50D050', 'Cu': '#C88033', 'Zn': '#7D80B0',
            'Ga': '#C28F8F', 'Ge': '#668F8F', 'As': '#BD80E3', 'Se': '#FFA100', 'Br': '#A62929',
            'Kr': '#5CB8D1', 'Rb': '#702EB0', 'Sr': '#00FF00', 'Y': '#94FFFF', 'Zr': '#94E0E0',
            'Nb': '#73C2C9', 'Mo': '#54B5B5', 'Tc': '#3B9E9E', 'Ru': '#248F8F', 'Rh': '#0A7D8C',
            'Pd': '#006985', 'Ag': '#C0C0C0', 'Cd': '#FFD98F', 'In': '#A67573', 'Sn': '#668080',
            'Sb': '#9E63B5', 'Te': '#D47A00', 'I': '#940094', 'Xe': '#429EB0', 'Cs': '#57178F',
            'Ba': '#00C900', 'La': '#70D4FF', 'Ce': '#FFFFC7', 'Pr': '#D9FFC7', 'Nd': '#C7FFC7',
            'Pm': '#A3FFC7', 'Sm': '#8FFFC7', 'Eu': '#61FFC7', 'Gd': '#45FFC7', 'Tb': '#30FFC7',
            'Dy': '#1FFFC7', 'Ho': '#00FF9C', 'Er': '#00E675', 'Tm': '#00D452', 'Yb': '#00BF38',
            'Lu': '#00AB24', 'Hf': '#4DC2FF', 'Ta': '#4DA6FF', 'W': '#2194D6', 'Re': '#267DAB',
            'Os': '#266696', 'Ir': '#175487', 'Pt': '#D0D0E0', 'Au': '#FFD123', 'Hg': '#B8B8D0',
            'Tl': '#A6544D', 'Pb': '#575961', 'Bi': '#9E4FB5', 'Po': '#AB5C00', 'At': '#754F45',
            'Rn': '#428296', 'Fr': '#420066', 'Ra': '#007D00', 'Ac': '#70ABFA', 'Th': '#00BAFF',
            'Pa': '#00A1FF', 'U': '#008FFF', 'Np': '#0080FF', 'Pu': '#006BFF', 'Am': '#545CF2',
            'Cm': '#785CE3', 'Bk': '#8A4FE3', 'Cf': '#A136D4', 'Es': '#B31FD4', 'Fm': '#B31FBA',
            'Md': '#B30DA6', 'No': '#BD0D87', 'Lr': '#C70066', 'Rf': '#CC0059', 'Db': '#D1004F',
            'Sg': '#D90045', 'Bh': '#E00038', 'Hs': '#E6002E', 'Mt': '#EB0026'
        }

"""
This file contains functions and a classes for the input of the assembly. They doublecheck that the input is correct and convert it to the correct format.
"""
import ast
from constants.Paths import project_path, default_ligand_db_path
from constants.Periodic_Table import all_atomic_symbols
from pathlib import Path, PurePath
from typing import Union, Any, Optional, Tuple, List, Dict
from src01.io_custom import read_yaml

# Define the key names in the assembly input file
_verbose = 'verbose'
_optimization_movie = 'optimization_movie'
_concatenate_xyz = 'concatenate_xyz'
_output_path = 'Output_Path'
_batches = 'Batches'
_name = 'Name'
_input_path = 'Input_Path'
_max_num_complexes = 'MAX_num_complexes'
_isomers = 'Isomers'
_optimisation = 'Optimisation_Choice'
_random_seed = 'Random_Seed'
_total_charge = 'Total_Charge'
_topology = 'Topology'
_metal = 'Metal'
_element = 'element'
_oxidation_state = 'oxidation_state'
_spin = 'spin'


class AssemblyInput(object):
    """
    Class that contains the settings for the assembly.
    """

    # Define valid keys and their acceptable types in the input
    # Global settings
    valid_keys = {
                        _verbose: [int, str],
                        _optimization_movie: [bool, str],
                        _concatenate_xyz: [bool, str],
                        _output_path: [str, Path],
                        _batches: [list, tuple, dict]
                        }
    # Batch settings
    batches_valid_keys = {
                        _name: [str],
                        _input_path: [str, type(None)],
                        _max_num_complexes: [int, str],
                        _isomers: [str],
                        _optimisation: [str, bool],
                        _random_seed: [int, str],
                        _total_charge: [int, str],
                        _topology: [str, list, tuple],
                        _metal: [dict],
                        }
    # Metal settings in the batch settings
    metal_valid_keys = {
                        _element: [str],
                        _oxidation_state: [int, str],
                        _spin: [int, str],
                        }

    def __init__(self, path: Union[str, Path] = 'assembly_input.yml'):
        """
        Reads the global settings file and stores the settings in the class.
        """
        # Check if the file exists and make to Path
        self.path = self.ensure_assembly_input_file_present(path)

        # Set the batch name to None. It will be set later to the respective batch name when iterating over the batches so that the error messages can be more specific.
        self.batch_name = None

        # Read into dict
        self.global_settings = read_yaml(self.path)

        # Check the input and set the class variables
        self.verbose, self.optimization_movie, self.concatenate_xyz, self.Output_Path, self.Batches = self.check_and_return_global_settings()

    def check_and_return_global_settings(self) -> tuple[int, bool, bool, Path, list]:
        """
        Checks the global settings and returns them.
        """
        for key, types in self.valid_keys.items():
            if key not in self.global_settings.keys():
                self.raise_error(f"Key '{key}' not found in input file. Please add it.")

            self.check_correct_input_type(input=self.global_settings[key], types=types, varname=key)

        # Check if the output path exists and is a directory
        Output_Path = self.ensure_output_directory_present(self.global_settings[_output_path], varname=_output_path)

        verbose = self.get_int_from_input(self.global_settings[_verbose], varname=_verbose)
        optimization_movie = self.get_bool_from_input(self.global_settings[_optimization_movie], varname=_optimization_movie)
        concatenate_xyz = self.get_bool_from_input(self.global_settings[_concatenate_xyz], varname=_concatenate_xyz)
        Batches =  self.get_batches_from_input(self.global_settings[_batches])

        return verbose, optimization_movie, concatenate_xyz, Output_Path, Batches

    def get_batches_from_input(self, batches: Union[list, tuple, dict]):
        """
        Checks if the batches input is correct.
        """
        if isinstance(batches, (list, tuple)):
            batches = list(batches)
        elif isinstance(batches, dict):
            batches = [batches]
        else:
            self.raise_error(f"Input '{_batches}' must be a list/tuple of dicts or a dict, but is {type(batches)}.")

        return batches

    def check_batches_input(self):
        """
        Checks the batch settings for errors and raises errors if the settings are not correct. This check is done already here so that potential errors are raised before the assembly starts.
        """
        all_batch_names = []
        for batch_settings in self.Batches:
            batch_name, *_ = self.check_and_return_batch_settings(batch_settings)
            all_batch_names.append(batch_name)
        self.batch_name = None  # Reset the batch name to None so that the error messages are not specific to a batch

        # Check if the batch names are unique
        if len(all_batch_names) != len(set(all_batch_names)):
            self.raise_error(f"Batch names must be unique but are not. Batch names are: {all_batch_names}", varname=_name)

        return

    def return_batch_settings(self, batch_index: int):
        """
        Returns the settings for a single batch.
        """
        batch_settings = self.Batches[batch_index]
        return self.check_and_return_batch_settings(batch_settings)

    def check_and_return_batch_settings(self, batch_settings: dict):
        """
        Checks the batch settings for a single batch. Raises errors if the settings are not correct and returns the settings in the correct format.
        """
        # Get name of batch to make error messages more specific
        try:
            self.batch_name = str(batch_settings[_name])
        except KeyError:
            self.raise_error(f"Key '{_name}' not found in input file. Please add it.")

        # Check if all keys are present
        for key, types in self.batches_valid_keys.items():
            if key not in batch_settings.keys():
                self.raise_error(f"Key '{key}' not found in input file. Please add it.")
            varname = f"{_batches}->{key}"
            self.check_correct_input_type(input=batch_settings[key], types=types, varname=varname)

        # Here we take the batch inputs and format them correctly
        Ligand_json = self.get_ligand_db_path_from_input(batch_settings[_input_path])
        Max_Num_Assembled_Complexes = self.get_int_from_input(batch_settings[_max_num_complexes], varname=f'{_batches}->{_max_num_complexes}')
        Generate_Isomer_Instruction = self.get_isomers_from_input(batch_settings[_isomers])
        Optimisation_Instruction = self.get_bool_from_input(batch_settings[_optimisation], varname=f'{_batches}->{_optimisation}')
        Random_Seed = self.get_int_from_input(batch_settings[_random_seed], varname=f'{_batches}->{_random_seed}')
        Total_Charge = self.get_int_from_input(batch_settings[_total_charge], varname=f'{_batches}->{_total_charge}')
        metal_list = self.get_metal_from_input(batch_settings[_metal])
        topology_similarity = self.get_topology_from_input(batch_settings[_topology])

        return self.batch_name, Ligand_json, Max_Num_Assembled_Complexes, Generate_Isomer_Instruction, Optimisation_Instruction, Random_Seed, Total_Charge, metal_list, topology_similarity

    def get_ligand_db_path_from_input(self, ligand_db_path):
        """
        Checks the input for the ligand database path.
        """
        varname = f'{_batches}->{_input_path}'
        if ligand_db_path is None or ligand_db_path == '' or ligand_db_path == 'default':
            print(f"Ligand database path '{varname}' is not specified in batch {self.batch_name} in file '{self.path}'. Using full default ligand database.")
            ligand_db_path = default_ligand_db_path
        else:
            # Check if the file exists
            ligand_db_path = self.ensure_file_present(ligand_db_path, varname=varname)

        return Path(ligand_db_path)

    def get_metal_from_input(self, metal: dict):
        """
        Checks the metal input for correct input.
        """
        varname = f'{_batches}->{_metal}'
        for key, types in self.metal_valid_keys.items():
            if key not in metal.keys():
                self.raise_error(f"Key '{key}' is missing in '{_metal}'.", varname=varname)
            self.check_correct_input_type(input=metal[key], types=types, varname=f'{varname}->{key}')

        element, oxidation_state, spin = metal[_element], metal[_oxidation_state], metal[_spin]

        # Check the element
        if not is_chemical_element(element):
            self.raise_error(f"Input element '{element}' is not a valid chemical element symbol, e.g. 'Fe' for iron.", varname=f'{varname}->{_element}')
        element = str(element)

        # Check the oxidation state
        oxidation_state = self.get_int_from_input(oxidation_state, varname=f'{varname}->{_oxidation_state}')
        if oxidation_state > 0:
            oxidation_state = f"+{oxidation_state}"
        else:
            self.raise_error(f"Input oxidation state '{oxidation_state}' is not a positive integer > 0.", varname=f'{varname}->{_oxidation_state}')

        # Check the spin
        spin = self.get_int_from_input(spin, varname=f'{varname}->{_spin}')
        if spin < 0:
            self.raise_error(f"Input spin '{spin}' is not a positive integer >= 0.", varname=f'{varname}->{_spin}')
        spin = str(spin)

        return [element, oxidation_state, spin]

    def get_topology_from_input(self, topology: str):
        """
        Checks the topology input for correct input.
        """
        topology = str(topology)
        error_message = f"Topology '{topology}' is not in the correct format. It should be of format '(1,1,2,2)' for specifying denticities or of format (1,1,2,2)--(1,1,2,3) for specifying both denticities and similarities."
        varname = f'{_batches}->{_topology}'

        splits = topology.split('--')
        if len(splits) == 1:    # only denticities given
            try:
                denticities = ast.literal_eval(topology)
            except (ValueError, SyntaxError):
                # If the input is not a list of integers, then we raise an error
                self.raise_error(error_message, varname=varname)
            similarities = list(range(1, len(denticities) + 1)) # Default similarities are 1, 2, 3, i.e. all are different

        elif len(splits) == 2:      # both denticities and similarities are given
            try:
                denticities, similarities = splits[0], splits[1]
                denticities, similarities = ast.literal_eval(denticities), ast.literal_eval(similarities)
            except (ValueError, SyntaxError):  # If the input is not a list of integers, then we raise an error
                self.raise_error(error_message, varname=varname)

        else:        # more than one '--' in the input
            self.raise_error(error_message, varname=varname)

        # Check that denticities and similarities are either lists or tuples and make them to lists for the rest of the code
        if not isinstance(denticities, (list, tuple)):
            self.raise_error(error_message, varname=varname)
        if not isinstance(similarities, (list ,tuple)):
            self.raise_error(error_message, varname=varname)
        denticities, similarities = list(denticities), list(similarities)

        # Check that denticities and similarities have the same length and that they are positive integers
        if len(denticities) != len(similarities):
            self.raise_error(f"Topology '{topology}' has different number of denticities and similarities.", varname=varname)
        if not all(isinstance(denticity, int) and denticity > 0 for denticity in denticities):
            self.raise_error(f"Topology '{topology}' has denticities that are not positive integers.", varname=varname)
        if not all(isinstance(similarity, int) and similarity > 0 for similarity in similarities):
            self.raise_error(f"Topology '{topology}' has similarities that are not positive integers.", varname=varname)

        output_topology = str(denticities) + '--' + str(similarities)
        return output_topology

    def get_isomers_from_input(self, isomers):
        """
        Checks if the isomers input is correct.
        """
        possible_values = ['Generate Lowest Energy', 'Generate All']
        isomers = str(isomers)
        if isomers not in possible_values:
            self.raise_error(f" Input value '{isomers}' not recognized. It must be one of {possible_values} (case sensitive).", varname=f'{_batches}->{_isomers}')

        return isomers
    
    def check_correct_input_type(self, input, types: list, varname: str) -> Any:
        """
        Checks if the input is of the correct type and raises an error if not.
        """
        if not any(isinstance(input, type) for type in types):
            input_type = type(input).__name__
            types = tuple([type.__name__ for type in types])
            self.raise_error(f"Input '{input}' is not any the allowed types {types} but of type '{input_type}'", varname=varname)

        return input
    
    def ensure_output_directory_present(self, path: Union[str, Path], varname: str) -> Path:
        """
        Checks if the output directory exists and creates it if it does not.
        """
        try:
            path = Path(path)
        except TypeError:
            self.raise_error(message=f"Output directory '{path}' is not a valid string.", varname=varname)

        path = path.resolve()   # get absolute path
        if not path.is_dir():
            print(f"Output directory '{path}' does not exist. Creating it now.")
            path.mkdir(parents=True, exist_ok=True)

        return path
    
    
    def ensure_assembly_input_file_present(self, path: Union[str, Path]) -> Path:
        """
        Checks if the path to the input file exist and the path is a file. Raises an error if not.
        """
        try:
            path = Path(path)
        except TypeError:
            raise TypeError(f"The input filepath '{path}' is not a valid string.")

        if not path.is_file():
            raise FileNotFoundError(f"The input filepath '{path}' either doesn't exist or is not a file.")

        path = path.resolve()   # get absolute path
        return path

    def ensure_file_present(self, path: Union[str, Path], varname: str) -> Path:
        """
        Checks if the path to the file exist and the path is a file.
        """
        try:
            path = Path(path)
        except TypeError:
            self.raise_error(message=f"The input filepath '{path}' is not a valid string.", varname=varname)

        if not path.is_file():
            self.raise_error(message=f"The input filepath '{path}' either doesn't exist or is not a file.", varname=varname)

        path = path.resolve()  # get absolute path
        return path
    
    
    def get_bool_from_input(self, input: Union[str, bool], varname: str) -> bool:
        """
        Returns a bool from a string or bool input.
        """
        meaning = str(input).lower()  # Convert bool/string to lowercase string
        if meaning == 'true':
            meaning = True
        elif meaning == 'false':
            meaning = False
        else:
            self.raise_error(message=f"Input '{input}' can not be recognized as bool. Valid inputs are the words 'True' or 'False' in uppercase/lowercase.", varname=varname)

        return meaning
    
    def get_int_from_input(self, input: Union[str, int], varname: str) -> int:
        """
        Returns an int from a string or int input.
        """
        if isinstance(input, int):
            output = input
        elif isinstance(input, bool):
            self.raise_error(message=f"Input '{input}' is of type bool, but int is expected.", varname=varname)
        else:
            try:
                output = int(input)
            except ValueError:
                self.raise_error(message=f"Input '{input}' can not be recognized as int.", varname=varname)

        return output
    
    def raise_error(self, message: str, varname: str='', file=None):
        """
        Raises an AssemblyInputError with the given message and the path to the input file.
        """
        if file is None:
            file = self.path
        batch_name = self.batch_name or ''  # If batch_name is None, set it to empty string. Otherwise enrich the error message with the batch name.

        raise AssemblyInputError(message=message, varname=varname, file=file, batch_name=batch_name)


def is_chemical_element(element: str) -> bool:
    """
    Checks if the input is a valid chemical element.
    """
    return element in all_atomic_symbols


class AssemblyInputError(Exception):
    """
    Exception raised for errors in the input.
    """
    def __init__(self, message: str, varname: str='', file: str='', batch_name:str =''):
        if varname != '':
            varname = f" for key '{varname}'"
        if batch_name != '':
            batch_name = f" in batch '{batch_name}'"
        if file != '':
            file = f" in input file '{file}'"
        total_message = f"Invalid input{varname}{batch_name}{file}: {message}"
        super().__init__(total_message)
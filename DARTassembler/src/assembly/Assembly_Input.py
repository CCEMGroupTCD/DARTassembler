"""
This file contains functions and a classes for the input of the assembly. They doublecheck that the input is correct and convert it to the correct format.
"""
import difflib
import warnings
from copy import deepcopy
import random

import yaml
import ase
import pandas as pd

from DARTassembler.src.constants.Periodic_Table import all_atomic_symbols
from DARTassembler.src.ligand_extraction.composition import Composition
from pathlib import Path
from typing import Union, Any, Tuple
from DARTassembler.src.ligand_extraction.io_custom import read_yaml
from DARTassembler.src.ligand_extraction.utilities_Molecule import stoichiometry2atomslist
from DARTassembler.src.metalig.metalig_utils import get_correct_ligand_db_path_from_input

allowed_topologies = ['mer-3-2-1', 'mer-4-1-1', '5-1', '2-1-1', '2-2']  # list all allowed topologies here, others will be rejected
allowed_verbosities = [0, 1, 2, 3]

# Define the key names in the assembly input file
# Global settings
_verbose = 'verbosity'    # verbosity
_optimization_movie = 'ffmovie'
_concatenate_xyz = 'concatenate_xyz'
_output_path = 'output_directory'
_complex_name_length = 'complex_name_length'
_same_isomer_names = 'same_isomer_names'
# Batch settings
_batches = 'batches'
_name = 'name'
_input_path = 'ligand_db_file'
_max_num_complexes = 'max_num_complexes'
_ligand_choice = 'ligand_choice'
_isomers = 'isomers'
_optimisation = 'forcefield'
_random_seed = 'random_seed'
_total_charge = 'total_charge'
_topology = 'geometry'
_metal = 'metal'
_element = 'metal_center'
_oxidation_state = 'metal_oxidation_state'
_spin = 'spin'
_complex_name_appendix = 'complex_name_appendix'
_geometry_modifier_filepath = 'geometry_modifier_filepath'
_bidentate_rotator = 'bidentate_rotator'
# _ligand_filters_path = 'ligand_filters_path'
_gaussian_path = 'gaussian_input_filepath'
# Others
_same_ligand_keyword = 'same_ligand_as_previous'


# Define names for the settings in the ligand filters file:
# Global
_ligand_db_path = 'input_db_file'
_output_ligand_db_path = 'output_db_file'
_output_ligands_info = 'output_ligands_info'
_filters = 'filters'
# Filters
_filter = 'filter'
_graph_hash_wm = 'graph_IDs'     # graph_hash_with_metal
_denticities = 'apply_to_denticities'
_denticities_of_interest = 'denticities'
_remove_ligands_with_neighboring_coordinating_atoms = 'remove_ligands_with_adjacent_coordinating_atoms'
_only_confident_charges = 'only_confident_charges'
_remove_ligands_with_beta_hydrogens = 'remove_ligands_with_beta_hydrogens'
_strict_box_filter= 'strict_box_filter'
_filter_even_odd_electron_count = 'remove_even_odd_electron_count'
_dentfilters = 'denticity_dependent_filters'
_atm_neighbors = 'atomic_neighbors'
_atom = 'atom'
_neighbors = 'neighbors'
# Denticity dependent filters
_acount = 'number_of_atoms'
_acount_min = 'min'
_acount_max = 'max'

_stoichiometry = 'stoichiometry'

_ligcomp = 'ligand_composition'
_ligcomp_atoms_of_interest = 'elements'
_ligcomp_instruction = 'instruction'

_ligand_charges = 'ligand_charges'

_metals_of_interest = 'metal_ligand_binding_history'

_coords = 'coordinating_atoms_composition'
_coords_atoms_of_interest = 'elements'
_coords_instruction = 'instruction'

_mw = 'molecular_weight'
_mw_min = 'min'
_mw_max = 'max'
_interatomic_distances = 'interatomic_distances'
_planarity = 'planarity'
_occurrences = 'occurrences'
_min = 'min'
_max = 'max'
_remove_missing_bond_orders = 'remove_ligands_with_missing_bond_orders'
_smarts_filter = 'smarts'
_smarts = 'smarts'
_should_be_present = 'should_contain'
_include_metal = 'include_metal'

assembler_gbl_defaults = {
    _verbose: 2,
    _optimization_movie: False,
    _concatenate_xyz: True,
    _output_path: 'DART',
    _complex_name_length: 8,
    _same_isomer_names: True,
}
def get_random_random_seed(start=1_000, end=9_999) -> int:
    return random.randint(start, end)

assembler_batch_defaults = {
    _optimisation: False,
    _random_seed: get_random_random_seed(),
    _complex_name_appendix: None,
    _geometry_modifier_filepath: None,
    _bidentate_rotator: 'auto',
}
ligandfilters_gbl_defaults = {
    _ligand_db_path: 'metalig',
    _output_ligand_db_path: None,
    _output_ligands_info: True,
}
ligandfilters_filter_defaults = {
    _graph_hash_wm: {},
    _denticities_of_interest: {},
    _remove_ligands_with_neighboring_coordinating_atoms: {_denticities: None},
    _only_confident_charges: {},
    _remove_ligands_with_beta_hydrogens: {_denticities: None},
    _strict_box_filter: {},
    _filter_even_odd_electron_count: {},
    _acount: {_denticities: None},
    _stoichiometry: {_denticities: None},
    _ligcomp: {_denticities: None},
    _ligand_charges: {_denticities: None},
    _metals_of_interest: {_denticities: None},
    _coords: {_denticities: None},
    _mw: {_denticities: None},
    _interatomic_distances: {_denticities: None},
    _planarity: {_denticities: None},
    _occurrences: {_denticities: None},
    _remove_missing_bond_orders: {_denticities: None},
    _atm_neighbors: {_denticities: None},
    _smarts_filter: {_denticities: None},
}

class BaseInput(object):

    def __init__(self, path: Union[str, Path]):
        self.path = path

    @classmethod
    def get_settings_from_input_file(cls, path: Union[str, Path]):
        """
        Alternative constructor that reads the input from a file.
        """
        path = cls.ensure_assembly_input_file_present(path)
        global_settings = read_yaml(path)

        return global_settings

    def read_yaml_and_fix_indentation_issues(self, path, indentations: dict = {}) -> dict:
        """
        Reads a yaml file and fixes indentation issues.
        """
        with open(path, 'r') as file:
            text = file.read()

        for key, value in indentations.items():
            text = text.replace(key, value)

        return yaml.safe_load(text)

    def check_correct_input_type(self, input, types: list, varname: str) -> Any:
        """
        Checks if the input is of the correct type and raises an error if not.
        """
        if not any(isinstance(input, type) for type in types):
            input_type = type(input).__name__
            types = tuple([type.__name__ for type in types])
            self.raise_error(f"Input '{input}' is not any the allowed types {types} but of type '{input_type}'",
                             varname=varname)

        return input

    def ensure_output_directory_valid(self, path: Union[str, Path], varname: str) -> Path:
        """
        Checks if the output directory is valid.
        """
        try:
            path = Path(path)
        except TypeError:
            self.raise_error(message=f"Output directory '{path}' is not a valid string.", varname=varname)

        path = path.resolve()  # get absolute path

        return path

    @classmethod
    def ensure_assembly_input_file_present(cls, path: Union[str, Path]) -> Path:
        """
        Checks if the path to the input file exist and the path is a file. Raises an error if not.
        """
        try:
            path = Path(path)
        except TypeError:
            raise TypeError(f"The input filepath '{path}' is not a valid string.")

        if not path.is_file():
            raise FileNotFoundError(f"The input filepath '{path}' either doesn't exist or is not a file.")

        path = path.resolve()  # get absolute path
        return path

    def get_path_from_input(self, path: Union[str, Path], varname: str, allow_none=False) -> Union[Path,None]:
        """
        Checks if the path to the file exist and the path is a file.
        """
        if allow_none and path is None:
            return None

        try:
            path = Path(path)
        except TypeError:
            self.raise_error(message=f"The input filepath '{path}' is not a valid path.", varname=varname)

        path = path.resolve()   # get absolute path
        return path

    def ensure_file_present(self, path: Union[str, Path], varname: str, allow_none:bool= False) -> Path:
        """
        Checks if the path to the file exist and the path is a file.
        """
        path = self.get_path_from_input(path, varname=varname, allow_none=allow_none)

        if not path.is_file():
            self.raise_error(message=f"The input filepath '{path}' either doesn't exist or is not a file.",
                             varname=varname)

        return path

    def get_bool_from_input(self, input: Union[str, bool], varname: str, allow_none=False) -> Union[bool,None]:
        """
        Returns a bool from a string or bool input.
        """
        if allow_none and input is None:
            return None

        meaning = str(input).lower()  # Convert bool/string to lowercase string
        if meaning == 'true':
            meaning = True
        elif meaning == 'false':
            meaning = False
        else:
            self.raise_error(
                message=f"Input '{input}' can not be recognized as bool. Valid inputs are the words 'True' or 'False' in uppercase/lowercase.",
                varname=varname)

        return meaning

    def get_int_from_input(self, input: Union[str, int], varname: str, allow_none = False, allow_str: str = None, allow=None) -> Union[int,None]:
        """
        Returns an int from a string or int input.
        """
        if allow_none and input is None:    # If input is None, return None
            return None

        if allow_str is not None and input == allow_str:  # If input is a specific string, return this string
            return input

        if isinstance(input, int):
            output = input
        elif isinstance(input, bool):
            self.raise_error(message=f"Input '{input}' is of type bool, but int is expected.", varname=varname)
        else:
            try:
                output = int(input)
            except ValueError:
                self.raise_error(message=f"Input '{input}' can not be recognized as int.", varname=varname)

        if not allow is None:
            if not output in allow:
                self.raise_error(message=f"Input '{input}' is not allowed. Allowed options: {allow}", varname=varname)

        return output

    def get_float_from_input(self, input: Union[str, int], varname: str, allow_none=False) -> Union[float,None]:
        """
        Returns an int from a string or int input.
        """
        if allow_none and input is None:
            return None

        if isinstance(input, float):
            output = input
        elif isinstance(input, bool):
            self.raise_error(message=f"Input '{input}' is of type bool, but float is expected.", varname=varname)
        else:
            try:
                output = float(input)
            except ValueError:
                self.raise_error(message=f"Input '{input}' can not be recognized as float.", varname=varname)

        return output

    def get_list_of_chemical_elements_from_input(self, input: Union[str, list, tuple], varname: str, allow_none=False) -> Union[list,None]:
        """
        Returns a list of elements from a string, list or tuple input.
        """
        if allow_none and input is None:
            return None

        if isinstance(input, str): # Convert stoichiometry to list of elements
            input = stoichiometry2atomslist(input)

        input = list(input)
        for i in range(len(input)):
            el = input[i]
            if not is_chemical_element(str(el)):
                self.raise_error(message=f"Input '{el}' can not be recognized as chemical element.", varname=varname)
            input[i] = str(el)

        return input

    def get_list_of_ints_from_input(self, input: Union[str, list, tuple], varname: str, allow_none=False) -> Union[list,None]:
        """
        Returns a list of ints from a string, list or tuple input.
        """
        if allow_none and input is None:
            return None

        if isinstance(input, (str,int)):
            input = [input]

        input = list(input)
        for i in range(len(input)):
            input[i] = self.get_int_from_input(input=input[i], varname=varname)

        return input

    def get_list_of_str_from_input(self, input: Union[str, list, tuple], allow_none=False) -> Union[list,None]:
        """
        Returns a list of str from a string, list or tuple input.
        """
        if allow_none and input is None:
            return None

        if isinstance(input, (str,int)):
            input = [input]

        input = list(input)
        for i in range(len(input)):
            input[i] = str(input[i])

        return input

    def get_instruction_from_input(self, input: Union[str, list, tuple], varname: str, valid_instructions: Union[list,None] = None, allow_none=False) -> Union[str,None]:
        if allow_none and input is None:
            return None

        default_valid_instructions = ['must_contain_and_only_contain', 'must_at_least_contain', 'must_exclude', 'must_only_contain_in_any_amount']
        if valid_instructions is None:
            valid_instructions = default_valid_instructions

        input = str(input).lower()
        if input not in valid_instructions:
            self.raise_error(message=f"Input '{input}' is not a valid instruction. Valid instructions are {valid_instructions}.", varname=varname)

        return input

    def check_input_types(self, valid_keys: dict, settings: dict):
        # Check for keys that are not allowed
        for key in settings.keys():
            if key not in valid_keys:
                self.raise_error(f"Key '{key}' is not a valid key. Please check the documentation.", varname=key)

        # Check if all necessary keys are specified
        for key, types in valid_keys.items():
            real_keys = tuple(settings.keys())
            if key not in real_keys:
                self.raise_error(f"Key '{key}' not found in input file, please add it. All keys found are {real_keys}.")
            self.check_correct_input_type(input=settings[key], types=types, varname=key)

        return

    def check_dict_is_fully_specified(self, d: dict, varname: str):
        """
        Checks that all values in the dict are either all None or all not None.
        """
        n_none = len([value for value in d.values() if value is None])
        n_total = len(d)
        if n_none != 0 and n_none != n_total:
            self.raise_error(message=f"Unspecified values in dict '{varname}' even though some other values are specified. Please specify either all values to use the filter or none to disable it.", varname=varname)

    def check_if_settings_not_recognized(self, actual_settings, valid_settings: dict):
        """
        Check if there are any unrecognized settings and raises a warning if so.
        """
        for actual_key in actual_settings:
            if not actual_key in valid_settings:
                self.raise_error(message=f"The provided word '{actual_key}' is not a valid DART setting. Please check the documentation", varname=actual_key)

        for batch in actual_settings[_batches]:
            for actual_key in batch:
                if not actual_key in valid_settings[_batches]:
                    self.raise_error(message=f"The provided word '{actual_key}' is not a valid DART setting. Please check the documentation.", varname=actual_key)

        return

    def raise_warning(self, message: str, varname: str = '', file=None):
        """
        Raises a warning with the given message and the path to the input file.
        """
        if file is None:
            file = self.path
        batch_name = self.batch_name if hasattr(self, 'batch_name') and self.batch_name is not None else ''

        if varname != '':
            varname = f" for key '{varname}'"
        if batch_name != '':
            batch_name = f" in batch '{batch_name}'"
        if file != '':
            file = f" in input file '{file}'"
        total_message = f"Invalid input{varname}{batch_name}{file}:\n\t\t{message}"

        warnings.warn(total_message, UserWarning)

    def raise_error(self, message: str, varname: str = '', file=None):
        """
        Raises an AssemblyInputError with the given message and the path to the input file.
        """
        if file is None:
            file = self.path

        if hasattr(self, 'batch_name') and self.batch_name is not None:
            batch_name = f" in batch '{self.batch_name}'"
        elif hasattr(self, 'current_denticity') and self.current_denticity is not None:
            batch_name = f" in denticity '{self.current_denticity}'"
        elif hasattr(self, 'filtername') and self.filtername is not None:
            batch_name = f" in filter '{self.filtername}'"
        else:
            batch_name = ''

        raise AssemblyInputError(message=message, varname=varname, file=file, batch_name=batch_name)


class LigandFilterInput(BaseInput):

    # Allowed keys and types for the input file.
    # Every list of allowed types must either be [dict] or include type(None).
    valid_keys = {
        _ligand_db_path: [str, Path, type(None)],
        _output_ligand_db_path: [str, Path, type(None)],
        _output_ligands_info: [bool, str],
        _filters: [list, tuple],
    }

    filter_keys = {
        _graph_hash_wm: {
            _graph_hash_wm: [str, list, tuple]},
        _denticities_of_interest: {
            _denticities_of_interest: [list, tuple]},
        _remove_ligands_with_neighboring_coordinating_atoms: {
            _remove_ligands_with_neighboring_coordinating_atoms: [bool, str],
            _denticities: [list, tuple, type(None), int]},
        # _only_confident_charges: {                    # filter removed from options and made mandatory
        #     _only_confident_charges: [bool, str]},
        _remove_ligands_with_beta_hydrogens: {
            _remove_ligands_with_beta_hydrogens: [bool, str],
            _denticities: [list, tuple, type(None), int]},
        # _strict_box_filter: {
        #     _strict_box_filter: [bool, str]},
        # _filter_even_odd_electron_count: {            # filter removed from options because all ligands with confident charges have even electron count
        #     _filter_even_odd_electron_count: [str]},
        _acount: {
            _acount_min: [int, str, type(None)],
            _acount_max: [int, str, type(None)],
            _denticities: [list, tuple, type(None), int],
            },
        _stoichiometry: {
            _stoichiometry: [str],
            _denticities: [list, tuple, type(None), int],
            },
        _ligcomp: {
            _ligcomp_atoms_of_interest: [list, tuple, str],
            _ligcomp_instruction: [str],
            _denticities: [list, tuple, type(None), int],
            },
        _ligand_charges: {
            _ligand_charges: [list, tuple],
            _denticities: [list, tuple, type(None), int],
            },
        _metals_of_interest: {
            _metals_of_interest: [list, tuple],
            _denticities: [list, tuple, type(None), int],
            },
        _coords: {
            _coords_atoms_of_interest: [list, tuple, str],
            _coords_instruction: [str],
            _denticities: [list, tuple, type(None), int],
            },
        _mw: {
            _mw_min: [float, str, type(None)],
            _mw_max: [float, str, type(None)],
            _denticities: [list, tuple, type(None), int],
            },
        _interatomic_distances: {
            _min: [float, str, type(None)],
            _max: [float, str, type(None)],
            _denticities: [list, tuple, type(None), int],
            },
        _planarity: {
            _min: [float, str, type(None)],
            _max: [float, str, type(None)],
            _denticities: [list, tuple, type(None), int],
            },
        _occurrences: {
            _min: [int, str, type(None)],
            _max: [int, str, type(None)],
            _denticities: [list, tuple, type(None), int],
            },
        _remove_missing_bond_orders: {
            _remove_missing_bond_orders: [bool, str],
            _denticities: [list, tuple, type(None), int],
            },
        _atm_neighbors: {
            _atom: [str],
            _neighbors: [list, tuple, str],
            _denticities: [list, tuple, type(None), int],
            },
        _smarts_filter: {
            _smarts: [str],
            _should_be_present: [bool, str],
            _include_metal: [bool, str],
            _denticities: [list, tuple, type(None), int],
            },
        }


    def __init__(self, path: Union[str, Path]):
        """
        Class for reading and checking the input file for the ligand filter. The input file should be a yaml file.
        """
        self.path = None
        super().__init__(path)

        self.raw_input_settings = self.get_settings_from_input_file(path)

        self.ligand_db_path = None
        self.output_ligand_db_path = None
        self.filters = None
        self.settings = self.set_and_check_settings()    # Set all settings and check if they are valid

    def _get_null_value_of_filter(self, allowed_types: list):
        """
        Returns the null value of a filter depending on the allowed types. There are only two cases: If the filter is a dict containing other filters, the null value must be set to an empty dicts, so that one can iterate over this filter. If the filter is not a dict, the null value must be None, so that one can check if the filter is set or not.
        """
        if allowed_types == [dict]:
            default = {}
        elif type(None) in allowed_types:
            default = None
        else:
            raise Warning(
                f"Implementation Error: Allowed types must either be [dict] or include type(None).")

        return default


    def set_and_check_settings(self) -> dict:

        settings = self.raw_input_settings

        # Add default values to settings
        for key, value in ligandfilters_gbl_defaults.items():
            if key not in settings:
                settings[key] = deepcopy(value)

        self.check_input_types(valid_keys=self.valid_keys, settings=settings)
        self.ligand_db_path = self.check_ligand_db_path(settings[_ligand_db_path])
        self.output_ligand_db_path = self.get_path_from_input(path=settings[_output_ligand_db_path], varname=_output_ligand_db_path, allow_none=True)
        self.output_filtered_ligands = self.get_bool_from_input(input=settings[_output_ligands_info], varname=_output_ligands_info)
        self.output_ligand_db_path = self.output_ligand_db_path or self.get_output_ligand_db_path()

        # Check all input types
        self.filters = self.check_filters(all_filters=settings[_filters])
        out_settings = {_ligand_db_path: self.ligand_db_path,
                        _output_ligand_db_path: self.output_ligand_db_path,
                        _filters: self.filters
                        }

        return out_settings

    def check_filters(self, all_filters: list) -> dict:
        out_settings = []
        for full_filter in all_filters:

            try:
                self.filtername = str(full_filter[_filter])
            except KeyError:
                self.raise_error(f"Key '{_filter} is missing in filter {full_filter}.")

            # Check for valid filter name
            valid_filter_names = tuple(self.filter_keys.keys())
            if not self.filtername in valid_filter_names:
                similar_word = get_closest_word(word=self.filtername, words=valid_filter_names)
                similar_string = f"Did you mean '{similar_word}'? " if similar_word != '' else ''
                self.raise_error(f"Filter '{self.filtername}' is not a valid filter. {similar_string}Valid filter names are {valid_filter_names}.", varname=_filter)

            # Add default values to filter
            default_values = ligandfilters_filter_defaults[self.filtername]
            for key, value in default_values.items():
                if key not in full_filter:
                    full_filter[key] = deepcopy(value)

            filter_values = {key: value for key, value in full_filter.items() if key != _filter}
            self.check_input_types(valid_keys=self.filter_keys[self.filtername], settings=filter_values)

            out_filter_settings = {_filter: self.filtername}
            if self.filtername == _denticities_of_interest:
                out_filter_settings[_denticities_of_interest] = self.check_denticities_of_interest(settings=filter_values)
            elif self.filtername == _graph_hash_wm:
                out_filter_settings[_graph_hash_wm] = self.get_list_of_str_from_input(input=filter_values[_graph_hash_wm])
            elif self.filtername == _only_confident_charges:
                out_filter_settings[_only_confident_charges] = self.get_bool_from_input(input=filter_values[_only_confident_charges], varname=_only_confident_charges)
            elif self.filtername == _strict_box_filter:
                out_filter_settings[_strict_box_filter] = self.get_bool_from_input(input=filter_values[_strict_box_filter], varname=_strict_box_filter)
            elif self.filtername == _filter_even_odd_electron_count:
                out_filter_settings[_filter_even_odd_electron_count] = self.check_even_odd_electron_count_input(settings=filter_values)

            # Denticity dependent filters
            elif self.filtername == _stoichiometry:
                out_filter_settings[_stoichiometry] = self.check_stoichiometry_input(stoichiometry=filter_values[_stoichiometry])
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_stoichiometry}:{_denticities}', allow_none=True)
            elif self.filtername == _ligcomp:
                out_filter_settings[_ligcomp_atoms_of_interest] = self.get_list_of_chemical_elements_from_input(input=filter_values[_ligcomp_atoms_of_interest], varname=f'{_ligcomp}:{_ligcomp_atoms_of_interest}')
                out_filter_settings[_ligcomp_instruction] = self.get_instruction_from_input(input=filter_values[_ligcomp_instruction], varname=f'{_ligcomp}:{_ligcomp_instruction}')
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_ligcomp}:{_denticities}', allow_none=True)
            elif self.filtername == _metals_of_interest:
                out_filter_settings[_metals_of_interest] = self.get_list_of_chemical_elements_from_input(input=filter_values[_metals_of_interest], varname=f'{_metals_of_interest}')
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_metals_of_interest}:{_denticities}', allow_none=True)
            elif self.filtername == _ligand_charges:
                out_filter_settings[_ligand_charges] = self.get_list_of_ints_from_input(input=filter_values[_ligand_charges], varname=f'{_ligand_charges}', allow_none=True)
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_ligand_charges}:{_denticities}', allow_none=True)
            elif self.filtername == _coords:
                out_filter_settings[_coords_atoms_of_interest] = self.get_list_of_chemical_elements_from_input(input=filter_values[_coords_atoms_of_interest], varname=f'{_coords}:{_coords_atoms_of_interest}')
                out_filter_settings[_coords_instruction] = self.get_instruction_from_input(input=filter_values[_coords_instruction], varname=f'{_coords}:{_coords_instruction}')
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_coords}:{_denticities}', allow_none=True)
            elif self.filtername == _remove_ligands_with_neighboring_coordinating_atoms:
                out_filter_settings[_remove_ligands_with_neighboring_coordinating_atoms] = self.get_bool_from_input(
                input=filter_values[_remove_ligands_with_neighboring_coordinating_atoms],
                varname=_remove_ligands_with_neighboring_coordinating_atoms)
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_remove_ligands_with_neighboring_coordinating_atoms}:{_denticities}', allow_none=True)
            elif self.filtername == _remove_ligands_with_beta_hydrogens:
                out_filter_settings[_remove_ligands_with_beta_hydrogens] = self.get_bool_from_input(
                    input=filter_values[_remove_ligands_with_beta_hydrogens],
                    varname=_remove_ligands_with_beta_hydrogens)
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_remove_ligands_with_beta_hydrogens}:{_denticities}', allow_none=True)
            elif self.filtername == _mw:
                out_filter_settings.update(self.check_min_max_input(filter_values=filter_values, filter_name=_mw))
            elif self.filtername == _acount:
                out_filter_settings.update(self.check_min_max_input(filter_values=filter_values, filter_name=_acount))
            elif self.filtername == _interatomic_distances:
                out_filter_settings.update(self.check_min_max_input(filter_values=filter_values, filter_name=_interatomic_distances))
            elif self.filtername == _planarity:
                out_filter_settings.update(self.check_min_max_input(filter_values=filter_values, filter_name=_planarity))
            elif self.filtername == _occurrences:
                out_filter_settings.update(self.check_min_max_input(filter_values=filter_values, filter_name=_occurrences))
            elif self.filtername == _remove_missing_bond_orders:
                out_filter_settings[_remove_missing_bond_orders] = self.get_bool_from_input(input=filter_values[_remove_missing_bond_orders], varname=_remove_missing_bond_orders)
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_remove_missing_bond_orders}:{_denticities}', allow_none=True)
            elif self.filtername == _atm_neighbors:
                out_filter_settings[_atom] = self.get_list_of_chemical_elements_from_input(input=filter_values[_atom], varname=f'{_atm_neighbors}:{_atom}')[0]
                out_filter_settings[_neighbors] = self.get_list_of_chemical_elements_from_input(input=filter_values[_neighbors], varname=f'{_atm_neighbors}:{_neighbors}')
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_atm_neighbors}:{_denticities}', allow_none=True)
            elif self.filtername == _smarts_filter:
                out_filter_settings[_smarts] = filter_values[_smarts]
                out_filter_settings[_should_be_present] = self.get_bool_from_input(input=filter_values[_should_be_present], varname=_should_be_present)
                out_filter_settings[_include_metal] = self.get_bool_from_input(input=filter_values[_include_metal], varname=_include_metal)
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_smarts_filter}:{_denticities}', allow_none=True)
            else:
                self.raise_error(f"Filter '{self.filtername}' is not a valid filter.", varname=_filter)

            out_settings.append(out_filter_settings)

        return out_settings

    def get_output_ligand_db_path(self):
        """
        Returns default path to the output ligand database.
        """
        path = Path('filtered_ligand_db.jsonlines').resolve()

        return path

    def check_min_max_input(self, filter_values: dict, filter_name: str) -> dict:
        outsettings = {
            _min: self.get_float_from_input(input=filter_values[_min], varname=f'{filter_name}:{_min}', allow_none=True),
            _max: self.get_float_from_input(input=filter_values[_max], varname=f'{filter_name}:{_max}', allow_none=True),
            _denticities: self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{filter_name}:{_denticities}', allow_none=True)
            }

        return outsettings

    def check_stoichiometry_input(self, stoichiometry: str) -> str:
        """
        Checks the stoichiometry input. Accepts the stoichiometry in format where if an element is present with a quantity of 1, the 1 can optionally be omitted.
        @return: Stoichiometry in format where if an element is present with a count of 1, the 1 is not omitted.
        """
        stoichiometry = str(stoichiometry)
        stoichiometry = Composition(stoichiometry).get_stoichiometry()

        return stoichiometry



    def check_ligand_db_path(self, path) -> Path:
        path = get_correct_ligand_db_path_from_input(path)
        if path is None:
            self.raise_error(f'Invalid ligand database filepath.', varname=_ligand_db_path)
        self.ensure_file_present(path=path, varname=_ligand_db_path)

        return Path(path)

    def check_denticities_of_interest(self, settings) -> Union[list, None]:
        input = settings[_denticities_of_interest]
        if isinstance(input, (tuple, list)):
            input = self.get_list_of_ints_from_input(input=input, varname=f'denticity in {_denticities_of_interest}')

        return input

    def check_even_odd_electron_count_input(self, settings) -> Union[str, None]:
        input = settings[_filter_even_odd_electron_count]
        if isinstance(input, str):
            input = input.lower()

        if not (input is None or input in ['even', 'odd']):
            self.raise_error(message=f"Input '{input}' can not be recognized as (None, even, odd).", varname=_filter_even_odd_electron_count)

        return input




class AssemblyInput(BaseInput):
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
                        _batches: [list, tuple, dict],
                        _complex_name_length: [int, str],
                        _same_isomer_names: [bool, str],
                        }
    # Batch settings
    batches_valid_keys = {
                        _name: [str],
                        _input_path: [str, list, tuple, type(None)],
                        _max_num_complexes: [int, str],
                        _element: [str],
                        _oxidation_state: [int, str],
                        _isomers: [str],
                        _optimisation: [str, bool],
                        _random_seed: [int, str],
                        _total_charge: [int, str],
                        _topology: [str, list, tuple],
                        _complex_name_appendix: [str, type(None)],
                        _geometry_modifier_filepath: [str, type(None)],
                        _bidentate_rotator: [str],
                        }
    total_keys = deepcopy(valid_keys)
    total_keys.update({
        _batches: batches_valid_keys,
        })

    def __init__(self, path: Union[str, Path] = 'assembly_input.yml'):
        """
        Reads the global settings file and stores the settings in the class.
        """
        # Read into dict
        super().__init__(path)
        self.global_settings = self.get_settings_from_input_file(path=self.path)

        # Set the batch name to None. It will be set later to the respective batch name when iterating over the batches so that the error messages can be more specific.
        self.batch_name = None

        # Check the input and set the class variables
        self.verbose = None
        self.optimization_movie = None
        self.concatenate_xyz = None
        self.Output_Path = None
        self.Batches = None
        self.complex_name_length = None
        self.check_and_set_global_settings()
        self.check_if_settings_not_recognized(actual_settings=self.global_settings, valid_settings=self.total_keys)
        self.check_batches_input()

    def check_and_set_global_settings(self):
        """
        Checks the global settings and sets them as attributes.
        """
        for key, types in self.valid_keys.items():
            if key not in self.global_settings.keys():
                if key not in assembler_gbl_defaults.keys():
                    self.raise_error(f"Mandatory key '{key}' not found in input file. Please add it.")
                else:
                    self.global_settings[key] = assembler_gbl_defaults[key]

            self.check_correct_input_type(input=self.global_settings[key], types=types, varname=key)

        self.verbose = self.get_int_from_input(self.global_settings[_verbose], varname=_verbose, allow=allowed_verbosities)
        self.optimization_movie = self.get_bool_from_input(self.global_settings[_optimization_movie], varname=_optimization_movie)
        self.concatenate_xyz = self.get_bool_from_input(self.global_settings[_concatenate_xyz], varname=_concatenate_xyz)
        self.Batches =  self.get_batches_from_input(self.global_settings[_batches])
        self.complex_name_length = self.get_int_from_input(self.global_settings[_complex_name_length], varname=_complex_name_length)
        self.same_isomer_names = self.get_bool_from_input(self.global_settings[_same_isomer_names], varname=_same_isomer_names)

        # Check if the output path exists and is a directory.
        self.Output_Path = self.ensure_output_directory_valid(self.global_settings[_output_path], varname=_output_path)

        return

    def get_ligandfilters_from_input(self, ligandfilters_path: Union[str, Path, None]) -> Union[None, dict]:
        """
        Reads the ligand filters from the ligand filters input file.
        """
        if ligandfilters_path is None:
            return None

        ligandfilters = LigandFilterInput(path=ligandfilters_path).get

        return ligandfilters

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
                if key in assembler_batch_defaults.keys():
                    batch_settings[key] = assembler_batch_defaults[key]
                else:
                    self.raise_error(f"Mandatory key '{key}' not found in input file. Please add it.")

            varname = f"{_batches}->{key}"
            self.check_correct_input_type(input=batch_settings[key], types=types, varname=varname)

        # Here we take the batch inputs and format them correctly
        Ligand_json, topology_similarity = self.get_ligand_db_path_and_topologies_from_input(batch_settings[_input_path], topology=batch_settings[_topology])
        Max_Num_Assembled_Complexes = self.get_int_from_input(batch_settings[_max_num_complexes], varname=f'{_batches}->{_max_num_complexes}', allow_str='all')
        Generate_Isomer_Instruction = self.get_isomers_from_input(batch_settings[_isomers])
        Optimisation_Instruction = self.get_bool_from_input(batch_settings[_optimisation], varname=f'{_batches}->{_optimisation}')
        Random_Seed = self.get_int_from_input(batch_settings[_random_seed], varname=f'{_batches}->{_random_seed}')
        Total_Charge = self.get_int_from_input(batch_settings[_total_charge], varname=f'{_batches}->{_total_charge}')
        metal_list = self.get_metal_from_input(element=batch_settings[_element], oxidation_state=batch_settings[_oxidation_state])
        complex_name_appendix = batch_settings[_complex_name_appendix] or ''
        geometry_modifier_filepath = self.get_geometry_modifier_from_input(batch_settings[_geometry_modifier_filepath])
        bidentate_rotator = self.get_bidentate_rotator_from_input(batch_settings[_bidentate_rotator], varname=f'{_batches}->{_bidentate_rotator}')

        if isinstance(Ligand_json, list):
            similarity = topology_similarity.split('--')[1].lstrip('[').rstrip(']').split(', ')
            n_diff_ligands = len(set(similarity))
            if not len(Ligand_json) == n_diff_ligands:
                self.raise_error(f"Input '{_input_path}' is a list of paths and must have the same length as the number of different ligands specified in the similarity list at the end of the topology. Yet, the topology {topology_similarity} specifies {n_diff_ligands} different ligands, but {len(Ligand_json)} paths were given.", varname=f'{_batches}->{_input_path}')

        return self.batch_name, Ligand_json, Max_Num_Assembled_Complexes, Generate_Isomer_Instruction, Optimisation_Instruction, Random_Seed, Total_Charge, metal_list, topology_similarity, complex_name_appendix, geometry_modifier_filepath, bidentate_rotator

    def get_bidentate_rotator_from_input(self, bidentate_rotator: str, varname: str):
        """
        Checks the input for the bidentate rotator.
        """
        if not bidentate_rotator in ['auto', 'slab', 'horseshoe']:
            self.raise_error(f"Input '{bidentate_rotator}' for '{varname}' is not valid. Valid inputs are 'auto', 'slab' and 'horseshoe'.")

        return bidentate_rotator

    def get_random_seed_from_input(self, random_seed: Union[int, str, None], varname: str) -> Union[int, None]:
        """
        Checks the input for the random seed.
        """
        if random_seed is None:
            random_seed = get_random_random_seed()
        else:
            self.get_int_from_input(random_seed, varname=varname)

        return random_seed

    def get_geometry_modifier_from_input(self, path):
        if path is None:
            return None

        path = self.ensure_file_present(path, varname=f'{_batches}->{_geometry_modifier_filepath}')

        try:
            mols = ase.io.read(path, index=':', format='xyz')
        except Exception as e:
            self.raise_error(f"Error reading geometry modifier file '{path}': {e}")

        if not len(mols) == 2:
            self.raise_error(f"Geometry modifier file '{path}' must contain exactly two concatenated geometries, but contains {len(mols)} geometries.")

        old_geometry, new_geometry = mols
        if len(old_geometry) != len(new_geometry):
            self.raise_error(f"Geometry modifier file '{path}' must contain two geometries with the same number of atoms, but the two geometries have {len(old_geometry)} and {len(new_geometry)} atoms, respectively.")

        old_atoms = list(old_geometry.get_chemical_symbols())
        new_atoms = list(new_geometry.get_chemical_symbols())
        if not old_atoms == new_atoms:
            self.raise_error(f"Geometry modifier file '{path}' must contain two geometries with the same elements in the same order, but the elements differ: {old_atoms} and {new_atoms}.")

        return path

    def get_ligand_db_path_and_topologies_from_input(self, ligand_db_path, topology) -> Tuple[Union[Path, list], str]:
        """
        Checks the input for the ligand database path and the topology.
        @param: ligand_db_path: Path to the ligand database. If None, the default ligand database will be used. If a list, must be of same length as the topology.
        """
        _, topology_list = self.get_topology_from_input(topology)
        n_ligands = len(topology_list)
        varname = f'{_batches}->{_input_path}'

        # Allow specifying list with length 1 instead of single path
        if isinstance(ligand_db_path, (list, tuple)) and len(ligand_db_path) == 1:
            assert n_ligands != 1
            ligand_db_path = ligand_db_path[0]

        if ligand_db_path is None or isinstance(ligand_db_path, (str, Path)):   # Single input ligand database path is given
            output_ligand_db_path = get_correct_ligand_db_path_from_input(ligand_db_path)
            if output_ligand_db_path is None:
                self.raise_error(f"Invalid ligand database filepath.", varname=varname)
            similarity_list = list(range(1, n_ligands + 1))

        elif isinstance(ligand_db_path, (list, tuple)):   # List of paths and keywords is given
            # If the input is a list, it can be a list of either ligand db paths or same_ligand_keywords. Here we map the ligand db input to the historical way, in which we have a similarity list specifying same ligands and the ligand db list has the reduced length of how many different similarities there are.
            if len(ligand_db_path) != n_ligands: # new format: ligand db list must be same length as topology
                self.raise_error(f"Inconsistent number of denticities and ligand databases. Input is a list of ligand db paths of length {len(ligand_db_path)} while the number of specified denticities in the topology is {n_ligands}. Please provide input where the number of ligand database paths is the same as the number of denticities in the topology.")

            # Reduce ligand db list to old format, in which there are no keywords and the number of elements is the number of different similarities
            output_ligand_db_path = []
            for input_path in ligand_db_path:
                if input_path == _same_ligand_keyword:
                    continue
                path = get_correct_ligand_db_path_from_input(input_path)
                if path is None:
                    self.raise_error(f"Invalid ligand database filepath.", varname=varname)
                path = Path(self.ensure_file_present(path, varname=varname))
                output_ligand_db_path.append(path)

            # The first entry cannot be a keyword since the keywords always reference to the previous ligand db path.
            if ligand_db_path[0] == _same_ligand_keyword:
                self.raise_error(f"Invalid first entry '{_same_ligand_keyword}' in ligand database path list. Please provide a list where the first entry is a ligand database path. Only later elements can be this keyword since it always refers to the previous element in the list.")

            # Build the similarity list such that it maps to the historic format, i.e. [1, 2, 2] for one different and two same ligands
            similarity = 0
            similarity_list = []
            for path in  ligand_db_path:
                if path != _same_ligand_keyword:
                    similarity += 1
                similarity_list.append(similarity)

        full_topology_str = str(topology_list) + '--' + str(similarity_list)
        return output_ligand_db_path, full_topology_str

    def get_topology_from_input(self, topology: str) -> tuple[str, list]:
        """
        Checks the topology input for correct input.
        @topology: Input topology in the format '3-2-1'.
        @returns: Tuple(topology, topology_list): A tuple of a string of the format '[3, 2, 1]' for specifying denticities and the same thing as list.
        """
        varname = f'{_batches}->{_topology}'

        topology = str(topology)
        if not topology in allowed_topologies:
            self.raise_error(f"Invalid topology '{topology}'. Supported topologies are {allowed_topologies}.", varname=varname)

        # Remove 'mer-' from the topology string for backwards compatibility with old code
        if topology.startswith('mer-'):
            topology = topology[4:]

        denticities = [int(dent) for dent in topology.split('-')]

        # error_message = f"Topology '{topology}' is not in the correct format. It must be a list of denticities in the format '3-2-1'."
        # try:
        #     denticities = ast.literal_eval(topology)
        # except (ValueError, SyntaxError):
        #     # If the input is weird, then we raise an error
        #     self.raise_error(error_message, varname=varname)
        #
        # # Check that denticities are either lists or tuples and make them to lists for the rest of the code
        # if not isinstance(denticities, (list, tuple)):
        #     self.raise_error(error_message, varname=varname)
        # denticities = list(self.get_list_of_ints_from_input(input=denticities, varname=varname))
        #
        # if not any(sorted(denticities) == sorted(top) for top in allowed_topologies):
        #     self.raise_error(f"Invalid topology '{topology}'. This topology is not supported. Supported topologies are {allowed_topologies}.", varname=varname)

        # Validity checks
        # - Check that denticities are positive integers
        if not all(isinstance(denticity, int) and denticity > 0 for denticity in denticities):
            self.raise_error(f"Invalid topology '{topology}'. Please provide a list of positive integers.", varname=varname)
        # - Check that same integers are clustered together in the list.
        occurrences = pd.Series(denticities).value_counts()
        for dent, occ in occurrences.items():
            dent_indices = [i for i, d in enumerate(denticities) if d == dent]
            min_index, max_index = min(dent_indices), max(dent_indices)
            if not max_index - min_index + 1 == occ:
                self.raise_error(f"Invalid topology '{topology}'. Please provide a list of integers where same integers are clustered together. For example, the topology 2-2-1 is good, while 2-1-2 is bad because not all '2' appear in series.", varname=varname)

        output_topology = str(denticities)
        return output_topology, denticities

    def get_metal_from_input(self, element: str, oxidation_state: int) -> list:
        """
        Checks the metal input for correct input.
        """
        # Check the element
        if not is_chemical_element(element):
            self.raise_error(f"Input element '{element}' is not a valid chemical element symbol, e.g. 'Fe' for iron.", varname=_element)
        element = str(element)

        # Check the oxidation state
        oxidation_state = self.get_int_from_input(oxidation_state, varname=_oxidation_state)
        if oxidation_state > 0:
            oxidation_state = f"+{oxidation_state}"
        else:
            self.raise_error(f"Input oxidation state '{oxidation_state}' is not a positive integer > 0.", varname=_oxidation_state)

        # Check the spin
        # spin = self.get_int_from_input(spin, varname=_spin)
        # if spin < 0:
        #     self.raise_error(f"Input spin '{spin}' is not a positive integer >= 0.", varname=f'{varname}->{_spin}')
        # spin = str(spin)

        return [element, oxidation_state]

    def get_isomers_from_input(self, isomers):
        """
        Checks if the isomers input is correct.
        """
        values = {'lowest_energy': 'Generate Lowest Energy', 'all': 'Generate All'}

        possible_values = list(values.keys())
        isomers = str(isomers)
        if isomers not in possible_values:
            self.raise_error(f" Input value '{isomers}' not recognized. It must be one of {possible_values} (case sensitive).", varname=f'{_batches}->{_isomers}')

        isomers = values[isomers]   # Convert input format to the format used in the assembly code

        return isomers



def is_chemical_element(element: str) -> bool:
    """
    Checks if the input is a valid chemical element.
    """
    return element in all_atomic_symbols

def get_dict_tree_as_string(d: dict, sep: str = ':') -> list[str]:
    """
    Returns a list of all keys in the dict where the item is a string of the dictpath to the element in the format 'key1:key2:key3'
    """
    dict_tree = []
    for key, value in d.items():
        if isinstance(value, dict):
            dict_tree += [f'{key}{sep}{subkey}' for subkey in get_dict_tree_as_string(d=value, sep=sep)]
        else:
            dict_tree += [f'{key}']
    return dict_tree

def find_element_in_dict_from_key_path(d: dict, key_path: str, sep: str = ':') -> Any:
    keys = key_path.split(sep)
    rv = d
    for key in keys:
        rv = rv[key]
    return rv

def get_closest_word(word: str, words: Union[list,tuple]) -> str:
    """
    Returns the closest word in the list of words.
    """
    try:
        closest_word = difflib.get_close_matches(word, words, n=1)[0]
    except:
        closest_word = ''

    return closest_word

class AssemblyInputError(Exception):
    """
    Exception raised for errors in the input.
    """
    def __init__(self, message: str, varname: str='', file: str='', batch_name:str =''):
        if varname != '':
            varname = f" for key '{varname}'"
        file = Path(file).name
        if file != '':
            file = f" in input file '{file}'"
        total_message = f"\n\t--> Assembly Input Error{varname}{batch_name}{file}:\n\t\t{message}"
        super().__init__(total_message)

class LigandCombinationError(Exception):
    """
    Exception raised for errors in the choice of ligands.
    """
    def __init__(self, message: str, file: str='', batch_name:str =''):
        file = Path(file).name
        if file != '':
            file = f" in input file '{file}'"
        total_message = f"\n\t--> Ligand Combination Error{batch_name}{file}:\n\t\t{message}"
        super().__init__(total_message)

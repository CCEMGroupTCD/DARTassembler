"""
This file contains functions and a classes for the input of the assembly. They doublecheck that the input is correct and convert it to the correct format.
"""
import ast
import difflib
import warnings
from copy import deepcopy
import ase

from DARTassembler.src.constants.Paths import default_ligand_db_path
from DARTassembler.src.constants.Periodic_Table import all_atomic_symbols
from pathlib import Path
from typing import Union, Any
from DARTassembler.src.ligand_extraction.io_custom import read_yaml

# Define the key names in the assembly input file
# Global settings
_verbose = 'verbose'
_optimization_movie = 'optimization_movie'
_concatenate_xyz = 'concatenate_xyz'
_overwrite_output_path = 'overwrite_output_path'
_output_path = 'Output_Path'
_complex_name_length = 'complex_name_length'
# Batch settings
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
_complex_name_appendix = 'complex_name_appendix'
_geometry_modifier_filepath = 'geometry_modifier_filepath'
_bidentate_rotator = 'bidentate_rotator'
# _ligand_filters_path = 'ligand_filters_path'


# Define names for the settings in the ligand filters file:
# Global
_ligand_db_path = 'input_ligand_db_path'
_output_ligand_db_path = 'output_ligand_db_path'
_filters = 'filters'
# Filters
_filter = 'filter'
_denticities = 'denticities'
_denticities_of_interest = 'denticities_of_interest'
_remove_ligands_with_neighboring_coordinating_atoms = 'remove_ligands_with_neighboring_coordinating_atoms'
_only_confident_charges = 'only_confident_charges'
_remove_ligands_with_beta_hydrogens = 'remove_ligands_with_beta_hydrogens'
_strict_box_filter= 'strict_box_filter'
_filter_even_odd_electron_count = 'remove_even_odd_electron_count'
_dentfilters = 'denticity_dependent_filters'
# Denticity dependent filters
_acount = 'atom_count'
_acount_min = 'min'
_acount_max = 'max'

_ligcomp = 'ligand_composition'
_ligcomp_atoms_of_interest = 'atoms_of_interest'
_ligcomp_instruction = 'instruction'

_ligand_charges = 'ligand_charges'

_metals_of_interest = 'metals_of_interest'

_coords = 'coordinating_atoms_composition'
_coords_atoms_of_interest = 'atoms_of_interest'
_coords_instruction = 'instruction'

_mw = 'molecular_weight'
_mw_min = 'min'
_mw_max = 'max'


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

    def get_int_from_input(self, input: Union[str, int], varname: str, allow_none=False) -> Union[int,None]:
        """
        Returns an int from a string or int input.
        """
        if allow_none and input is None:
            return None

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

        if isinstance(input, str):
            input = [input]

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
        for key, types in valid_keys.items():
            real_keys = tuple(settings.keys())
            if key not in real_keys:
                similar_word = get_closest_word(word=key, words=real_keys)
                similar_string = f"The closest key provided is '{similar_word}'. " if similar_word != '' else ''
                # Check if key is present in input file
                self.raise_error(f"Key '{key}' not found in input file, please add it. {similar_string}All keys found are {real_keys}.")
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
                self.raise_warning(message=f"Setting '{actual_key}' is not recognized and will be skipped.", varname=actual_key)

        for batch in actual_settings[_batches]:
            for actual_key in batch:
                if not actual_key in valid_settings[_batches]:
                    self.raise_warning(message=f"Setting '{actual_key}' is not recognized and will be skipped.",
                                       varname=actual_key)

                if actual_key == _metal:
                    for metal_key in batch[_metal]:
                        if not metal_key in valid_settings[_metal]:
                            self.raise_warning(message=f"Setting '{metal_key}' is not recognized and will be skipped.",
                                               varname=actual_key)

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
        total_message = f"Invalid input{varname}{batch_name}{file}: {message}"
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
        _filters: [list, tuple],
    }

    filter_keys = {
        _denticities_of_interest: {
            _denticities_of_interest: [list, tuple]},
        _remove_ligands_with_neighboring_coordinating_atoms: {
            _remove_ligands_with_neighboring_coordinating_atoms: [bool, str]},
        # _only_confident_charges: {                    # filter removed from options and made mandatory
        #     _only_confident_charges: [bool, str]},
        _remove_ligands_with_beta_hydrogens: {
            _remove_ligands_with_beta_hydrogens: [bool, str]},
        _strict_box_filter: {
            _strict_box_filter: [bool, str]},
        # _filter_even_odd_electron_count: {            # filter removed from options because all ligands with confident charges have even electron count
        #     _filter_even_odd_electron_count: [str]},
        _acount: {
            _acount_min: [int, str, type(None)],
            _acount_max: [int, str, type(None)],
            _denticities: [list, tuple, type(None)],
            },
        _ligcomp: {
            _ligcomp_atoms_of_interest: [list, tuple],
            _ligcomp_instruction: [str],
            _denticities: [list, tuple, type(None)],
            },
        _ligand_charges: {
            _ligand_charges: [list, tuple],
            _denticities: [list, tuple, type(None)],
            },
        _metals_of_interest: {
            _metals_of_interest: [list, tuple],
            _denticities: [list, tuple, type(None)],
            },
        _coords: {
            _coords_atoms_of_interest: [list, tuple],
            _coords_instruction: [str],
            _denticities: [list, tuple, type(None)],
            },
        _mw: {
            _mw_min: [float, str, type(None)],
            _mw_max: [float, str, type(None)],
            _denticities: [list, tuple, type(None)],
            },
        }


    def __init__(self, path: Union[str, Path]):
        """
        Class for reading and checking the input file for the ligand filter. The input file should be a yaml file.
        todo: Add that missing filter keys are added with value None so that one doesn't have to specify everything.
        """
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

        self.check_input_types(valid_keys=self.valid_keys, settings=settings)
        self.ligand_db_path = self.check_ligand_db_path(settings)
        self.output_ligand_db_path = self.get_path_from_input(path=settings[_output_ligand_db_path], varname=_output_ligand_db_path, allow_none=True)
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

            filter_values = {key: value for key, value in full_filter.items() if key != _filter}
            self.check_input_types(valid_keys=self.filter_keys[self.filtername], settings=filter_values)

            out_filter_settings = {_filter: self.filtername}
            if self.filtername == _denticities_of_interest:
                out_filter_settings[_denticities_of_interest] = self.check_denticities_of_interest(settings=filter_values)
            elif self.filtername == _remove_ligands_with_neighboring_coordinating_atoms:
                out_filter_settings[_remove_ligands_with_neighboring_coordinating_atoms] = self.get_bool_from_input(input=filter_values[_remove_ligands_with_neighboring_coordinating_atoms], varname=_remove_ligands_with_neighboring_coordinating_atoms)
            elif self.filtername == _only_confident_charges:
                out_filter_settings[_only_confident_charges] = self.get_bool_from_input(input=filter_values[_only_confident_charges], varname=_only_confident_charges)
            elif self.filtername == _remove_ligands_with_beta_hydrogens:
                out_filter_settings[_remove_ligands_with_beta_hydrogens] = self.get_bool_from_input(input=filter_values[_remove_ligands_with_beta_hydrogens], varname=_remove_ligands_with_beta_hydrogens)
            elif self.filtername == _strict_box_filter:
                out_filter_settings[_strict_box_filter] = self.get_bool_from_input(input=filter_values[_strict_box_filter], varname=_strict_box_filter)
            elif self.filtername == _filter_even_odd_electron_count:
                out_filter_settings[_filter_even_odd_electron_count] = self.check_even_odd_electron_count_input(settings=filter_values)
            # Denticity dependent filters
            elif self.filtername == _acount:
                out_filter_settings[_acount_min] = self.get_int_from_input(input=filter_values[_acount_min], varname=f'{_acount}:{_acount_min}', allow_none=True)
                out_filter_settings[_acount_max] = self.get_int_from_input(input=filter_values[_acount_max], varname=f'{_acount}:{_acount_max}', allow_none=True)
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_acount}:{_denticities}', allow_none=True)
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
            elif self.filtername == _mw:
                out_filter_settings[_mw_min] = self.get_float_from_input(input=filter_values[_mw_min], varname=f'{_mw}:{_mw_min}', allow_none=True)
                out_filter_settings[_mw_max] = self.get_float_from_input(input=filter_values[_mw_max], varname=f'{_mw}:{_mw_max}', allow_none=True)
                out_filter_settings[_denticities] = self.get_list_of_ints_from_input(input=filter_values[_denticities], varname=f'{_mw}:{_denticities}', allow_none=True)

            out_settings.append(out_filter_settings)

        return out_settings

    def get_output_ligand_db_path(self):
        """
        Returns a path to the output ligand database. If the path already exists, it will be renamed to avoid overwriting.
        """
        path = Path('filtered_ligand_db.json').resolve()

        idx = 1
        while path.exists():
            path = Path(f'{path}_{idx}')
            idx += 1

        return path

    def check_ligand_db_path(self, settings) -> Path:
        path = settings[_ligand_db_path]
        if path is None and default_ligand_db_path.exists():
            path = default_ligand_db_path
        else:
            self.raise_error(f'Input ligand db path is not specified and no db found at default location. Please specify a valid path to the input ligand db.', varname=_ligand_db_path)
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
                        _overwrite_output_path: [bool, str],
                        _complex_name_length: [int, str],
                        }
    # Batch settings
    batches_valid_keys = {
                        _name: [str],
                        _input_path: [str, list, tuple, type(None)],
                        _max_num_complexes: [int, str],
                        _isomers: [str],
                        _optimisation: [str, bool],
                        _random_seed: [int, str],
                        _total_charge: [int, str],
                        _topology: [str, list, tuple],
                        _metal: [dict],
                        _complex_name_appendix: [str,type(None)],
                        _geometry_modifier_filepath: [str, type(None)],
                        _bidentate_rotator: [str],
                        # _ligand_filters_path: [str, Path, type(None)],
                        }
    # Metal settings in the batch settings
    metal_valid_keys = {
                        _element: [str],
                        _oxidation_state: [int, str],
                        _spin: [int, str],
                        }
    total_keys = deepcopy(valid_keys)
    total_keys.update({
        _batches: batches_valid_keys,
        _metal: metal_valid_keys,
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
        self.overwrite_output_path = None
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
                self.raise_error(f"Key '{key}' not found in input file. Please add it.")

            self.check_correct_input_type(input=self.global_settings[key], types=types, varname=key)

        self.verbose = self.get_int_from_input(self.global_settings[_verbose], varname=_verbose)
        self.overwrite_output_path = self.get_bool_from_input(self.global_settings[_overwrite_output_path], varname=_overwrite_output_path)
        self.optimization_movie = self.get_bool_from_input(self.global_settings[_optimization_movie], varname=_optimization_movie)
        self.concatenate_xyz = self.get_bool_from_input(self.global_settings[_concatenate_xyz], varname=_concatenate_xyz)
        self.Batches =  self.get_batches_from_input(self.global_settings[_batches])
        self.complex_name_length = self.get_int_from_input(self.global_settings[_complex_name_length], varname=_complex_name_length)

        # Check if the output path exists and is a directory. Must come after self.overwrite_output_path is set.
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
        complex_name_appendix = batch_settings[_complex_name_appendix] or ''
        geometry_modifier_filepath = self.get_geometry_modifier_from_input(batch_settings[_geometry_modifier_filepath])
        bidentate_rotator = self.get_bidentate_rotator_from_input(batch_settings[_bidentate_rotator], varname=f'{_batches}->{_bidentate_rotator}')

        if isinstance(Ligand_json, list):
            similarity = topology_similarity.split('--')[1].lstrip('[').rstrip(']').split(', ')
            n_diff_ligands = len(set(similarity))
            if not len(Ligand_json) == n_diff_ligands:
                self.raise_error(f"Input '{_input_path}' is a list of paths and must have the same length as the number of different denticities specified in the topology. Yet, the topology {topology_similarity} has {n_diff_ligands} different denticities, but {len(Ligand_json)} paths were given.", varname=f'{_batches}->{_input_path}')

        return self.batch_name, Ligand_json, Max_Num_Assembled_Complexes, Generate_Isomer_Instruction, Optimisation_Instruction, Random_Seed, Total_Charge, metal_list, topology_similarity, complex_name_appendix, geometry_modifier_filepath, bidentate_rotator

    def get_bidentate_rotator_from_input(self, bidentate_rotator: str, varname: str):
        """
        Checks the input for the bidentate rotator.
        """
        if not bidentate_rotator in ['auto', 'slab', 'horseshoe']:
            self.raise_error(f"Input '{bidentate_rotator}' for '{varname}' is not valid. Valid inputs are 'auto', 'slab' and 'horseshoe'.")

        return bidentate_rotator

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

    def get_ligand_db_path_from_input(self, ligand_db_path) -> Union[Path, list]:
        """
        Checks the input for the ligand database path.
        @param: ligand_db_path: Path to the ligand database. If None, the default ligand database will be used. If a list, must be of same length as the topology.
        """
        varname = f'{_batches}->{_input_path}'
        if ligand_db_path is None or ligand_db_path == '' or ligand_db_path == 'default':
            print(f"Ligand database path '{varname}' is not specified in batch {self.batch_name} in file '{self.path}'. Using full default ligand database.")
            ligand_db_path = Path(default_ligand_db_path)
        else:
            # Check if the file exists
            if isinstance(ligand_db_path, str):
                ligand_db_path = Path(self.ensure_file_present(ligand_db_path, varname=varname))
            elif isinstance(ligand_db_path, (list, tuple)):
                ligand_db_path = [Path(self.ensure_file_present(path, varname=varname)) for path in ligand_db_path]
            else:
                self.raise_error(f"Input '{varname}' must be a string or a list/tuple of strings, but is {type(ligand_db_path)}.")

        return ligand_db_path

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
        @returns: topology: A string of the format '[3, 2, 1]' or '[3, 2, 1]--[3, 2, 1]' for specifying denticities and similarities, respectively.
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
        total_message = f"\n\t--> Invalid input{varname}{batch_name}{file}:\n\t\t{message}"
        super().__init__(total_message)

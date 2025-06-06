from datetime import datetime
import numpy as np
from ase import Atoms
import warnings
import pandas as pd
from collections import namedtuple
from typing import Tuple
from DARTassembler.src.constants.chem import Element

def sort_dict_as(d: dict, order: list[str], strict: bool=False) -> dict:
    """
    Sorts a dictionary by the order of keys specified in the `order` list. If `strict` is True, it raises an error if the keys of the dictionary do not match the order.
    :param d: Dictionary to sort
    :param order: List of keys in the desired order
    :param strict: If True, raises an error if the keys of the dictionary do not match the order. If False, ignores keys that are not in the order.
    :return: Sorted dictionary
    """
    if strict:
        diff_set = set(order).symmetric_difference(set(d.keys()))
        if len(diff_set) > 0:
            raise ValueError(f'Dictionary keys {d.keys()} do not match the order {order}. Difference: {diff_set}')

    real_keys = list(d.keys())
    sorted_keys = sorted(real_keys, key=lambda x: order.index(x) if x in order else len(order))
    sorted_dict = {key: d[key] for key in sorted_keys}
    assert set(real_keys) == set(
        sorted_dict.keys()), f'Sorting failed. Real keys: {real_keys}, sorted keys: {sorted_dict.keys()}'

    return sorted_dict

def get_stable_sorted_value_counts(values: list, ascending: bool=False) -> dict:
    """
    Returns a dictionary of value counts sorted by the value in ascending order. The sorting is stable, i.e. the order of equal values is preserved by ordering first by the values and second by the alphabetical order of the keys.
    :param values: list of values
    :param ascending: Whether to sort the indices in ascending order (True) or descending order (False)
    :return: Dictionary of value counts sorted by the value in ascending order
    """
    # Using mergesort is important to make this sorting stable
    return pd.Series(values).value_counts().sort_index().sort_values(ascending=ascending, kind='mergesort').to_dict()

def angle_between_ab_ac_vectors(a: list, b: list, c: list, degrees: bool=True) -> float:
    """
    Calculate the angle at point 'a' formed by line segments 'ab' and 'ac'.
    """
    a, b, c = np.array(a), np.array(b), np.array(c)
    if a.shape != b.shape or b.shape != c.shape:
        raise ValueError("All points must have the same dimensions.")

    ab = b - a
    ac = c - a
    cosine_angle = np.dot(ab, ac) / (np.linalg.norm(ab) * np.linalg.norm(ac))
    cosine_angle = np.clip(cosine_angle, -1.0, 1.0)  # clip for numerical stability
    angle = np.arccos(cosine_angle)

    if degrees:
        angle = np.degrees(angle)
    return angle

def angles_of_triangle(point1: list, point2: list, point3: list, degrees: bool=True) -> Tuple[float, float, float]:
    """
    Calculate the angles of a triangle formed by three points.
    """
    angle1 = np.degrees(angle_between_ab_ac_vectors(point1, point2, point3, degrees))
    angle2 = np.degrees(angle_between_ab_ac_vectors(point2, point1, point3, degrees))
    angle3 = np.degrees(angle_between_ab_ac_vectors(point3, point1, point2, degrees))

    return angle1, angle2, angle3

def series2namedtuple(s, name='S'):
    return namedtuple(name, s.index)(*s)

def unroll_dict_into_columns(df, dict_col: str, prefix: str, delete_dict: bool=False) -> pd.DataFrame:
    dict_data = [d for d in df[dict_col]]
    df_dict_data = pd.DataFrame(dict_data, index=df.index)
    df_dict_data = df_dict_data.rename(columns={col: prefix + col for col in df_dict_data})

    df = df.join(df_dict_data, validate='1:1')
    if delete_dict:
        df = df.drop(columns=dict_col)

    return df

def make_None_to_NaN(val):
    if val is None:
        return np.nan
    else:
        return val

def update_dict_with_warning_inplace(dict_to_update, dict_with_information, update_properties: list=None):
    if update_properties is None:
        update_properties = list(dict_with_information.keys())

    for prop in update_properties:
        different_values =  prop in dict_to_update and \
                            dict_with_information[prop] != dict_to_update[prop]
        if different_values:
            warnings.warn(f'Overwriting dictionary with property {prop} which already existed.')
        dict_to_update[prop] = dict_with_information[prop]

    return

def sort_dict_recursively_inplace(d: dict) -> None:
    """
    Sorts all dictionaries in a dictionary of dictionaries recursively by the key.
    """
    d = {key: d[key] for key in sorted(d.keys())}
    for value in d.values():
        if isinstance(value, dict):
            sort_dict_recursively_inplace(value)

    return

def sorted_dict_of_dicts(d: dict, debug=False) -> dict:

    """
    Sorts dictionaries and recursively sorts dictionaries of dictionaries to infinite order.
    :param d: dictionary to sort
    :return: sorted dictionary
    """
    sorted_d = {}
    keys = sorted(d.keys())

    for key in keys:
        value = d[key]

        if (isinstance(value, dict) and (len(value) > 1)):
            value = sorted_dict_of_dicts(value)

        sorted_d[key] = value

    if debug:
        assert (len(d) == len(sorted_d) and all([val == d[key] for key, val in
                                             sorted_d.items()])), 'Sorted dictionary is different than original one, there must be a bug.'
    return sorted_d

def call_method_on_object(obj, method: str):
    return getattr(obj, method)

def is_between(value: float, range: list, include_left=True, include_right=True) -> bool:
    if len(range) == 0:
        return True
    if not sorted(range) == list(range):
        raise ValueError(f'`range` is not sorted!: {range}')

    elif len(range) == 2:
        if include_left and include_right:
            between = range[0] <= value <= range[1]
        elif include_left:
            between = range[0] <= value < range[1]
        elif include_right:
            between = range[0] < value <= range[1]
        else:
            between = range[0] < value < range[1]
    else:
        ValueError(f'Length of `range` must be either 0 or 2 but is {len(range)}.')

    return between

def identify_metal_in_ase_mol(mol: Atoms):
    metals = [el for el in mol.get_chemical_symbols() if Element(el).is_metal]
    assert len(metals) == 1, f"Molecule has more than one metal: {metals}"
    return metals[0]

def identify_metal_in_atoms_list(atoms: list):
    metals = [el for el in atoms if Element(el).is_metal]
    assert len(metals) == 1, "Molecule seems to be not a single metal complex, metal identification failed"

    return metals[0]

def coordinates_to_xyz_str(coordinates: dict):
    str_ = f"{len(list(coordinates.keys()))}\n \n"
    for coord in coordinates.values():
        str_ += f"{coord[0]} \t {coord[1][0]} \t {coord[1][1]} \t {coord[1][2]} \n"

    return str_

def atomic_props_dict_to_lists(prop: dict, flatten=False) -> list:
    """
    Converts an atomic property (e.g. partial charge, coordinates) from format {0: ['La', [1, 2, 3]], ...} into three lists of indices, atoms, values.
    :param prop: atomic property (e.g. partial charge, coordinates) in format {0: ['La', [1, 2, 3]], ...}
    :return: lists of indices, atoms, values
    """
    indices = []
    atoms = []
    values = []
    for idx, l in prop.items():
        indices.append(idx)
        atoms.append(l[0])
        values.append(l[1])
    if not flatten:
        return indices, atoms, values
    else:
        values = np.array(values).T.tolist()
        return indices, atoms, *values

def original_xyz_indices_to_indices_wo_metal(orig_coordinates: dict) -> list:
    """
    Converts the indices of atoms of the original xyz files to the new indices when the metal atom is deleted from the xyz. That means all atoms after the metal shift one index up.
    :param orig_coordinates: Coordinates of original xyz file
    :return: Dictionary mapping original xyz indices to new wo_metal indices
    """
    orig_indices = [(idx, l[0]) for idx, l in orig_coordinates.items()]
    
    orig_to_wo_metal_indices = {}
    counter = 0
    for orig_idx, el in orig_indices:
        
        is_metal = Element(el).is_metal
        if is_metal:
            metal_idx = orig_idx
            continue
        
        orig_to_wo_metal_indices[orig_idx] = counter
        counter += 1
    
    # Double checking.
    assert len(orig_to_wo_metal_indices) == len(orig_indices) - 1
    for orig_idx, new_idx in orig_to_wo_metal_indices.items():
        if orig_idx < metal_idx:
            assert orig_idx == new_idx
        elif orig_idx == metal_idx:
            assert not metal_idx in orig_to_wo_metal_indices
        else:
            assert orig_idx == new_idx + 1
    
    return orig_to_wo_metal_indices

def convert_atomic_props_from_original_xyz_indices_to_indices_wo_metal(atomic_props, orig_coords):
    """
    Changes the indices of an atomic property dictionary in order to match the new indices when the metal is deleted in the xyz.
    """
    orig_to_wo_metal_indices = original_xyz_indices_to_indices_wo_metal(orig_coords)
    
    atomic_props_new_idc = {}
    for orig_idx, prop in atomic_props.items():
        if orig_idx in orig_to_wo_metal_indices:
            new_idx = orig_to_wo_metal_indices[orig_idx]
            atomic_props_new_idc[new_idx] = prop
    
    all_new_elements = [props[0] for props in atomic_props_new_idc.values()]
    assert not any([Element(el).is_metal for el in all_new_elements]), 'Found metal in ligand? Or Index Error.'
    return atomic_props_new_idc

def get_all_atomic_properties_in_long_array(atoms: list, coords: list, atomic_props: dict):
    all_atomic_props = np.array(coords).T.tolist()  # transpose list
    for name, prop in atomic_props.items():
        _, _, prop_values = atomic_props_dict_to_lists(prop, flatten=True)
        all_atomic_props.append(prop_values)
    all_atomic_props.append(atoms)
    all_atomic_props = np.array(all_atomic_props, dtype='object').T
    
    return all_atomic_props

def get_all_atomic_properties_with_modified_coordinates_wo_metal_in_long_array(atoms, coords, atomic_props, modified_coordinates):
    not_metal_idx = [not Element(el).is_metal for el in atoms]
    
    all_atomic_props = get_all_atomic_properties_in_long_array(atoms, coords, atomic_props)
    all_atomic_props_wo_metal = all_atomic_props[not_metal_idx, :]
    
    _, wo_metal_atoms, *wo_metal_coords = atomic_props_dict_to_lists(modified_coordinates, flatten=True)
    # replace original coordinates with modified coordinates
    all_atomic_props_wo_metal[:, 0:3] = np.array(wo_metal_coords).T
    
    assert wo_metal_atoms == all_atomic_props_wo_metal[:, 4].tolist()
    return all_atomic_props_wo_metal

def flatten_list(l: list) -> list:
    return [item for sublist in l for item in sublist]

def get_duration_string(start: datetime, without_microseconds=True) -> str:
    end = datetime.now()
    duration = end - start

    if without_microseconds:
        duration = str(duration).split('.')[0]
    else:
        duration = str(duration)

    return duration

def check_equal(val1, val2) -> bool:
    """
    Checks if two values are equal. Respects None, np.nan and empty arrays. None and np.nan are compared as equal.
    :return: bool
    """
    # Handle the case that either value is None or np.nan. Use pd.isna() so that both None and np.nan are treated as equal.
    nan_val1 = pd.isna(val1)
    nan_val2 = pd.isna(val2)
    if nan_val1 and nan_val2:
        return True
    if nan_val1 or nan_val2:
        return False

    # Handle the case that both values are empty arrays. They are equal.
    if isinstance(val1, np.ndarray) and isinstance(val2, np.ndarray):
        if len(val1) == 0 and len(val2) == 0:
            return True
        return np.array_equal(val1, val2)

    # Handle the case of floats. Use np.isclose to compare floats.
    if isinstance(val1, float) and isinstance(val2, float):
        return np.isclose(val1, val2)

    return val1 == val2
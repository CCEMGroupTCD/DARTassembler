import numpy as np
from pymatgen.core.periodic_table import Element as Pymatgen_Element
from src01.constants import metals_in_pse
from ase import Atoms
import networkx as nx
import warnings
from copy import deepcopy
import pandas as pd

def unroll_dict_into_columns(df, dict_col: str, prefix: str):
    dict_data = [d for d in df[dict_col]]
    df_dict_data = pd.DataFrame(dict_data, index=df.index)
    df_dict_data = df_dict_data.rename(columns={col: prefix + col for col in df_dict_data})
    df = df.join(df_dict_data, validate='1:1')

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

def sorted_dict_of_dicts(d: dict) -> dict:
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

    assert (len(d) == len(sorted_d) and all([val == d[key] for key, val in
                                             sorted_d.items()])), 'Sorted dictionary is different than original one, there must be a bug.'
    return sorted_d


def call_method_on_object(obj, method: str):
    return getattr(obj, method)


def identify_metal_in_ase_mol(mol: Atoms):

    metals = set(mol.get_atomic_numbers()).intersection(set(metals_in_pse))
    assert len(metals) == 1, "Molecule seems to be not a single metal complex, metal identification failed"

    return Pymatgen_Element.from_Z(metals.pop()).symbol


def coordinates_to_xyz_str(coordinates: dict):
    """
    returns a string that can be written into an .xyz file
    """
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
        
        is_metal = Pymatgen_Element(el).is_metal
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
    :param atomic_props:
    :param orig_coords:
    :return:
    """
    orig_to_wo_metal_indices = original_xyz_indices_to_indices_wo_metal(orig_coords)
    
    atomic_props_new_idc = {}
    for orig_idx, prop in atomic_props.items():
        if orig_idx in orig_to_wo_metal_indices:
            new_idx = orig_to_wo_metal_indices[orig_idx]
            atomic_props_new_idc[new_idx] = prop
    
    all_new_elements = [props[0] for props in atomic_props_new_idc.values()]
    assert not any([Pymatgen_Element(el).is_metal for el in all_new_elements]), 'Found metal in ligand? Or Index Error.'
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
    not_metal_idx = [not Pymatgen_Element(el).is_metal for el in atoms]
    
    all_atomic_props = get_all_atomic_properties_in_long_array(atoms, coords, atomic_props)
    all_atomic_props_wo_metal = all_atomic_props[not_metal_idx, :]
    
    _, wo_metal_atoms, *wo_metal_coords = atomic_props_dict_to_lists(modified_coordinates, flatten=True)
    # replace original coordinates with modified coordinates
    all_atomic_props_wo_metal[:, 0:3] = np.array(wo_metal_coords).T
    
    assert wo_metal_atoms == all_atomic_props_wo_metal[:, 4].tolist()
    return all_atomic_props_wo_metal


def flatten_list(l: list) -> list:
    return [item for sublist in l for item in sublist]
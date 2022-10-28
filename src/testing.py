"""
All kinds of utitily functions specifically for testing.
"""
from src.utilities import atomic_props_dict_to_lists
import numpy as np
from mendeleev import element

def all_same(*inputs: tuple) -> bool:
    inp0 = inputs[0]
    all_inputs_same = True
    for inp in inputs:
        if inp != inp0:
            all_inputs_same = False
            break
            
    return all_inputs_same

def all_same_number_of_elements(*inputs: tuple) -> bool:
    lenghts = [len(i) for i in inputs]
    all_same_length = len(set(lenghts)) == 1
    
    return all_same_length


def check_partial_charges_add_up(mol_id, partial_charges: dict, total_charge: int, rtol=1e-4, atol: float=1e-4):
    pc_idc, pc_atoms, pc_values = atomic_props_dict_to_lists(partial_charges,
                                                             flatten=True)
    assert np.isclose(sum(pc_values), total_charge, rtol=rtol, atol=atol), f'{mol_id}: Partial charges added up are not equal to total charge.'
    
    return


def check_if_metal_entries_consistent(global_props: dict, mol_id: str, other_metal_entries: list):
    metal_entries = other_metal_entries
    
    if 'metal_str_if_exists' in global_props:
        metal_oxi_str = global_props['metal_str_if_exists']
        if isinstance(metal_oxi_str, str):  # skip if NaN
            metal_name = metal_oxi_str.split('(')[0]
            metal_entries.append(metal_name)
    
    norm_metal_entries = [element(metal.capitalize()).symbol if isinstance(metal, str) else element(int(metal)).symbol for metal in
                          metal_entries]
    assert all_same(norm_metal_entries), f'{mol_id}: Inconsistent metals: {metal_entries}.'
    
    return
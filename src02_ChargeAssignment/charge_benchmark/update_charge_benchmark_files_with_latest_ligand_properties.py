"""
Note: The functionality of this script has been integrated in the script `merge_benchmark_charges.py`. Therefore, it is currently not useful for anything but kept just in case.

This script updates the charge benchmark datasets with the ligand properties of the latest complex_db version.
"""
import pandas as pd
from src01.io_custom import load_full_ligand_db
from pathlib import Path
from copy import deepcopy
from collections import MutableMapping
from contextlib import suppress
import numpy as np

def get_value_of_nested_dict(dic, keys):
    for key in keys:
        dic = dic[key]
    return dic

def make_nested_dict_inplace(dic, keys, value):
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value

    return

def delete_keys_from_nested_dict(dictionary, keys):
    for key in keys:
        with suppress(KeyError):
            del dictionary[key]
    for value in dictionary.values():
        if isinstance(value, MutableMapping):
            delete_keys_from_nested_dict(value, keys)

def set_keys_of_dict(d: dict, id_cols: list):
    new_d = {}
    duplicate_keys = []
    for values in d.values():
        new_values = {key: val for key, val in values.items() if not key in id_cols}
        keys = [values[id_col] for id_col in id_cols]

        try:
            # Check if keys already exist in dict and if so, add to list
            get_value_of_nested_dict(dic=new_d, keys=keys)
            duplicate_keys.append(keys)
        except KeyError:
            make_nested_dict_inplace(dic=new_d, keys=keys, value=new_values)

    return new_d

def get_properties_of_benchmark_dataset_from_full_ligand(full_ligand: dict, old_ligand:dict, update_columns: list) -> dict:
    old_ligand = deepcopy(old_ligand)
    full_ligand_props = list(full_ligand.keys())

    assert all([col in full_ligand_props for col in update_columns]), f'Some of the properties in `update_columns` ({update_columns}) can not be found in the ligand properties ({full_ligand_props}).'
    for prop in update_columns:
        old_ligand[prop] = full_ligand[prop]

    return old_ligand

if __name__ == '__main__':

    db_version = '1.3'

    latest_full_ligand_db_path = f'../../data/tmQMG_Jsons/tmQM_Ligands_full_v{db_version}.json'
    charge_benchmark_datasets = [f'../../database/ligand_charges/charge_benchmark/all_ligand_charges_with_high_confidence_v{db_version}.csv']
    id_cols = ['CSD_code', 'stoichiometry']
    update_columns = ['unique_name', 'name', 'graph_hash']

    assert latest_full_ligand_db_path.endswith(f'_v{db_version}.json'), f'Specified version number {db_version} does\'t match with version number in `latest_complex_db_path`.'
    full_ligand_db = load_full_ligand_db(latest_full_ligand_db_path)
    full_ligand_dict = set_keys_of_dict(d=full_ligand_db, id_cols=id_cols)

    for charge_benchmark_path in charge_benchmark_datasets:
        df = pd.read_csv(charge_benchmark_path)
        df_old = deepcopy(df)

        ligand_props = [prop for prop in df.columns if not prop in id_cols]

        benchmark_ligands = df.to_dict(orient='index')
        for lig in benchmark_ligands.values():

            # If the CSD code or the stoichiometry cannot be found in the ligand database, write NaN to the benchmark csv.
            try:
                latest_lig = full_ligand_dict[lig['CSD_code']][lig['stoichiometry']]
                new_lig_props = {key: latest_lig[key] for key in update_columns}
            except KeyError:
                new_lig_props = {col: np.nan for col in update_columns}

            lig.update(new_lig_props)

        df = pd.DataFrame.from_dict(benchmark_ligands, orient='index')

        pass



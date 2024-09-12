
import pandas as pd
from copy import deepcopy
import gc
import json
from tqdm import tqdm
import numpy as np
from pathlib import Path

from DARTassembler.src.constants.Paths import project_path
from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand
from DARTassembler.src.ligand_extraction.utilities_Molecule import get_standardized_stoichiometry_from_atoms_list
from DARTassembler.src.ligand_extraction.io_custom import load_complex_db, load_full_ligand_db, load_unique_ligand_db, save_unique_ligand_db, save_complex_db, save_full_ligand_db
from DARTassembler.src.linear_charge_solver.linear_charge_solver.linear_charge_solver import LinearChargeSolver
from typing import Union

def get_charges_of_unique_ligands(all_complexes: Union[str, Path, dict], max_iterations=None) -> pd.DataFrame:
    """
    So far uses only the linear charge solver.
    :param all_complexes: path to a json of all complexes with ligands with unique ligand names
    :return: dataframe with the charge of each unique ligand
    """
    save_dir = str(project_path().extend('data', 'linear_charge_solver', 'output'))

    solver = LinearChargeSolver(
                                all_complexes=all_complexes,
                                save_dir=save_dir,
                                max_iterations=max_iterations
                                )
    df_ligands = solver.calculate_unique_ligand_charges(output_uncertain_charges_as_nan=False)

    return df_ligands


def update_ligand_with_charge_inplace(lig: RCA_Ligand, charges: dict):
    new_keys = ['pred_charge', 'pred_charge_confidence', 'pred_charge_is_confident']
    assert not any([hasattr(lig, key) for key in new_keys])

    uname = lig.unique_name
    if uname in charges:
        charge = charges[uname]['pred_charge']
        confidence = charges[uname]['confidence']
        is_confident = charges[uname]['is_confident']
        pred_charge_exact = charges[uname]['pred_charge_exact']
    else:
        # This means that this unique ligand was not part of the output of the charge assignment because its complexes were filtered out in the charge assignment due to lack of oxidation state or similar issues.
        charge = np.nan
        confidence = np.nan
        is_confident = False
        pred_charge_exact = np.nan

    lig.pred_charge = charge
    lig.pred_charge_is_confident = is_confident
    lig.global_props.update({
                                'LCS_pred_charge': charge,
                                'LCS_pred_charge_confidence': confidence,
                                'LCS_pred_charge_is_confident': is_confident,
                                'LCS_pred_charge_exact': pred_charge_exact
                                })

    return


def update_databases_with_charges(df_ligand_charges: pd.DataFrame, data_store_path: Union[str, Path]):
    complex_db_path = Path(data_store_path, 'complex_db.json')
    ligands_db_path = Path(data_store_path, 'tmQM_Ligands_full.json')
    unique_ligands_db_path = Path(data_store_path, 'tmQM_Ligands_unique.json')

    df_ligand_charges = df_ligand_charges.set_index('unique_name')
    charges = df_ligand_charges.to_dict(orient='index')

    unique_ligands = load_unique_ligand_db(path=unique_ligands_db_path)
    not_intersecting_ligands = set(unique_ligands.keys()).symmetric_difference(set(charges.keys()))
    print(f'Number of unique ligands not outputted by the charge assignment: {len(not_intersecting_ligands)}')
    for ulig in unique_ligands.values():
        update_ligand_with_charge_inplace(ulig, charges)
    save_unique_ligand_db(db=unique_ligands, path=unique_ligands_db_path)

    ligands = load_full_ligand_db(ligands_db_path)
    for lig in ligands.values():
        update_ligand_with_charge_inplace(lig, charges)
    save_full_ligand_db(db=ligands, path=ligands_db_path)

    complexes = load_complex_db(path=complex_db_path)
    for c in complexes.values():
        for lig in c['ligands']:
            update_ligand_with_charge_inplace(lig, charges)
    save_complex_db(db=complexes, path=complex_db_path)

    return


def unique_ligands_from_Ligand_batch_json_files(n=10,
                                                data_store_path: str = str(project_path().extend('data', 'tmQMG_Jsons_test'))
                                                ):
    gc.collect()

    # first we generate the full ligand dict
    ligand_dict = {}
    # expected_length = 0
    for i in range(n):
        with open(f"{data_store_path}/Lig_Batch_{i}.json", "r") as handle:
            dict_ = json.load(handle)
            ligand_dict.update(dict_)
            # expected_length += len(dict_)

    #
    name_ghash_dict = {}
    # we append the dict piecewise by the graph hash
    for lig_key, lig_dict in tqdm(ligand_dict.items(), desc="Extracting graph hashs"):
        # new_ligand_dict[lig_key] = lig_dict["graph_hash"]
        name_ghash_dict[lig_key] = lig_dict["graph_hash"]

    grouped_ligands_by_hash = {}
    for k, v in tqdm(name_ghash_dict.items(), desc="grouping ligands by hash"):
        grouped_ligands_by_hash.setdefault(v, []).append(k)

    unique_ligand_dict = {}
    for same_ligands in tqdm(grouped_ligands_by_hash.values(), desc="filling unique ligand dict"):
        name = same_ligands[0]
        uname = 'unq_' + name
        for ligand_name in same_ligands:
            ligand_dict[ligand_name]["unique_name"] = uname

        unique_ligand = deepcopy(ligand_dict[name])

        denticities = [ligand_dict[ligand_name]['denticity'] for ligand_name in same_ligands]
        metals = [ligand_dict[ligand_name]['original_metal_symbol'] for ligand_name in same_ligands]

        # Add useful statistical information of all ligands for this unique ligand
        count_denticities = pd.Series(denticities).value_counts().sort_values(ascending=False).to_dict()
        count_metals = pd.Series(metals).value_counts().sort_values(ascending=False).to_dict()
        unique_ligand_infos = {
                                'occurrences': len(same_ligands),
                                'count_denticities': count_denticities,
                                'count_metals': count_metals,
                                'n_denticities': len(count_denticities),
                                'n_metals': len(count_metals),
        }
        unique_ligand.update(unique_ligand_infos)
        unique_ligand['all_ligand_names'] = same_ligands

        for ligand_name in same_ligands:
            # update ligands with unique_ligand information for easier debugging
            ligand_dict[ligand_name]['unique_ligand_information'] = unique_ligand_infos
            # update ligands with bool flag if this ligand is the ligand chosen as unique ligand
            ligand_dict[ligand_name]['is_chosen_unique_ligand'] = ligand_name == name

        # Delete attribute original metal from unique_ligand since it is confusing and no real attribute of a unique ligand
        del unique_ligand['original_metal']
        del unique_ligand['original_metal_symbol']

        unique_ligand_dict[uname] = unique_ligand

    #
    # Now we can dump ligand_dict to json, which is a huge dict containing ALL ligands with their unique name
    # and also unique_ligand_dict which only contains unique ligands
    with open(f"{data_store_path}/tmQM_Ligands_full.json", "w+") as file:
        json.dump(ligand_dict, file)
        print(f"Full ligand database saved to {f'{data_store_path}/tmQM_Ligands_full.json'}")
    with open(f"{data_store_path}/tmQM_Ligands_unique.json", "w+") as file:
        json.dump(unique_ligand_dict, file)
        print(f"Unique ligand database saved to {f'{data_store_path}/tmQM_Ligands_unique.json'}")

    return


def update_complex_db_with_ligands(complex_json: str, ligand_json: str, save_complex_db_path: Union[str, Path]):
    """
    This function reads in the initial json file that was input of the ligand extraction process and also reads in the json with all extracted ligands. The dictionary for each complex is then updated with dictionaries of the ligands and again saved.
    This function is just very much hacked together and will be subject of refactoring everything into a single class soon.
    :param complex_json: path to the initial input json of the ligand extraction process
    :param ligand_json: path to the json of all extracted ligands
    :param save_complex_db_path: path to where the with ligands enriched complex db should be saved as json
    :return:
    """
    print('Start updating complex db with ligands.')

    CSD_global_props_paths = ['../Bin_(Old)/databases(raw)/tmQM_raw/raw_data/CSD.csv', '../database/tmQM/raw_data/CSD.csv']
    CSD_global_props_path = None
    for path in CSD_global_props_paths:
        if Path(path).exists():
            CSD_global_props_path = path
    if CSD_global_props_path is None:
        raise ValueError(f'CSD.csv could not be found at any of the given paths {CSD_global_props_paths}.')
    CSD_global_props = pd.read_csv(CSD_global_props_path, index_col=0).set_index('CSD_code').to_dict(orient='index')

    with open(complex_json, 'r') as file:
        all_complexes = json.load(file)
    with open(ligand_json, 'r') as file:
        all_ligands_list = json.load(file)

    all_csd_codes = np.unique([lig['CSD_code'] for lig in all_ligands_list.values()]).tolist()
    all_ligands = {c_id: {} for c_id in all_csd_codes}
    for lig_id, lig in all_ligands_list.items():
        c_id = lig['CSD_code']
        all_ligands[c_id][lig_id] = lig

    # Delete superfluous complexes which had been filtered out but are still in the input complexes.
    all_complex_csd_codes = list(all_complexes.keys())
    deleted_csd_codes = [c_id for c_id in all_complex_csd_codes if not c_id in all_csd_codes]
    for c_id in deleted_csd_codes:
        del all_complexes[c_id]

    for c_id, c in tqdm(all_complexes.items(), desc='Updating complex db'):
        ligands = all_ligands[c_id]

        if c_id in CSD_global_props:
            csd_global = CSD_global_props[c_id]
        else:
            csd_global = {'name': np.nan, 'metal_nr_if_exists': np.nan, 'metal_name': np.nan}

        del c['graph_dict']

        metals = [lig['original_metal_symbol'] for lig in ligands.values()]
        assert len(np.unique(metals)) == 1, 'different original metals for ligands from the same complex.'
        c['mol_id'] = c_id
        c['stoichiometry'] = get_standardized_stoichiometry_from_atoms_list(c['atomic_props']['atoms'])
        c['metal'] = metals[0]
        c['metal_oxi_state'] = csd_global['metal_nr_if_exists']
        c['total_q'] = c['global_props']['charge']
        c['Metal_q'] = c['atomic_props']['partial_charge'][c['atomic_props']['atoms'].index(c['metal'])]
        c['global_props']['CSD_iupac_name'] = csd_global['name']

        c['ligands'] = []
        for lig_id, lig in ligands.items():
            del lig['graph_dict']
            c['ligands'].append(lig)

        assert c['global_props']['CSD_code'] == c_id
        assert c['total_q'] == c['global_props']['charge']

        assert np.isclose(sum(c['atomic_props']['partial_charge']), c['total_q'], atol=1e-3), 'Formal and partial charges are not consistent.'
        if not isinstance(csd_global['metal_name'], float):
            assert csd_global['metal_name'] == c['metal']

    with open(save_complex_db_path, 'w') as file:
        json.dump(all_complexes, file)

    print(f'Saved complex db with ligand information to {save_complex_db_path}.')





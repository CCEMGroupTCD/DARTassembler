import pandas as pd
import json
from collections import Counter, defaultdict
import numpy as np
from pathlib import Path
from copy import deepcopy
from tqdm import tqdm

from pymatgen.core.composition import Composition
from itertools import combinations
def find_min_diff(arr):
    arr = sorted(arr)
    diff = 10 ** 20     # inf

    # Find the min diff by comparing adjacent
    # pairs in sorted array
    for i in range(len(arr) - 1):
        if arr[i + 1] - arr[i] < diff:
            diff = arr[i + 1] - arr[i]

    # Return min diff
    return diff

def get_stoichiometry_from_ligand(ligand: dict) -> str:
    """
    Returns stoichiometry sorted by
    :param ligand:
    :return:
    """
    if 'stoichiometry' in ligand:
        return ligand['stoichiometry']

    atomic_props = ligand['atomic_props']
    if 'coordinates' in atomic_props:
        elements = atomic_props['coordinates']
    elif 'partial_charge' in atomic_props:
        elements = atomic_props['partial_charge']
    else:
        raise ValueError('Cannot calculate stoichiometry from given ligand.')

    stoichiometry = Counter([val[0] for val in elements.values()])
    stoichiometry = ''.join(sorted([f'{el}{count}' for el, count in stoichiometry.items()]))

    return stoichiometry

def similar_but_not_same(lig, ligands, unique_id):
    for lig2 in ligands:
        same_stoich = lig['stoichiometry'] == lig2['stoichiometry']
        same_ligands = lig[unique_id] == lig2[unique_id]
        if same_stoich and not same_ligands:
            return True

    return False



if __name__ == '__main__':
    all_complexes_path = '../../../data/final_db_versions/complex_db_1.5_all.json'
    # all_complexes_path = '../../../data_output/CSD_MM_G_Jsons_test/complex_db_10k_1_final.json' # 10k complexes for testing
    batches = ['Cian', 'Marconi', 'Manting']
    n_ligands_per_batch = 500
    save_csvs_dir = '../../../test/debug/databases/charge_benchmark/CSD_ligands_to_benchmark'



    print(f'Start loading complex db from {all_complexes_path}.')
    with open(all_complexes_path, 'r') as file:
        all_complexes = json.load(file)
    print(f'Loaded all complexes.')

    all_same_ligands = defaultdict(list)
    for c_id, c in all_complexes.items():
        for lig in c['ligands']:
            all_same_ligands[lig['unique_name']].append(lig['name'])

    all_ligands = []
    for c_id, c in tqdm(all_complexes.items(), desc='Build full ligand db'):
        # Skip complexes without metal oxidation state since we cannot calculate their charges anyway.
        if np.isnan(c['metal_oxi_state']):
            continue

        for lig in c['ligands']:
            stoichiometry = get_stoichiometry_from_ligand(lig)
            same_ligands = [name for name in all_same_ligands[lig['unique_name']] if not name == lig['name']]
            all_ligands.append({
                                'CSD_code': c_id,
                                'metal': c['metal'],
                                'stoichiometry': stoichiometry,
                                'occurrences': sum([lig['unique_name'] == lig2['unique_name'] for lig2 in c['ligands']]),
                                'connected': lig['denticity'] > 0,
                                'charge': np.nan,
                                'denticity': np.nan,
                                'issue_detected': np.nan,
                                'comment': np.nan,
                                'confidence': np.nan,
                                'same_ligands': same_ligands,
                                'unique_name': lig['unique_name'],
                                'orig_d': lig['denticity'],
                                'too_similar': similar_but_not_same(lig=lig, ligands=c['ligands'], unique_id='unique_name'),
                                # 'metal_oxi_state': c['metal_oxi_state']
                                })

    df = pd.DataFrame(all_ligands)
    df_all = deepcopy(df)

    # Drop complexes with multiple ligands with the same stoichiometry so that each ligand can be identified by the complex and the stoichiometry. If the ligands are exactly the same it doesn't need to be dropped because it will have the same charge anyway.
    df = df[~df['too_similar']]
    df = df.drop_duplicates(subset=['unique_name'], keep='first')        # must be last filter
    df = df.sample(frac=1, random_state=1).reset_index(drop=True)

    # for debugging
    df_benchmarked = df.head(n_ligands_per_batch*len(batches))
    denticities = df_benchmarked['orig_d'].value_counts()

    for i, name in enumerate(batches):
        start_idx = i * n_ligands_per_batch
        end_idx = start_idx + n_ligands_per_batch
        df_batch = df.iloc[start_idx:end_idx,:]
        df_batch = df_batch.drop(columns=['unique_name', 'orig_d', 'too_similar'])

        out_filename = Path(save_csvs_dir, f'{name}_ligand_charges_CSD.csv')
        df_batch.to_csv(out_filename, index=False)
        print(f'Saved benchmark csv to {out_filename}.')

    print('Done!')



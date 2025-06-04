import pandas as pd
import json
from collections import Counter
import numpy as np
from pathlib import Path
from copy import deepcopy

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

already_benchmarked_ligands = ['unq_CSD-OCEFOU-01-a', 'unq_CSD-DIRTED-02-a', 'unq_CSD-WUVPAE-03-a', 'unq_CSD-VASYOF-04-a', 'unq_CSD-XIKROZ-03-a', 'unq_CSD-IWAKAS-02-a', 'unq_CSD-GISHUM-02-a', 'unq_CSD-FUNJEE-04-a', 'unq_CSD-QIKRIL-06-a', 'unq_CSD-KEHLOZ-01-a', 'unq_CSD-QEZCOO-02-a', 'unq_CSD-EGUDAM-01-a', 'unq_CSD-IRUCII-03-a', 'unq_CSD-FEMTUN-02-a', 'unq_CSD-FENDAE-02-a', 'unq_CSD-BAPTIY-04-a', 'unq_CSD-OFUNUZ-02-a', 'unq_CSD-MIPDAS-03-a', 'unq_CSD-XALYUE-02-a', 'unq_CSD-JOQYIA-03-a', 'unq_CSD-GEKDUV-02-a', 'unq_CSD-YIRSIE-02-a', 'unq_CSD-GIDLEM-06-a', 'unq_CSD-BOQRIL-03-a', 'unq_CSD-RAKROL-03-a', 'unq_CSD-GERYIO-02-a', 'unq_CSD-ITELIE-04-a', 'unq_CSD-CAQPAM-04-a', 'unq_CSD-CUTXEU-01-b', 'unq_CSD-UTECAZ-01-a', 'unq_CSD-BAJGOJ-04-a', 'unq_CSD-MUQZUU-08-a', 'unq_CSD-RUGWEW-01-b', 'unq_CSD-SATJIJ-01-a', 'unq_CSD-ZIGVES-03-a', 'unq_CSD-TAQTOV-03-a', 'unq_CSD-ZITNOF-04-a', 'unq_CSD-TIQZOJ-01-a', 'unq_CSD-BIYPEH-02-b', 'unq_CSD-YUTKIH-04-a', 'unq_CSD-MUMLAH-02-a', 'unq_CSD-CIGBID-06-a', 'unq_CSD-JUPHIL-02-a', 'unq_CSD-IPTBNI-02-a', 'unq_CSD-BOWSEM-01-a']

output_ligands_in_same_format = ['unq_CSD-OCEFOU-01-a', 'unq_CSD-DIRTED-02-a', 'unq_CSD-WUVPAE-03-a', 'unq_CSD-VASYOF-04-a', 'unq_CSD-XIKROZ-03-a', 'unq_CSD-IWAKAS-02-a', 'unq_CSD-GISHUM-02-a', 'unq_CSD-FUNJEE-04-a', 'unq_CSD-QIKRIL-06-a', 'unq_CSD-KEHLOZ-01-a', 'unq_CSD-QEZCOO-02-a', 'unq_CSD-EGUDAM-01-a', 'unq_CSD-IRUCII-03-a', 'unq_CSD-FEMTUN-02-a', 'unq_CSD-FENDAE-02-a', 'unq_CSD-BAPTIY-04-a', 'unq_CSD-OFUNUZ-02-a', 'unq_CSD-MIPDAS-03-a', 'unq_CSD-XALYUE-02-a', 'unq_CSD-JOQYIA-03-a', 'unq_CSD-GEKDUV-02-a', 'unq_CSD-YIRSIE-02-a', 'unq_CSD-GIDLEM-06-a', 'unq_CSD-BOQRIL-03-a', 'unq_CSD-RAKROL-03-a', 'unq_CSD-GERYIO-02-a', 'unq_CSD-ITELIE-04-a', 'unq_CSD-CAQPAM-04-a', 'unq_CSD-CUTXEU-01-b', 'unq_CSD-UTECAZ-01-a', 'unq_CSD-BAJGOJ-04-a', 'unq_CSD-MUQZUU-08-a', 'unq_CSD-RUGWEW-01-b', 'unq_CSD-SATJIJ-01-a', 'unq_CSD-ZIGVES-03-a', 'unq_CSD-TAQTOV-03-a', 'unq_CSD-ZITNOF-04-a', 'unq_CSD-TIQZOJ-01-a', 'unq_CSD-BIYPEH-02-b', 'unq_CSD-YUTKIH-04-a', 'unq_CSD-MUMLAH-02-a', 'unq_CSD-CIGBID-06-a', 'unq_CSD-JUPHIL-02-a', 'unq_CSD-IPTBNI-02-a', 'unq_CSD-BOWSEM-01-a', 'unq_CSD-LULWUK-02-a', 'unq_CSD-AKODAF-05-a', 'unq_CSD-WAVWAS-01-a', 'unq_CSD-VUSSOT-02-a', 'unq_CSD-SEQJAA-04-a', 'unq_CSD-TEWNEP-02-b', 'unq_CSD-REKWEM-04-a', 'unq_CSD-ENENIW-02-a', 'unq_CSD-VEDGET-04-a', 'unq_CSD-XEVDAE-02-a', 'unq_CSD-QOJZUN-01-a', 'unq_CSD-CUTXEU-01-b', 'unq_CSD-CANLEJ-03-a', 'unq_CSD-FITXAI-06-a', 'unq_CSD-EWIVIQ-04-a', 'unq_CSD-KEZBEZ-08-a', 'unq_CSD-DARTIC-01-a', 'unq_CSD-DEFJOP-02-a', 'unq_CSD-ASOXIQ-02-a', 'unq_CSD-QAQREI-04-a', 'unq_CSD-MEGGUB-03-a', 'unq_CSD-IDUGOE-01-a', 'unq_CSD-BIFWUI-01-a']


if __name__ == '__main__':
    all_complexes_path = '../../data/linear_charge_fitting/all_complexes_tmQMg.json'
    batches = ['Cian', 'Marconi', 'Manting']
    n_ligands_per_batch = 500
    save_csvs_dir = '../../data/linear_charge_fitting'


    with open(all_complexes_path, 'r') as file:
        all_complexes = json.load(file)

    all_ligands = []
    for c_id, c in all_complexes.items():
        for lig in c['ligands']:
            stoichiometry = get_stoichiometry_from_ligand(lig)
            all_ligands.append({
                                'CSD_code': c_id,
                                'stoichiometry': stoichiometry,
                                'metal': c['metal'],
                                'charge': np.nan,
                                'denticity': np.nan,
                                'confidence': np.nan,
                                'comment': np.nan,
                                'unique_name': lig['unique_name'],
                                'name': lig['name'],
                                'graph_hash': lig['graph_hash'],
                                # 'orig_d': lig['denticity']
                            })

    df = pd.DataFrame(all_ligands)
    df_all = deepcopy(df)

    df = df.drop_duplicates(subset=['CSD_code', 'stoichiometry'], keep=False)       # must be first filter
    df = df.drop_duplicates(subset=['unique_name'], keep='first')                   # must be second last filter
    df = df[~df['unique_name'].isin(already_benchmarked_ligands)]                   # must be last filter
    df = df.sample(frac=1, random_state=0).reset_index(drop=True)

    for i, name in enumerate(batches):
        start_idx = i * n_ligands_per_batch
        end_idx = start_idx + n_ligands_per_batch
        df_batch = df.iloc[start_idx:end_idx,:]

        out_filename = Path(save_csvs_dir, f'{name}_ligand_charges.csv')
        df_batch.to_csv(out_filename, index=False)

    # doublechecking
    # df_benchmarked = df.head(n_ligands_per_batch*len(batches))
    # denticities = df_benchmarked['orig_d'].value_counts()

    # output another csv with speficic ligands in the same format for already benchmarked ligands.
    df_spec = df_all[df_all['unique_name'].isin(output_ligands_in_same_format)]
    df_spec = df_spec.drop_duplicates(subset=['CSD_code', 'stoichiometry'], keep=False)
    df_spec = df_spec.drop_duplicates(subset='unique_name')
    df_spec = df_spec.set_index(keys='unique_name')
    df_spec = df_spec.reindex(output_ligands_in_same_format)
    df_spec = df_spec.reset_index()
    out_filename = Path(save_csvs_dir, f'Cian_already_assigned_ligand_charges.csv')
    df_spec.to_csv(out_filename, index=False)



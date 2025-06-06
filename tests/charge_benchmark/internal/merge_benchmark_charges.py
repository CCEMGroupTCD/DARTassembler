import pandas as pd
from pathlib import Path
import numpy as np
from tests.charge_benchmark.internal.charge_benchmark import update_ligands_with_information_from_ligand_db

if __name__ == '__main__':

    benchmark_charge_dir = '../../../../dev/test/debug/databases/charge_benchmark'
    benchmark_charge_filenames = {
            'C1': 'C1.csv',
            'C2': 'C2.csv',
            'C3': 'C3.csv',
            'Man': 'Man.csv',
            'Mar': 'Mar.csv',
        }

    # Update the resulting merged dfs with ligand information from the latest full ligand db.
    # Ligands are identified by 'CSD_code' and 'stoichiometry' (hard coded).
    db_version = '1.4_all'
    latest_full_ligand_db_path = f'../../data/final_db_versions/full_ligand_db_v{db_version}.json'
    update_properties = ['unique_name', 'name', 'graph_hash', 'donor_elements']

    df_all = pd.DataFrame()
    for name, filename in benchmark_charge_filenames.items():
        df0 = pd.read_csv(Path(benchmark_charge_dir, filename))
        df0['author'] = name
        df_all = pd.concat((df_all, df0), axis=0)
    df_all = df_all.reset_index(drop=True)

    df = df_all[df_all['charge'].notna()]
    assert df.drop(columns='comment').notna().all().all(), f'Any of the columns {df.columns} is NaN which should not be NaN'

    df.loc[df['comment'] == '-', 'comment'] = np.nan

    # Add ligand information from the ligand db to each ligand in the charge benchmark csv.
    assert latest_full_ligand_db_path.endswith(f'_v{db_version}.json'), f'Specified version number {db_version} does\'t match with version number in `latest_complex_db_path`.'
    df = df.drop(columns=update_properties, errors='ignore')
    df = update_ligands_with_information_from_ligand_db(
                                                            df_benchmark=df,
                                                            latest_full_ligand_db_path=latest_full_ligand_db_path,
                                                            update_properties=update_properties
                                                        )
    if 'charge' in df:
        df['prediction_error'] = df['charge'] - df['charge']

    df['high_confidence'] = (df['confidence'] == 3) & df['comment'].isna() & df['unique_name'].notna()
    df_confident = df[df['high_confidence']]
    df_confident = df_confident.drop(columns='high_confidence')


    n_duplicates = df.loc[df['graph_hash'].notna(), 'graph_hash'].duplicated().sum()
    if n_duplicates > 0:
        print(f'WARNING: {n_duplicates} duplicates of graph hashes found in the data! They are not excluded as of now.')

    df.to_csv(Path(benchmark_charge_dir, f'all_ligand_charges_v{db_version}.csv'), index=False)
    df_confident.to_csv(Path(benchmark_charge_dir, f'all_ligand_charges_with_high_confidence_v{db_version}.csv'), index=False)

    n_high_confidence = sum(df['high_confidence'])
    print(f'Saved two output csv files, one with all charges with {len(df)} entries and one with only highly confident charges with {n_high_confidence} entries.')

    print('Done!')


    # old_df = pd.read_csv('/Users/timosommer/PhD/projects/RCA/projects/CreateTMC/database/databases/charge_benchmark/all_ligand_charges_new.csv')
    # df = df.reset_index(drop=True)
    # drop_cols = []
    # # not_nan = df['unique_name'].notna()
    # pd.testing.assert_frame_equal(df.drop(columns=drop_cols), old_df[df.columns].drop(columns=drop_cols))
    # print('All good!')





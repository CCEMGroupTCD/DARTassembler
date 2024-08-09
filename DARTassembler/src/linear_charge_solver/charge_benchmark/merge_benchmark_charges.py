import pandas as pd
from pathlib import Path
import numpy as np
from DARTassembler.src.ligand_extraction.io_custom import load_full_ligand_db

def set_dict_keys_to_csd_code_and_stoichiometry(d: dict) -> dict:
    """
    Uses the values of the given dict `d` to setup a new dictionary in the format {CSD_code: {stoichiometry: values}}. This is useful because that way each ligand can be looked up by CSD_code and stoichiometry which has much better performance than searching the whole dictionary for the correct CSD_code and stoichiometry in a loop.
    Additionally, if the original dictionary has duplicates of the keys CSD_code and stoichiometry (which happens frequently) these entries are excluded. This is because these entries should not be used in the charge benchmark anyway and to make sure that these entries are really not used, they are removed.
    """
    new_d = {}
    duplicate_keys = []

    for values in d.values():
        csd_code, stoichiometry, uname = values['global_props']['CSD_code'], values['stoichiometry'], values['unique_name']

        if not csd_code in new_d:
            new_d[csd_code] = {}

        if not stoichiometry in new_d[csd_code]:
            new_d[csd_code][stoichiometry] = values
        else:
            duplicate_keys.append([csd_code, stoichiometry, uname])


    # Remove ligands with not unique CSD_code, stoichiometry, unique name keys.
    duplicate_keys = pd.DataFrame(duplicate_keys).drop_duplicates().values.tolist()
    for dupl_csd_code, dupl_stoichiometry, _ in duplicate_keys:
        try:
            del new_d[dupl_csd_code][dupl_stoichiometry]
        except KeyError:
            pass

    return new_d

def update_ligands_with_information_from_ligand_db(df_benchmark: pd.DataFrame, latest_full_ligand_db_path: str, update_properties: list) -> pd.DataFrame:
    full_ligand_db = load_full_ligand_db(latest_full_ligand_db_path)
    df_full_ligand_db = pd.DataFrame.from_dict(full_ligand_db, orient='index')

    benchmark_ligands = df_benchmark.to_dict(orient='index')
    for lig in benchmark_ligands.values():

        ligand_df = df_full_ligand_db[df_full_ligand_db['global_props'].apply(lambda gbl: gbl['CSD_code'] == lig['CSD_code'])]
        ligand_df = ligand_df[ligand_df['stoichiometry'] == lig['stoichiometry']]

        same_stoi_but_different_ligands = ligand_df['unique_name'].nunique() > 1
        lig_not_in_db = len(ligand_df) == 0

        if same_stoi_but_different_ligands or lig_not_in_db:
            # If the CSD code or the stoichiometry cannot be found in the ligand database, write NaN to the benchmark csv.
            new_lig_props = {col: np.nan for col in update_properties}
        else:
            full_ligand = ligand_df.iloc[0,:].to_dict()

            # Doublecheck that the ligand has all relevant properties given.
            props_not_in_lig = [prop for prop in update_properties if prop not in full_ligand]
            assert not props_not_in_lig, f'Missing properties {props_not_in_lig} in ligand {full_ligand["name"]}.'

            new_lig_props = {key: full_ligand[key] for key in update_properties}

        lig.update(new_lig_props)

    df_benchmark = pd.DataFrame.from_dict(benchmark_ligands, orient='index')

    return df_benchmark



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
    update_properties = ['unique_name', 'name', 'graph_hash', 'local_elements']

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
    if 'pred_charge' in df:
        df['prediction_error'] = df['charge'] - df['pred_charge']

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





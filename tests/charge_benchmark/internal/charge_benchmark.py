from pathlib import Path
import pandas as pd
from typing import Union
import numpy as np
from sklearn.metrics import accuracy_score
from copy import  deepcopy
from DARTassembler.src.constants.paths import project_path, charge_benchmark_dir
from DARTassembler.src.misc.io import load_full_ligand_db


class ChargeBenchmark:

    def __init__(self, true_charge_name: str = 'charge', only_confident: bool = True):
        self.true_charge_name = true_charge_name
        self.only_confident = only_confident # If True, consider only charges where the benchmark charge is confident by the author

        self.benchmark_charge_dir = charge_benchmark_dir
        self.benchmark_charge_filenames = {
            'C1': 'C1.csv',
            'C2': 'C2.csv',
            'Man': 'Man.csv',
            'Mar': 'Mar.csv',
        }
        self.update_properties = ['unique_name', 'name', 'graph_hash', 'donor_elements', 'charge', 'has_confident_charge']

    def merge_benchmark_datasets(self, latest_full_ligand_db_path: str):
        """
        Update the resulting merged dfs with ligand information from the latest full ligand db. Ligands are identified by 'CSD_code' and 'stoichiometry' (hard coded).
        """
        df_all = pd.DataFrame()
        important_columns = ['CSD_code', 'stoichiometry', 'metal', 'charge', 'n_donors', 'issue_detected', 'confidence', 'comment', 'author']
        for name, filename in self.benchmark_charge_filenames.items():
            df0 = pd.read_csv(Path(self.benchmark_charge_dir, filename))
            df0['author'] = name
            df0 = df0[important_columns]
            df_all = pd.concat((df_all, df0), axis=0)
        df_all = df_all.reset_index(drop=True)

        self.df_all = df_all[df_all['charge'].notna()]
        assert self.df_all.drop(
            columns=['issue_detected', 'comment']).notna().all().all(), f'Any of the columns {self.df_all.columns} is NaN which should not be NaN'

        self.df_all.loc[self.df_all['issue_detected'] == '-', 'issue_detected'] = np.nan
        self.df_all.loc[self.df_all['comment'] == '-', 'comment'] = np.nan

        # Add ligand information from the ligand db to each ligand in the charge benchmark csv.
        self.df_all = update_ligands_with_information_from_ligand_db(
            df_benchmark=self.df_all,
            latest_full_ligand_db_path=latest_full_ligand_db_path,
            update_properties=self.update_properties
        )
        if 'charge' in self.df_all:
            self.df_all['prediction_error'] = self.df_all['charge'] - self.df_all['charge']
            # Drop rows where no prediction was made because the ligand was not in the ligand db. This is due to benchmarked ligands which belong to complexes which were filtered out, e.g. because they were missing the metal oxidation state.
            self.df_all = self.df_all[self.df_all['charge'].notna()]

        if self.only_confident:
            self.df_all['high_confidence'] = (self.df_all['confidence'] == 3) & self.df_all['issue_detected'].isna() & self.df_all['unique_name'].notna()
            self.df_confident = self.df_all[self.df_all['high_confidence']]
            self.df_confident = self.df_confident.drop(columns='high_confidence')
        else:
            self.df_confident = deepcopy(self.df_all)


        # n_duplicates = self.df_all.loc[self.df_all['graph_hash'].notna(), 'graph_hash'].duplicated().sum()
        # if n_duplicates > 0:
        #     print(
        #         f'WARNING: {n_duplicates} duplicates of graph hashes found in the data! They are not excluded as of now.')

        return self.df_all

    def calculate_scores_of_charge_benchmark(self, full_ligand_db: Union[str, Path], expected_charge_interval=(-4, 1)):

        full_ligand_db = Path(full_ligand_db)
        self.merge_benchmark_datasets(latest_full_ligand_db_path=full_ligand_db)

        delimiter = '============================'
        print(f'\nSCORES:\n{delimiter}')

        try:
            frac_not_scale_invariant = sum(~self.df_confident['charge_scale_invariant']) / len(self.df_confident)
            print(f'Frac. not scale invariant: {frac_not_scale_invariant:.4g}')
        except KeyError:
            pass

        y_pred = deepcopy(self.df_confident['charge'])
        y_true = deepcopy(self.df_confident[self.true_charge_name])
        valid = y_true.notna() & y_pred.notna()

        y_true_all, y_pred_all = y_true[valid], y_pred[valid]
        acc = accuracy_score(y_true=y_true_all, y_pred=y_pred_all)
        print(f'Total accuracy (n={len(y_true_all)}): {acc:.4g}')

        confident = valid & self.df_confident['has_confident_charge']
        y_true_conf, y_pred_conf = y_true[confident], y_pred[confident]
        acc = accuracy_score(y_true=y_true_conf, y_pred=y_pred_conf)
        print(f'Confident accuracy (n={len(y_true_conf)}): {acc}')

        y_true_no_conf, y_pred_no_conf = y_true[~confident], y_pred[~confident]
        acc = accuracy_score(y_true=y_true_no_conf, y_pred=y_pred_no_conf)
        print(f'Non-confident accuracy (n={len(y_true_no_conf)}): {acc:.4g}')

        if len(confident) > 0:
            print(f'Frac confident predictions: {sum(confident) / len(confident):.4g}')
        else:
            print(f'Frac confident predictions: 0.0')

        print(delimiter + '\n')

        return


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

    ligand_db_version = 'v1.7'
    only_confident = False  # If True, consider only charges where the benchmark charge is confident by the author
    out_csv = Path('..', f'benchmark_ligand_charges_{ligand_db_version}.csv')   # Output file

    charge_benchmark = ChargeBenchmark(only_confident=only_confident)
    complex_db = project_path().extend('data', 'final_db_versions', f'complex_db_{ligand_db_version}.json') # Very big, therefore only local, not on github
    charge_benchmark.calculate_scores_of_charge_benchmark(full_ligand_db=complex_db)

    #%% Write benchmark csv to file
    df_out = charge_benchmark.df_confident
    drop_columns = ['issue_detected', 'author', 'confidence', 'metal', 'n_donors', 'name', 'graph_hash', 'comment']
    int_columns = ['true_charge', 'charge', 'prediction_error']
    rename_columns = {'charge': 'true_charge', 'donor_elements': 'donors'}
    order_columns = ['CSD_code', 'stoichiometry', 'donors', 'true_charge', 'charge', 'has_confident_charge', 'prediction_error']
    index = 'unique_name'
    df_out = df_out.set_index(index)
    df_out = df_out.sort_values(by=['has_confident_charge'], ascending=False)
    df_out = df_out.drop(columns=drop_columns)
    df_out = df_out.rename(columns=rename_columns)
    df_out[int_columns] = df_out[int_columns].astype(int)
    df_out = df_out[order_columns]
    df_out['donors'] = df_out['donors'].apply(lambda x: '-'.join(sorted(x)))
    df_out.to_csv(out_csv, index=True)



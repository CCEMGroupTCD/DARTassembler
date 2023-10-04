from pathlib import Path
import pandas as pd
from typing import Union
import numpy as np
from sklearn.metrics import accuracy_score
from copy import  deepcopy

from src02_ChargeAssignment.charge_benchmark.merge_benchmark_charges import \
    update_ligands_with_information_from_ligand_db


class ChargeBenchmark:

    def __init__(self, true_charge_name: str):
        self.true_charge_name = true_charge_name

        self.benchmark_charge_dir = '../test/debug/databases/charge_benchmark'
        self.benchmark_charge_filenames = {
            'Cian1': 'Cian_already_assigned_ligand_charges.csv',
            'Cian2': 'Cian_BenchMark_corrected_170223.csv',
            'Manting': 'Manting_Corrected_170223.csv',
            'Marconi1': 'Marconi_corrected_by_Cian_020223.csv',
            'Cian3': 'Cian_ligand_charges_CSD_060423.csv'
        }
        self.update_properties = ['unique_name', 'name', 'graph_hash', 'local_elements', 'pred_charge', 'pred_charge_is_confident']

    def merge_benchmark_datasets(self, latest_full_ligand_db_path: str):
        """
        Update the resulting merged dfs with ligand information from the latest full ligand db. Ligands are identified by 'CSD_code' and 'stoichiometry' (hard coded).
        """
        df_all = pd.DataFrame()
        important_columns = ['CSD_code', 'stoichiometry', 'metal', 'charge', 'denticity', 'issue_detected', 'confidence', 'comment', 'author']
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
        if 'pred_charge' in self.df_all:
            self.df_all['prediction_error'] = self.df_all['charge'] - self.df_all['pred_charge']

        self.df_all['high_confidence'] = (self.df_all['confidence'] == 3) & self.df_all['issue_detected'].isna() & self.df_all['unique_name'].notna()
        self.df_confident = self.df_all[self.df_all['high_confidence']]
        self.df_confident = self.df_confident.drop(columns='high_confidence')

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
            print(f'Frac. not scale invariant: {frac_not_scale_invariant:.2g}')
        except KeyError:
            pass

        y_pred = deepcopy(self.df_confident['pred_charge'])
        y_true = deepcopy(self.df_confident[self.true_charge_name])
        valid = y_true.notna() & y_pred.notna()

        y_true_all, y_pred_all = y_true[valid], y_pred[valid]
        acc = accuracy_score(y_true=y_true_all, y_pred=y_pred_all)
        print(f'Total accuracy (n={len(y_true_all)}): {acc:.2g}')

        confident = valid & self.df_confident['pred_charge_is_confident']
        y_true_conf, y_pred_conf = y_true[confident], y_pred[confident]
        acc = accuracy_score(y_true=y_true_conf, y_pred=y_pred_conf)
        print(f'Confident accuracy (n={len(y_true_conf)}): {acc:.2g}')

        y_true_no_conf, y_pred_no_conf = y_true[~confident], y_pred[~confident]
        acc = accuracy_score(y_true=y_true_no_conf, y_pred=y_pred_no_conf)
        print(f'Non-confident accuracy (n={len(y_true_no_conf)}): {acc:.2g}')

        if len(confident) > 0:
            print(f'Frac confident predictions: {sum(confident) / len(confident):.2g}')
        else:
            print(f'Frac confident predictions: 0.0')

        print(delimiter + '\n')

        return

if __name__ == '__main__':

    ligand_db_version = 'v1.7'

    charge_benchmark = ChargeBenchmark(true_charge_name='charge')
    charge_benchmark.calculate_scores_of_charge_benchmark(f'/Users/timosommer/PhD/projects/RCA/projects/DART/data/final_db_versions/complex_db_{ligand_db_version}.json')
"""
Implements a linear equation solver based on linear fitting for calculating charges.
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from sklearn.linear_model import LinearRegression, Ridge, SGDRegressor, Lasso
from copy import deepcopy
from sklearn.metrics import accuracy_score, r2_score
from datetime import datetime
from tqdm import tqdm
import matplotlib.pyplot as plt
# import statsmodels.api as sm
# from statsmodels.regression.quantile_regression import QuantReg
from sklearn.preprocessing import StandardScaler
from collections import Counter
from typing import Union
from itertools import product

from DARTassembler.src.ligand_extraction.io_custom import load_json

pd.options.mode.chained_assignment = None
from scipy.sparse import dok_array


class LinearChargeSolver:
    
    def __init__(self,
                 all_complexes: Union[str, dict],
                 save_dir: str,
                 n_test: Union[bool, int] = False,
                 ligand_data: Union[None, str, list] = None,
                 true_charge_name: str= '',
                 use_sparse_matrix: bool = True,
                 max_iterations: Union[None, int] = None,
                 lin_reg_alpha: float = 0.000000001,
                 lin_reg_max_iter: int = 50_000
                 ) -> None:
        """
        @param all_complexes: Either a dictionary or a path to a json of the complex database.
        @param save_dir: Where to save output
        @param n_test: Set to an integer to reduce the number of complexes for a test run.
        @param ligand_data: Ligand information to keep, mostly useful for debugging.
        @param true_charge_name: If ligand_data has a column for a ground truth charge, the name of this column to calculate scores.
        @param use_sparse_matrix: If True, uses a sparse matrix for X to drastically reduce runtime and memory.
        """
        self.save_dir = save_dir
        self.n_test = n_test
        self.y_name = 'charge_bias'
        self.true_charge_name = true_charge_name
        self.use_sparse_matrix = use_sparse_matrix
        self.max_iterations = max_iterations
        self.lin_reg_alpha = lin_reg_alpha
        self.lin_reg_max_iter = lin_reg_max_iter

        df_add_ligand_data = self.get_ligand_data(ligand_data)

        self.all_complexes = self.get_all_complexes(all_complexes, self.n_test)

        self.all_complexes_names = list(self.all_complexes.keys())
        self.all_complexes_indices = {c: idx for idx, c in enumerate(self.all_complexes_names)}

        self.unique_ligands = self.get_unique_ligands_props(add_ligand_data=df_add_ligand_data)
        self.has_true_charge = any([true_charge_name in lig for lig in self.unique_ligands.values()])
        self.unique_ligand_names = list(self.unique_ligands.keys())
        self.unique_ligand_indices = {lig: idx for idx, lig in enumerate(self.unique_ligand_names)}
        
        self.n_ligands = len(self.unique_ligand_names)
        self.n_complexes = len(self.all_complexes)

    def get_ligand_data(self, ligand_data: Union[None,str,list], id_col: str='graph_hash') -> pd.DataFrame:
        """
        Returns a df with data for each ligand which will be appended to the ligand df calculated in this linear solver. The input can either be a string to a csv or a list of string to a csv.
        :param ligand_data: None, path or list of paths to csv with ligand data.
        :return: pd.DataFrame with all ligand data merged.
        """
        df_ligand_data = pd.DataFrame()
        if ligand_data is None:
            return df_ligand_data

        if isinstance(ligand_data, str):
            ligand_data = [ligand_data]

        for i, csv_path in enumerate(ligand_data):
            df = pd.read_csv(csv_path)

            if not id_col in df:
                raise ValueError(f'Could not find column {id_col} in csv {csv_path} to identify ligands.')

            if df_ligand_data.empty:
                df_ligand_data = deepcopy(df)
            else:
                df_ligand_data = df_ligand_data.join(df, how='outer', on=id_col, rsuffix=f'_{i}')
                raise NotImplementedError('This if path has never been checked so please check if it works.')

        for col in ['graph_hash', 'unique_name']:
            assert not df_ligand_data.duplicated(subset=col).any(), f'Duplicates in column {col} in provided ligand df. Please provide a csv without duplicates.'
        return df_ligand_data

    def filter_input_complexes_inplace(self):
        init_n_complexes = len(self.all_complexes)
        self.removed_complexes_ids = []
        for c_id, c in self.all_complexes.items():
            
            has_no_ligands = c['ligands'] == []
            has_no_oxi_state = np.isnan(c['metal_oxi_state'])
            
            if has_no_ligands or has_no_oxi_state:
                self.removed_complexes_ids.append(c_id)
            
        for c_id in self.removed_complexes_ids:
            self.all_complexes.pop(c_id)
            
        end_n_complexes = len(self.all_complexes)
        print(f'Removed {init_n_complexes - end_n_complexes} from input complexes due to missing OS or ligands.')
        
        return
    
    def get_all_complexes(self, all_complexes: str, n_test) -> dict:

        try:
            self.all_complexes = load_json(all_complexes)
            print('Reading in complexes.')
        except TypeError:
            self.all_complexes = all_complexes

        self.filter_input_complexes_inplace()
        
        if n_test != False:
            test_complex_names = list(self.all_complexes.keys())[0:n_test]
            self.all_complexes = {c_id: c for c_id, c in self.all_complexes.items() if c_id in test_complex_names}
        
        return self.all_complexes
    
    def get_unique_ligands_props(self, add_ligand_data: pd.DataFrame) -> dict:
        all_ligand_keys = set()

        if add_ligand_data.empty:
            add_ligand_data = {}
        else:
            add_ligand_data = add_ligand_data.set_index(keys='unique_name').to_dict(orient='index')

        self.unique_ligands = {}
        for c_id, c in self.all_complexes.items():
            for lig in c['ligands']:
    
                name = lig["name"]
                uname = lig['unique_name']
                # graph_hash = lig['graph_hash']
                n_protons = lig['n_protons']
                
                stoichiometry = lig['stoichiometry']
                
                if not uname in self.unique_ligands:
                    self.unique_ligands[uname] = {
                                        # 'graph_hash': graph_hash,
                                        'stoichiometry': stoichiometry,
                                        'n_protons': n_protons
                                        }

                    if uname in add_ligand_data.keys():
                        new_info = add_ligand_data[uname]
                        assert all([new_info[key] == value for key, value in self.unique_ligands[uname].items()]), f'Values of dict-items with same key in both dictionaries do not match for unique ligand {uname}.'
                        self.unique_ligands[uname].update(new_info)
                        all_ligand_keys.update(self.unique_ligands[uname].keys())

        # Add NaN for each entry if not given otherwise to have a dictionary where all keys exists for all ligands.
        for lig, key in product(self.unique_ligands.values(), all_ligand_keys):
            if not key in lig:
                lig[key] = np.nan

        return self.unique_ligands
    
    def calculate_scores(self, max_abs_charge: int=10):
        if not self.has_true_charge:
            print(f'No true charge given in self.df_ligands. Cannot calculate accuracy and will skip this function.')
            return

        if not self.has_true_charge:
            print()

        delimiter = '============================'
        print(f'\nSCORES:\n{delimiter}')

        y_pred = deepcopy(self.df_ligands['pred_charge'])
        n_crazy = sum(y_pred.abs() > max_abs_charge)
        print(f'Num. pred. charges >{max_abs_charge}: {n_crazy}')

        try:
            frac_not_scale_invariant = sum(~self.df_ligands['charge_scale_invariant']) / len(self.df_ligands)
            print(f'Frac. not scale invariant: {frac_not_scale_invariant:.2g}')
        except KeyError:
            pass

        if self.has_true_charge:
            y_true = deepcopy(self.df_ligands[self.true_charge_name])
            valid = y_true.notna() & y_pred.notna()

            y_true_all, y_pred_all = y_true[valid], y_pred[valid]
            acc = accuracy_score(y_true=y_true_all, y_pred=y_pred_all)
            print(f'Accuracy (n={len(y_true_all)}): {acc:.2g}')

            confident = valid & self.df_ligands['is_confident']
            y_true_conf, y_pred_conf = y_true[confident], y_pred[confident]
            acc = accuracy_score(y_true=y_true_conf, y_pred=y_pred_conf)
            print(f'Confident accuracy (n={len(y_true_conf)}): {acc:.2g}')

        print(delimiter + '\n')

        return

    def get_tabular_data(self):
        other_data = {}

        if self.use_sparse_matrix:
            X = dok_array((self.n_complexes, self.n_ligands), dtype=int)
        else:
            X = np.zeros(shape=(self.n_complexes, self.n_ligands), dtype=int)
        
        for c_id in self.all_complexes_names:   # Don't use dict.items() to make sure order is correct.
            c = self.all_complexes[c_id]
            c_idx = self.all_complexes_indices[c_id]
            
            # Our versions of the tmQM, the tmQMg and the CSD have different names for the total charge.
            if 'total_charge' in c:
                self.charge_input_name = 'total_charge'
            elif 'total_q' in c:
                self.charge_input_name = 'total_q'
            elif 'charge' in c:
                self.charge_input_name = 'charge'
            else:
                raise ValueError('Could not identify total charge for each complex in input json.')
            
            other_data[c_id] = {
                'total_q': c[self.charge_input_name],
                'metal_oxi_state': c['metal_oxi_state'],
                self.y_name: c[self.charge_input_name] - c['metal_oxi_state']
                                }
            
            for lig in c['ligands']:
                uname = lig['unique_name']
                lig_idx = self.unique_ligand_indices[uname]
                X[c_idx, lig_idx] += 1
    
        df_other_data = pd.DataFrame.from_dict(data=other_data, orient='index')
        y = df_other_data[self.y_name].to_numpy()
        
        return X, y, df_other_data

    def get_complexes_with_converged_ligands(self, X_new, y_new, pred_charges_new, n_ligand_protons):
        # Exclude complexes if the sum of their ligands is not consistent with total charge and oxidation state
        total_ligand_charges = np.round(pred_charges_new)
        total_charges = (X_new*total_ligand_charges).sum(axis=1)
        consistent = y_new == total_charges

        # Exclude ligands with odd number of electrons
        odd_ligand_electrons = (-total_ligand_charges + n_ligand_protons)%2 == 1
        only_ligands_with_even_n_electrons = (X_new*odd_ligand_electrons).sum(axis=1) == 0

        converged = consistent & (only_ligands_with_even_n_electrons)

        return converged

    def iter_exclude_bad_ligands_and_refit(self, X, y, pred_charges):
        X_new, y_new = deepcopy(X), deepcopy(y)
        pred_charges_new = deepcopy(pred_charges)
        n_ligand_protons = np.array([self.unique_ligands[uname]['n_protons'] for uname in self.unique_ligand_names], dtype=int)
        for i in range(self.max_iterations):
            converged = self.get_complexes_with_converged_ligands(X_new=X_new, y_new=y_new, pred_charges_new=pred_charges_new, n_ligand_protons=n_ligand_protons)

            if not all(converged):
                print(f'\t- Not converged complexes: {100*sum(~converged)/len(converged):.2g}%')
                print(f'\t\t- Partial R2 score: {self.get_R2_score(X_new, y_new, pred_charges_new):.4g}')
                X_new, y_new = X_new[converged], y_new[converged]
                lin_reg = Lasso(fit_intercept=False, alpha=self.lin_reg_alpha, max_iter=self.lin_reg_max_iter)
                lin_reg.fit(X=X_new, y=y_new)

                # Do not update ligands with only excluded complexes since these will all have values of 0
                ligand_still_exists_after_converging = X_new.sum(axis=0) != 0
                pred_charges[ligand_still_exists_after_converging] = lin_reg.coef_[ligand_still_exists_after_converging]
            else:
                print(f'\t- Converged after {i} iteration(s).')
                print(f'\t\t- Partial R2 score: {self.get_R2_score(X_new, y_new, pred_charges_new):.4g}')
                break

        return pred_charges, ligand_still_exists_after_converging

    def get_predicted_charges_under_scaling(self, X, y, scale_charge_by: int=500):
        pseudo_lin_reg = Lasso(fit_intercept=False, alpha=self.lin_reg_alpha, max_iter=self.lin_reg_max_iter)
        pseudo_charges_per_ligand = np.sum(X, axis=1) * scale_charge_by
        y_pseudo = deepcopy(y) + pseudo_charges_per_ligand
        pseudo_lin_reg.fit(X=X, y=y_pseudo)
        scaled_predicted_charges = pseudo_lin_reg.coef_ - scale_charge_by

        return scaled_predicted_charges

    def get_fitted_charges(self, X, y):
        print('Do the fit.')
        start = datetime.now()

        # self.lin_reg = LinearRegression(fit_intercept=False)
        # # self.lin_reg = SGDRegressor(fit_intercept=False, loss='squared_error', max_iter=1_000, tol=1e-4, n_iter_no_change=100, epsilon=0, alpha=0)
        # self.lin_reg = Ridge(fit_intercept=False, alpha=0.0001)

        self.lin_reg = Lasso(fit_intercept=False, alpha=self.lin_reg_alpha, max_iter=self.lin_reg_max_iter)

        self.lin_reg.fit(X=X, y=y)
        pred_charges = self.lin_reg.coef_

        if not self.max_iterations is None:
            print('Repeat fit without failed complexes until convergence:')
            pred_charges, self.ligand_still_exists_after_converging = self.iter_exclude_bad_ligands_and_refit(X=X, y=y, pred_charges=pred_charges)
        else:
            self.ligand_still_exists_after_converging = np.array([True]*X.shape[1])

        try:
            print(f'Rank: {self.lin_reg.rank_}')
            self.singular_values = self.lin_reg.singular_
        except AttributeError:
            pass

        self.scaled_predicted_charges = self.get_predicted_charges_under_scaling(X=X, y=y)


        # self.lin_reg = sm.OLS(y, X, hasconst=False)
        # self.results = self.lin_reg.fit()
        # # print(self.results.summary())
        # pred_charges = self.results.params
        # self.standard_errors = self.results.bse
        # self.predicted_values = self.results.predict()
        # self.singular_values = self.lin_reg.wexog_singular_values

        
        end = datetime.now()
        print(f'Duration of ligand charge fit: {end - start}')
        return pred_charges

    def get_ligand_df(self):
        # print('Make ligand df.')
        self.df_ligands = {}
        
        for uname, pred_charge in zip(self.unique_ligand_names, self.pred_charges):
            lig = self.unique_ligands[uname]
            self.df_ligands[uname] = {
                                        'pred_charge_exact': pred_charge,
                                        'pred_charge': round(pred_charge)
                                        }

            # Occurences counts how often a ligand occurs in different complexes. Multiple occurrences in a single complex are counted only once.
            # If self.use_sparse_matrix is True, self.df will not have the full occupation matrix and the occurrences cannot be calculated.
            if not self.use_sparse_matrix:
                occurences = sum(self.df[uname] != 0)
                self.df_ligands[uname]['occurrences'] = occurences

            if self.has_true_charge:
                self.df_ligands[uname].update({
                                        'error': pred_charge - lig[self.true_charge_name],
                                        'abs_error': np.abs(pred_charge - lig[self.true_charge_name]),
                                        'int_error': round(pred_charge) - lig[self.true_charge_name],
                                        'abs_int_error': np.abs(round(pred_charge) - lig[self.true_charge_name]),
                                            })

            self.df_ligands[uname].update(lig)
            
        self.df_ligands = pd.DataFrame.from_dict(self.df_ligands, orient='index')
        
        if self.has_true_charge:
            self.df_ligands['correct'] = (self.df_ligands['pred_charge'] == self.df_ligands[self.true_charge_name])
            self.df_ligands.loc[self.df_ligands[self.true_charge_name].isna(), 'correct'] = None
        
        try:
            self.df_ligands['singular_value'] = self.singular_values
        except ValueError:
            print('More ligands than complexes, can not write singular values to ligands. This is no problem and improves when you use more complexes.')
        except AttributeError:
            pass

        try:
            self.df_ligands['scaled_pred_charge_exact'] = self.scaled_predicted_charges
            self.df_ligands['scaled_pred_charge'] = self.scaled_predicted_charges.round()
            self.df_ligands['charge_scale_invariant'] = self.df_ligands['scaled_pred_charge'] == self.df_ligands['pred_charge']
        except AttributeError:
            print('No scaled predicted charges found.')
        
        try:
            self.df_ligands['standard_error'] = self.standard_errors
            self.df_ligands['pvalues'] = self.results.pvalues
            self.df_ligands['tvalues'] = self.results.tvalues
            self.df_ligands['HC0_se'] = self.results.HC0_se
            self.df_ligands['HC1_se'] = self.results.HC1_se

        except AttributeError:
            pass
        
        self.df_ligands['dist_to_int'] = self.df_ligands['pred_charge_exact'] - self.df_ligands['pred_charge']
        self.df_ligands['abs_dist_to_int'] = self.df_ligands['dist_to_int'].abs()

        confidence, is_confident = self.get_confidence_of_prediction(self.df_ligands)
        self.df_ligands['confidence'] = confidence
        self.df_ligands['is_confident'] = is_confident
        
        return self.df_ligands
    
    def get_confidence_of_prediction(self, df_ligands, expected_charge_interval=(-4, 1), cutoff=0.6) -> tuple:
        # Scale confidence between 0 (very uncertain) and 1 (very confident).
        self.confidence = 1 - 2*deepcopy(df_ligands['abs_dist_to_int'])

        # Set uncertainty to 0 if the predicted charge is suspiciously high.
        charge_is_too_high = (df_ligands['pred_charge'] < expected_charge_interval[0]) | (df_ligands['pred_charge'] > expected_charge_interval[1])
        self.confidence.loc[charge_is_too_high] = 0

        self.confidence.loc[~self.ligand_still_exists_after_converging] = 0

        print(
f"""N confident ligands:
    - Convergence: {100 * sum(self.ligand_still_exists_after_converging) / len(self.confidence):.2g}%
    - High charge: {100 * sum(~charge_is_too_high) / len(self.confidence):.3g}%"""
            )
        try:
            charge_not_scale_invariant = ~self.df_ligands['charge_scale_invariant']
            self.confidence.loc[charge_not_scale_invariant] = 0
            print(f'\t- Scale invariance: {100 * sum(~charge_not_scale_invariant) / len(self.confidence):.2g}%')
        except KeyError:
            print('Could not find scaled predicted charges to compare with for uncertainty prediction.')

        self.is_confident = self.confidence >= cutoff
        print(f'  -> Total: {100 * sum(self.is_confident)/len(self.confidence):.2g}% ({sum(self.is_confident)})')

        return self.confidence, self.is_confident
    
    def check_standard_output(self, df_ligands):
        """Checks that the output dataframe adheres to the standard.
        """
        necessary_columns = ['unique_name', 'pred_charge', 'confidence']
        assert all([col in df_ligands.columns for col in necessary_columns]), f'The output of the ligand dataframe is not as described by the standard: It does\nt have one or more of the necessary columns ({necessary_columns}).'
        
        for col in ['unique_name']:
            assert not df_ligands[col].duplicated().any(), f'The output of the ligand dataframe is not as described by the standard: It has duplicates in the column {col} which should be unique.'
            
        return

    def get_R2_score(self, X, y, pred_charges):
        pred_charges = pred_charges.round()
        y_pred = (X*pred_charges).sum(axis=1)
        score = r2_score(y_true=y, y_pred=y_pred)

        return score

    def print_overall_scores(self):
        score = self.get_R2_score(self.X, self.y, self.pred_charges)
        print(f'R2 score all: {score:.4g}')

        try:
            conf_score = self.get_R2_score(self.X[:, self.is_confident], self.y, self.pred_charges[self.is_confident])
            print(f'R2 score confident: {conf_score:.4g}')
        except AttributeError:
            pass

        return

    def get_df_with_all_tabular_data(self, df_other_data: pd.DataFrame,  X: np.array) -> pd.DataFrame:
        if not X is None:
            df_X = pd.DataFrame(data=X, columns=self.unique_ligand_names, index=self.all_complexes_names)
            dfs = [df_other_data, df_X]
        else:
            dfs = [df_other_data]

        self.df = pd.concat(dfs, axis=1)
        
        try:
            self.df['predictions'] = self.predicted_values
        except AttributeError:
            pass
        
        return self.df
    
    def read_df_with_all_tabular_data_from_csv(self, filepath: str):
        self.df = pd.read_csv(filepath, index_col=0)
        
        return self.df
    
    def fit_charge_solver(self):
        
        self.X, self.y, df_other_data = self.get_tabular_data()
    
        self.pred_charges = self.get_fitted_charges(self.X, self.y)

        X = None if self.use_sparse_matrix else self.X
        self.get_df_with_all_tabular_data(df_other_data=df_other_data, X=X)
        
        return
    
    def plot_abs_error_vs_dist_to_int(self, plot_dir):
        if not self.has_true_charge:
            print(f'No true charge given in self.df_ligands. Skip plotting abs_error vs dist_to_int.')
            return

        small_error = (self.df_ligands['abs_error'].abs() < 10) & self.df_ligands['abs_error'].notna()
        df_ligands = self.df_ligands[small_error]
        
        df_ligands.plot(x='abs_dist_to_int', y='abs_error', kind='scatter')
        # plt.yscale('symlog')
        
        # Save plot.
        outpath = Path(plot_dir, 'error_vs_dist.png')
        outpath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, dpi=300)
        return

    def plot_error_vs_singular_value(self, plot_dir):
        if not self.has_true_charge:
            print(f'No true charge given in self.df_ligands. Skip plotting error vs singular_value.')
            return
        if not 'singular_value' in self.df_ligands:
            print(f'No singular value given in self.df_ligands. Skip plotting error vs singular_value.')
            return

        small_error = (self.df_ligands['error'].abs() < 10) & self.df_ligands['error'].notna()
        df_ligands = self.df_ligands[small_error]
        df_ligands['error'] = df_ligands['error'].abs()
    
        df_ligands.plot(x='singular_value', y='error', kind='scatter')
    
        # Save plot.
        outpath = Path(plot_dir, 'error_vs_singular_value.png')
        outpath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, dpi=300)
        return

    def plot_abs_int_error_vs_standard_error(self, plot_dir):
        if not 'standard_error' in self.df_ligands:
            print(f'No standard error given in self.df_ligands. Skip plotting abs_int_error vs standard_error')
            return
        if not self.has_true_charge:
            print(f'No true charge given in self.df_ligands. Skip plotting abs_int_error vs standard_error.')
            return

        small_error = (self.df_ligands['abs_int_error'].abs() < 10) & self.df_ligands['abs_int_error'].notna()
        df_ligands = self.df_ligands[small_error]
        df_ligands['abs_int_error'] = df_ligands['abs_int_error'].abs()
    
        df_ligands.plot(x='standard_error', y='abs_int_error', kind='scatter')
    
        # Save plot.
        outpath = Path(plot_dir, 'abs_int_error_vs_standard_error.png')
        outpath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, dpi=300)
        return
    
    def add_ligand_charge_to_all_complexes(self, df_ligands: pd.DataFrame):
        """
        Adds the ligand charge to each ligand in each complex in self.all_complexes.
        :return:
        """
        for c_id, c in self.all_complexes.items():
            for lig in c['ligands']:
                lig_name = lig['unique_name']
                lig_values = df_ligands.loc[lig_name]
                
                # Assign values to ligand, which will be also saved in self.all_complexes because the ligand here is exactly the same object as in self.all_complexes.
                lig['charge_LCS'] = float(lig_values.get('pred_charge'))
                lig['abs_dist_to_int'] = float(lig_values.get('abs_dist_to_int'))
                if 'standard_error' in lig_values:
                    lig['standard_error'] = float(lig_values.get('standard_error'))
        
        return
    
    def save_all_complexes_as_json(self):
        Path(self.save_dir).mkdir(parents=True, exist_ok=True)

        outpath = Path(self.save_dir, 'all_complexes_with_LCS_charges.json')
        print(f'Write all complexes to {outpath}')
        with open(outpath, 'w') as file:
            json.dump(self.all_complexes, file)
        print('Finished writing to json.')
        
        return
    
    def save_dfs_to_csv(self):
        Path(self.save_dir).mkdir(parents=True, exist_ok=True)
        
        outpath = Path(self.save_dir, 'df_ligands.csv')
        self.df_ligands.to_csv(outpath)
        
        outpath = Path(self.save_dir, 'df.csv')
        self.df.to_csv(outpath)
        
        return
    
    def print_statistics(self, corr_method='spearman'):
        corr = self.df_ligands.corr(method=corr_method, numeric_only=True)

        if 'standard_error' in self.df_ligands.columns and self.has_true_charge:
            print(f'Corr of abs_int_error and se:')
            print(f"{corr._get_value('abs_int_error', 'standard_error'):.2f}")
        else:
            print('Either standard error or true charge are not specified. Therefore skip printing statistics.')

        return
    
    def calculate_unique_ligand_charges(self, output_uncertain_charges_as_nan: bool=False) -> pd.DataFrame:
        """
        Runs the whole process for calculating ligand charges and returns a df with all of them.
        :return: pd.DataFrame with unique ligand names, hashes, charges and uncertainties.
        """
        self.fit_charge_solver()
        df_ligands = deepcopy(self.get_ligand_df())

        self.print_overall_scores()

        if output_uncertain_charges_as_nan:
            df_ligands.loc[~df_ligands['is_confident'], 'pred_charge'] = np.nan

        df_ligands = df_ligands.reset_index(names='unique_name')
        
        self.check_standard_output(df_ligands)
        return df_ligands
    

if __name__ == '__main__':
    
    n_test = False
    db_version = '1.4_all'
    all_complexes_path = f'../../data/final_db_versions/complex_db_v{db_version}.json'
    all_unique_ligands_path = [f'../../test/debug/databases/charge_benchmark/all_ligand_charges_with_high_confidence_v{db_version}.csv']
    save_dir = '../../../../data/linear_charge_solver/output/'
    plot_dir = Path(save_dir, 'plots')
    
    solver = LinearChargeSolver(all_complexes=all_complexes_path, save_dir=save_dir, n_test=n_test, ligand_data=all_unique_ligands_path, true_charge_name='charge')
    df_ligands = solver.calculate_unique_ligand_charges()
    solver.add_ligand_charge_to_all_complexes(df_ligands=solver.df_ligands)
    # solver.save_dfs_to_csv()
    # solver.save_all_complexes_as_json()
    
    acc = solver.calculate_scores()
    solver.plot_abs_error_vs_dist_to_int(plot_dir)
    solver.plot_abs_int_error_vs_standard_error(plot_dir)
    solver.plot_error_vs_singular_value(plot_dir)
    solver.print_statistics()
    
    
    
    # solver.make_plots(ligands=['unq_CSD-ANAQUC-01-a'])
    
    # Read existing df with full run (i.e. n_test = False)
    # full_run_df_path = Path('../../data/linear_charge_fitting/221111_output_testing=False_tmQM/df.csv')
    # full_run_df_ligands_path = Path('../../data/linear_charge_fitting/221111_output_testing=False_tmQM/df_ligands.csv')
    # if full_run_df_path.exists() and full_run_df_ligands_path.exists():
    #     full_df = pd.read_csv(full_run_df_path, index_col=0)
    #     full_df_ligands = pd.read_csv(full_run_df_ligands_path, index_col=0)
    
    # Check if output unchanged.
    if n_test == 2_000:
        print('Check output.')
        df_old = pd.read_csv(Path(save_dir, f'221129_df_n_test={n_test}.csv'), index_col=0)
        df_ligands_old = pd.read_csv(Path(save_dir, f'221129_df_ligands_n_test={n_test}.csv'), index_col=0)
        df_ligands_old = df_ligands_old.drop(columns=['true_charge', 'error', 'abs_error', 'int_error', 'abs_int_error', 'pred_charge_exact', 'singular_value', 'uncertainty'])


        pd.testing.assert_frame_equal(solver.df_ligands[df_ligands_old.columns], df_ligands_old, check_like=True)
        pd.testing.assert_frame_equal(solver.df[df_old.columns], df_old, check_like=True)
        print('All output is unchanged.')
    
    print('Done!')


    
    


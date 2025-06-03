"""
This script makes a benchmark to compare the MetaLig charges with the ones from the tmQMg-l and the Realigands.
"""
from collections import defaultdict
from pathlib import Path
import pandas as pd
import ast
import pysmiles
import logging
logging.getLogger('pysmiles').setLevel(logging.CRITICAL)  # Anything higher than warning
import ase

from DARTassembler.src.ligand_extraction.utilities_Molecule import stoichiometry2atomslist
from DARTassembler.src.ligand_extraction.utilities_graph import get_graph_hash

#%%
def get_graph_hash_from_smiles(smiles: str, smiles_metal_bond_node_idx_groups: list[list[int]] = None) -> str:
    """
    Create a graph with metal from the smiles and the smiles_metal_bond_node_idx_groups.
    """
    graph = pysmiles.read_smiles(smiles, explicit_hydrogen=True)
    if not smiles_metal_bond_node_idx_groups is None:
        # Add the metal to the graph
        metal_idx = len(graph.nodes)
        graph.add_node(metal_idx, element='Hg')

        for idx_group in smiles_metal_bond_node_idx_groups:
            for idx in idx_group:
                graph.add_edge(metal_idx, idx)

    graph_hash = get_graph_hash(graph, node_attr='element')

    return graph_hash

def match_ligands_with_metalig_by_complex_and_atoms(row, df, csd_complexes_to_ligands):
    all_csd_ids = [row.csd_ids] if isinstance(row.csd_ids, str) else row.csd_ids
    all_fitting_metalig_ligands = set()
    for csd_id in all_csd_ids:
        try:
            metalig_ligands_of_this_complex = csd_complexes_to_ligands[csd_id]
        except KeyError:
            continue
        metalig_ligands_of_this_complex = set(metalig_ligands_of_this_complex)

        # Get the atoms from the database
        db_atoms = row.atoms
        db_donors = row.donors
        try:    # If the atoms are a string, convert them to a list
            db_atoms = stoichiometry2atomslist(db_atoms)
        except TypeError:
            pass
        try:
            db_donors = stoichiometry2atomslist(db_donors)
        except TypeError:
            pass

        fitting_metalig_ligands = []
        for metalig_ligand in metalig_ligands_of_this_complex:
            metalig_atoms = df.loc[metalig_ligand, 'atoms']
            metalig_donors = df.loc[metalig_ligand, 'donors']
            same_atoms = sorted(db_atoms) == sorted(metalig_atoms)
            same_donors = sorted(db_donors) == sorted(metalig_donors)
            if same_atoms and same_donors:
                fitting_metalig_ligands.append(metalig_ligand)

        # Only one ligand in the CSD complex should fit, otherwise we need to ignore this ligand
        if len(fitting_metalig_ligands) == 1:
            all_fitting_metalig_ligands.add(fitting_metalig_ligands[0])

    all_fitting_metalig_ligands = list(all_fitting_metalig_ligands)
    if len(all_fitting_metalig_ligands) == 1:
        return all_fitting_metalig_ligands.pop()
    else:
        return None


#%%
if __name__ == '__main__':

    # Name the databases to compare. We also tested 'octlig', but it is not possible to match uniquely because the CSD IDs of the parent complexes are not given, therefore it is not included in the comparison and that code is outcommented.
    other_dbs = ['tmqmgl', 'cell2mol', 'realigands']
    print('Start script.')

    #%% Read in the MetaLig and preprocess it.
    df = pd.read_csv(Path('data', 'metalig.csv'))
    df = df.set_index('unique_name')
    df['donors'] = df['donors'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    df['atoms'] = df['atoms'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    df['all_csd_complex_ids'] = df['all_csd_complex_ids'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    df['donors_and_graph_hash'] = df['donors'].apply(lambda x: '-'.join(sorted(x))) + ':' + df['graph_hash']
    not_unique_donors_and_graph_hashes = set(df['donors_and_graph_hash'][df['donors_and_graph_hash'].duplicated()].unique())

    #%% Get dict mapping each CSD complex ID to the ligands in the MetaLig
    csd_complexes_to_ligands = defaultdict(list)
    for row in df.itertuples():
        for csd_id in row.all_csd_complex_ids:
            csd_complexes_to_ligands[csd_id].append(row.Index)
    csd_complexes_to_ligands = dict(csd_complexes_to_ligands)

    data = defaultdict(dict)
    for db in other_dbs:
        df_db = pd.read_csv(Path('data', db+'.csv'))
        data[db]['n_total'] = len(df_db)

        # Drop charges that are NaN
        n_nan = df_db['charge'].isna().sum()
        if n_nan > 0:
            print(f'{db}: Dropping {n_nan} ({100*n_nan/len(df_db):.2g}%) ligands with NaN charge.')
            df_db = df_db.dropna(subset=['charge'])
        df_db['charge'] = df_db['charge'].astype(int)

        # Preprocess the database df
        if db == 'tmqmgl':
            alternative_charges = df_db['is_alternative_charge'] != 0
            print(f'{db}: Dropping {alternative_charges.sum()} ({100 * alternative_charges.sum() / len(df_db):.2g}%) ligands with alternative charges.')
            df_db = df_db[~alternative_charges].drop(columns=['is_alternative_charge'])
            df_db['csd_ids'] = df_db['parent_metal_occurrences'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
            df_db['csd_ids'] = df_db['csd_ids'].apply(lambda x: [name.split('-')[0] for l in x.values() for name in l])
            df_db['atoms'] = df_db['stoichiometry']
            donor_cols = [col for col in df_db.columns if col.startswith('dentic_') or col.startswith('haptic_')]
            df_db['donors'] = df_db.apply(lambda row: ''.join([col.split('_')[-1]+f'{row[col]}' for col in donor_cols if row[col]>0]), axis=1)
        elif db == 'realigands':
            drop_ligands = ['U1*.xyz']
            df_db = df_db[~df_db['Name'].isin(drop_ligands)]
            df_db = df_db.rename(columns={'Name': 'name'})
            df_db['csd_ids'] = df_db['name'].apply(lambda x: x.split('_')[0].removeprefix('U'))
            df_db['donors'] = df_db['name'].apply(lambda x: x.split('_')[4].removesuffix('*.xyz'))
            df_db['atoms'] = df_db['atoms'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
        elif db == 'cell2mol':
            df_db['atoms'] = df_db['atoms'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
        elif db == 'octlig':
            raise NotImplementedError('The octlig database is not included in the comparison because it does not have the CSD IDs of the parent complexes.')
            # graph_hashs = {}
            # for row in tqdm(df_db.itertuples(), desc=f'Calculating graph_hash for {db}'):
            #     graph_hashs[row.SMILES] = get_graph_hash_from_smiles(row.SMILES)
            # df_db['graph_hash'] = df_db['SMILES'].map(graph_hashs)
            # df_db['donors'] = df_db['catoms'].apply(lambda x: '-'.join(sorted(x.split(','))))
            # df_db['donors_and_graph_hash'] = df_db['donors'] + ':' + df_db['graph_hash']
            # # Make the graph hash pd.NA if the graph hash in the df is not unique to not match these
            # df_db['not_unique_donors_and_graph_hash'] = df_db['donors_and_graph_hash'].isin(not_unique_donors_and_graph_hashes)
            # df_db['donors_and_graph_hash'] = df_db['donors_and_graph_hash'].where(~df_db['not_unique_donors_and_graph_hash'], other=pd.NA)
            # 'unique_name' = 'donors_and_graph_hash'
        else:
            raise ValueError(f'DB {db} not recognized.')

        # Drop ligands with either missing csd_ids, atoms or donors
        missing = df_db['csd_ids'].isna() | df_db['atoms'].isna() | df_db['donors'].isna()
        print(f'{db}: Dropping {missing.sum()} ({100*missing.sum()/len(df_db):.2g}%) ligands with missing `csd_ids`, `atoms` or `donors`.')
        df_db = df_db[~missing]
        # For each ligand in the other database, find the matching ligand in the MetaLig by CSD IDs and atoms.
        df_db['unique_name'] = df_db.apply(
            lambda row: match_ligands_with_metalig_by_complex_and_atoms(row, df, csd_complexes_to_ligands), axis=1)
        # Drop ligands where no MetaLig match was found
        print(f'{db}: Dropping {df_db[df_db['unique_name'].isna()].shape[0]} ({100*df_db[df_db['unique_name'].isna()].shape[0]/len(df_db):.2g}%) ligands with no match in the MetaLig.')
        df_db = df_db.dropna(subset='unique_name')
        # Drop ligands where more than one MetaLig match was found
        duplicated = df_db.duplicated(subset='unique_name', keep=False)
        print(f'{db}: Dropping {duplicated.sum()} ({100*duplicated.sum()/len(df_db):.2g}%) ligands where more than one match was found.')
        df_db = df_db.drop_duplicates(subset='unique_name', keep=False)

        # Merge the ligands with the MetaLig
        assert df_db['unique_name'].is_unique, f'{db}: The column `unique_name` is not unique in the database. This should not happen.'
        assert df_db['unique_name'].notna().all(), f'{db}: The column `unique_name` has NaN values. This should not happen.'
        df_db = df_db[['charge', 'unique_name']]    # Only merge the columns we need
        df_db = df_db.set_index('unique_name')
        df = df.merge(df_db, left_index=True, right_index=True, how='left', suffixes=(None, '_' + db))
        diff_col = 'charge_diff_' + db
        same_col = 'same_charge_' + db
        df[diff_col] = df['charge'] - df['charge_' + db]
        df[same_col] = df[diff_col] == 0
        df[same_col] = df[same_col].where(df[diff_col].notna(), other=pd.NA)

        # Calculate accuracy and other data
        n_agreed = df[same_col].sum()
        n_matched = df[diff_col].notna().sum()
        data[db]['agreement'] = n_agreed / n_matched
        data[db]['n_diff'] = n_matched - n_agreed
        data[db]['n_matched'] = n_matched
        data[db]['n_compared'] = len(df_db)

    #%% Change order of the columns
    first_cols = ['stoichiometry', 'geometry', 'donors', 'n_haptic_atoms', 'n_ligand_instances', 'charge', 'charge_tmqmgl', 'charge_realigands', 'charge_cell2mol', 'charge_diff_tmqmgl', 'charge_diff_realigands', 'charge_diff_cell2mol',  'all_csd_complex_ids']
    df = df[first_cols + [col for col in df.columns if col not in first_cols]]

    # Print the results by converting to a DataFrame
    df_results = pd.DataFrame.from_dict(data, orient='index')
    df_results = df_results.sort_values('n_total', ascending=False)
    print(df_results)

    # Save a csv with all the ligands that have different charges for all databases.
    df_all_diffs = []
    for db in other_dbs:
        df_diff = df[df['same_charge_' + db] == False].copy()
        df_diff = df_diff.sample(frac=1, random_state=0)    # Shuffle the order of ligands
        df_diff = df_diff.reset_index(drop=False)
        df_diff['other_db'] = db
        df_diff['donors'] = df_diff['donors'].apply(lambda x: '-'.join(sorted(x)))
        df_diff['all_csd_complex_ids'] = df_diff['all_csd_complex_ids'].apply(lambda x: ', '.join(sorted(x)))
        df_diff = df_diff.rename(columns={'charge_' + db: 'other_charge', 'charge': 'metalig_charge', 'unique_name': 'metalig_name'})
        df_diff['other_charge'] = df_diff['other_charge'].astype(int)
        cols = ['metalig_name', 'stoichiometry', 'geometry', 'donors', 'n_haptic_atoms', 'n_ligand_instances', 'metalig_charge', 'other_charge', 'other_db', 'all_csd_complex_ids']
        df_diff = df_diff[cols]
        df_all_diffs.append(df_diff)
    df_all_diffs = pd.concat(df_all_diffs, ignore_index=True)
    df_all_diffs.to_csv(Path('data', 'diff_charges.csv'), index=False)

    # # %% Compare the ligands with an old df
    # old_df_all_diffs = pd.read_csv(Path('data', 'OLD_diff_charges.csv'))
    # for db in other_dbs:
    #     old_df_diffs = old_df_all_diffs[old_df_all_diffs['other_db'] == db].copy()
    #     new_df_diffs = df_all_diffs[df_all_diffs['other_db'] == db].copy()
    #     merged_diff = old_df_diffs.merge(new_df_diffs, on=['metalig_name', 'stoichiometry', 'geometry', 'donors', 'n_haptic_atoms', 'n_ligand_instances', 'metalig_charge', 'other_charge'], how='inner', indicator=True)
    #     print(f'{db}: {len(merged_diff)}/{len(new_df_diffs)} ligand entries are identical.')

    print('Done!')




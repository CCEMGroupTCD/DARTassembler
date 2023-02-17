"""
This script makes statistical plots of the unique ligand database.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from src01.io_custom import load_unique_ligand_db
from pathlib import Path
sns.set_theme()

if __name__ == '__main__':

    db_version = '1.3'
    save_plots_dir = f'../data/db_statistics/unique_ligand_statistics/v{db_version}'
    db_path = f'../data/final_db_versions/unique_ligand_db_v{db_version}.json'



    save_plots_dir = Path(save_plots_dir)
    save_plots_dir.mkdir(parents=True, exist_ok=True)
    data = pd.DataFrame.from_dict(load_unique_ligand_db(path=db_path), orient='index')

    #%% plot histograms
    props = {
                'denticity': {'xlog': False, 'ylog': True},
                'occurrences': {'xlog': True, 'ylog': True, 'bins': 20, 'discrete': False},
                'n_denticities': {'xlog': False, 'ylog': True},
                'n_metals': {'xlog': False, 'ylog': True},
                'pred_charge': {'xlog': False, 'ylog': True}
            }
    for prop, configs in props.items():
        plt.figure()
        assert prop in data
        bins = 'auto' if not 'bins' in configs else configs['bins']
        discrete = 'True' if not 'discrete' in configs else configs['discrete']
        sns.histplot(data=data, x=prop, discrete=discrete, log_scale=configs['xlog'], bins=bins)
        if configs['ylog']:
                plt.yscale('log')
        plt.title(f'Distribution of {prop} in unique ligands')
        save_path = Path(save_plots_dir, f'hist_{prop}.png')
        plt.savefig(fname=save_path, dpi=300)
        plt.close()

    #%% Make histogram of all coordinated elements
    all_ligand_to_metal_elements = data['local_elements'].apply(lambda x: np.unique(x).tolist()).apply(pd.Series).stack()
    elements_hist = all_ligand_to_metal_elements.value_counts().rename('Count').to_frame().reset_index(names='coordinated_elements')
    plt.figure()
    sns.barplot(data=elements_hist, x='coordinated_elements', y='Count', color='b')
    plt.yscale('log')
    plt.ylim(1)
    plt.title('Elements coordinated to the metal')
    save_path = Path(save_plots_dir, f'hist_coordinated_elements.png')
    plt.savefig(fname=save_path, dpi=300)
    plt.close()

    #%% Make histogram of all elements of the unique ligands
    all_unique_elements = data['atomic_props'].apply(lambda x: np.unique(x['atoms']).tolist()).apply(pd.Series).stack()
    all_elements_hist = all_unique_elements.value_counts().rename('Count').to_frame().reset_index(names='elements')
    plt.figure()
    sns.barplot(data=all_elements_hist, x='elements', y='Count', color='b')
    plt.yscale('log')
    plt.ylim(1)
    plt.title('Elements contained in unique ligands')
    save_path = Path(save_plots_dir, f'hist_elements.png')
    plt.savefig(fname=save_path, dpi=300)
    plt.close()

    #%% Make histogram of partial charges and coordinates
    props = ['x',  'y', 'z', 'partial_charge']
    for prop in props:
        all_props = data['atomic_props'].apply(lambda x: x[prop]).explode().rename(prop).to_frame()
        plt.figure()
        sns.histplot(data=all_props, x=prop)
        plt.yscale('log')
        plt.title(f'Distribution of {prop} in unique ligands')
        save_path = Path(save_plots_dir, f'hist_{prop}.png')
        plt.savefig(fname=save_path, dpi=300)
        plt.close()

    #%% Make pie plot of where unique ligands are lost/ filtered
    uligs_numbers = {
                        'no oxidation state\nof complex': data['pred_charge'].isna().sum(),
                        'charge not confident': ((~data['pred_charge_is_confident']) & data['pred_charge'].notna()).sum(),
                        'confident charge assigned': data['pred_charge_is_confident'].sum()
    }
    assert sum(uligs_numbers.values()) == len(data)
    plt.figure()
    colors = sns.color_palette()
    numbers = list(uligs_numbers.values())
    labels = list(uligs_numbers.keys())
    plt.pie(x=numbers, labels=labels, explode=(0,0,0.03), autopct=lambda pct: f'{pct/100*sum(numbers):.0f} ({pct:.0f}%)')
    plt.title('LCS unique ligand losses')
    save_path = Path(save_plots_dir, f'pie_ligand_charges.png')
    plt.savefig(fname=save_path, dpi=300)
    plt.close()

    #%% Occurrences of unique ligands for which no oxidation state was given.
    plt.figure()
    sns.histplot(data=data[data['pred_charge'].isna()], x='occurrences', discrete=True)
    plt.yscale('log')
    plt.title(f'Occurrences of unique ligands of complexes without OS')
    save_path = Path(save_plots_dir, f'hist_occurrences_without_os.png')
    plt.savefig(fname=save_path, dpi=300)
    plt.close()
    print('Done!')

    #%% Ratio of hydrogen to carbon
    plt.figure()
    data['#H'] = data['atomic_props'].apply(lambda x: sum([el == 'H' for el in x['atoms']]))
    data['#C'] = data['atomic_props'].apply(lambda x: sum([el == 'C' for el in x['atoms']]))
    data['#H/#C'] = data['#H'] / data['#C']
    sns.scatterplot(data=data, x='#C', y='#H/#C')
    plt.title(f'Ratio of H to C')
    save_path = Path(save_plots_dir, f'hist2D_HC_ratio.png')
    plt.savefig(fname=save_path, dpi=300)
    plt.close()
    print('Done!')


"""
This script makes statistical plots of the unique ligand database.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
from DARTassembler.src.ligand_extraction.utilities import unroll_dict_into_columns
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db
from pathlib import Path
sns.set_theme()
plt.rcParams['svg.fonttype'] = 'none' # for correct text rendering in some programs
import pysmiles


def flatten(l):
    return [item for sublist in l for item in sublist]

if __name__ == '__main__':

    db_version = '1.7'
    save_plots_dir = f'../../data/db_statistics/unique_ligand_statistics/v{db_version}'
    db_path = f'../../data/final_db_versions/unique_ligand_db_v{db_version}.json'
    exclude_unconnected_ligands = True
    exclude_uncertain_charges = True



    save_plots_dir = Path(save_plots_dir)
    save_plots_dir.mkdir(parents=True, exist_ok=True)
    data_dict = load_unique_ligand_db(path=db_path)
    data_raw = pd.DataFrame.from_dict(data_dict, orient='index')


    data = data_raw.copy()
    if exclude_unconnected_ligands:
        data = data[data['denticity'] > 0]
    if exclude_uncertain_charges:
        data = data[data['pred_charge_is_confident']]
    data = unroll_dict_into_columns(data, dict_col='global_props', prefix='gbl_', delete_dict=True)
    data = unroll_dict_into_columns(data, dict_col='stats', prefix='stats_', delete_dict=True)
    data['n_electrons'] = data['n_protons'] - data['pred_charge']
    data['odd_n_electrons'] = data['n_electrons'].apply(lambda n: n%2 == 1)

    data_same_graph = data.groupby('graph_hash', sort=False).agg(list)
    data_same_graph['n_same_graph_charges'] = data_same_graph.apply(lambda row: len(pd.unique([charge for conf, charge in zip(row['pred_charge_is_confident'], row['pred_charge']) if conf])), axis=1)



    #%% plot histograms
    props = {
                'denticity': {'xlog': False, 'ylog': False},
                'occurrences': {'xlog': True, 'ylog': True, 'bins': 20, 'discrete': False},
                'n_same_graph_denticities': {'xlog': False, 'ylog': True},
                'n_same_graphs': {'xlog': False, 'ylog': True},
                'stats_min_distance_to_metal': {'xlog': True, 'ylog': True, 'discrete': False},
                'n_metals': {'xlog': False, 'ylog': True},
                'pred_charge': {'xlog': False, 'ylog': False},
                'gbl_n_atoms': {'xlog': False, 'ylog': False},
                'gbl_molecular_weight': {'xlog': False, 'ylog': False},
                'gbl_LCS_pred_charge_confidence': {'xlog': False, 'ylog': False, 'discrete': False},
            }
    for prop, configs in props.items():
        plt.figure()
        assert prop in data, f'Property {prop} is not a ligand property.'
        bins = 'auto' if not 'bins' in configs else configs['bins']
        discrete = 'True' if not 'discrete' in configs else configs['discrete']
        sns.histplot(data=data, x=prop, discrete=discrete, log_scale=configs['xlog'], bins=bins)
        if configs['ylog']:
                plt.yscale('log')
        plt.title(f'Distribution of {prop} in unique ligands')
        save_path = Path(save_plots_dir, f'hist_{prop}.svg')
        plt.savefig(fname=save_path)
        plt.close()

    #%% Make histogram of denticities comparing with/without confident charge
    plt.figure()
    prop = 'denticity'
    ax = sns.histplot(data=data, x=prop, discrete=True, label='all', legend=False)
    sns.histplot(data=data[data['pred_charge_is_confident']], x=prop, discrete=True, label='confident charge', legend=False, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::2], labels[::2])
    plt.title(f'Distribution of {prop} in unique ligands with/without confident charges')
    save_path = Path(save_plots_dir, f'hist_{prop}_charge_comp.svg')
    plt.savefig(fname=save_path)
    plt.close()

    #%% Make histogram of all coordinated elements
    all_ligand_to_metal_elements = data['local_elements'].apply(lambda x: np.unique(x).tolist()).apply(pd.Series).stack()
    elements_hist = all_ligand_to_metal_elements.value_counts().rename('Count').to_frame().reset_index(names='coordinated_elements')
    plt.figure()
    sns.barplot(data=elements_hist, x='coordinated_elements', y='Count', color='b')
    plt.yscale('log')
    plt.ylim(1)
    plt.title('Elements coordinated to the metal')
    save_path = Path(save_plots_dir, f'hist_coordinated_elements.svg')
    plt.savefig(fname=save_path)
    plt.close()

    #%% Make histogram of all elements of the unique ligands
    all_unique_elements = data['atomic_props'].apply(lambda x: np.unique(x['atoms']).tolist()).apply(pd.Series).stack()
    all_elements_hist = all_unique_elements.value_counts().rename('Count').to_frame().reset_index(names='elements')
    plt.figure()
    sns.barplot(data=all_elements_hist, x='elements', y='Count', color='b')
    plt.yscale('log')
    plt.ylim(1)
    plt.title('Elements contained in unique ligands')
    save_path = Path(save_plots_dir, f'hist_elements.svg')
    plt.savefig(fname=save_path)
    plt.close()

    #%% Make histogram of partial charges and coordinates
    props = ['x',  'y', 'z', 'partial_charge']
    for prop in props:
        try:
            all_props = data['atomic_props'].apply(lambda x: x[prop]).explode().rename(prop).to_frame()
            plt.figure()
            sns.histplot(data=all_props, x=prop)
            plt.yscale('log')
            plt.title(f'Distribution of {prop} in unique ligands')
            save_path = Path(save_plots_dir, f'hist_{prop}.svg')
            plt.savefig(fname=save_path)
            plt.close()
        except KeyError:
            print(f'Property {prop} not in atomic_props, it is skipped.')

    #%% Make pie plot of where unique ligands are lost/ filtered
    pie_data = data_raw[data_raw['denticity'] > 0]
    uligs_numbers = {
                        'no oxidation state\nof complex': pie_data['pred_charge'].isna().sum(),
                        'charge not confident': ((~pie_data['pred_charge_is_confident']) & pie_data['pred_charge'].notna()).sum(),
                        'confident charge assigned': pie_data['pred_charge_is_confident'].sum()
    }
    assert sum(uligs_numbers.values()) == len(pie_data)
    plt.figure()
    colors = sns.color_palette()
    numbers = list(uligs_numbers.values())
    labels = list(uligs_numbers.keys())
    plt.pie(x=numbers, labels=labels, explode=(0,0,0.03), autopct=lambda pct: f'{pct/100*sum(numbers):.0f} ({pct:.0f}%)')
    plt.title('LCS unique ligand losses')
    save_path = Path(save_plots_dir, f'pie_ligand_charges.svg')
    plt.savefig(fname=save_path)
    plt.close()

    #%% Occurrences of unique ligands for which no oxidation state was given.
    plt.figure()
    sns.histplot(data=data[data['pred_charge'].isna()], x='occurrences', discrete=True)
    plt.yscale('log')
    plt.title(f'Occurrences of unique ligands of complexes without OS')
    save_path = Path(save_plots_dir, f'hist_occurrences_without_os.svg')
    plt.savefig(fname=save_path)
    plt.close()

    #%% Ratio of hydrogen to carbon
    plt.figure()
    data['#H'] = data['atomic_props'].apply(lambda x: sum([el == 'H' for el in x['atoms']]))
    data['#C'] = data['atomic_props'].apply(lambda x: sum([el == 'C' for el in x['atoms']]))
    data['#H/#C'] = data['#H'] / data['#C']
    sns.scatterplot(data=data, x='#C', y='#H/#C')
    plt.title(f'Ratio of H to C')
    save_path = Path(save_plots_dir, f'hist2D_HC_ratio.svg')
    plt.savefig(fname=save_path)
    plt.close()

    #%% Histogram of denticities with donor types
    _, ax = plt.subplots()
    top_n = 8
    denticities = [1, 2, 3, 4, 5]
    data['donors'] = data['local_elements'].apply(lambda el: '-'.join(sorted(el)))
    donors = data[data['denticity'].isin(denticities)]
    top_donors = donors.groupby('denticity')['donors'].agg(lambda x: pd.value_counts(x).index.tolist()[:top_n-1])
    donors['donors'] = donors.apply(lambda row: row['donors'] if row['donors'] in top_donors.loc[row['denticity']] else 'others', axis=1)
    donor_dict = donors.groupby('denticity')['donors'].agg(lambda x: pd.value_counts(x).to_dict()).to_dict()
    for dent, dons in donor_dict.items():
        order = [don for don in dons if not don == 'others'] + ['others']
        donor_dict[dent] = {don: dons[don] for don in order}
    donor_data = []
    for dent, dons in donor_dict.items():
        for name, count in dons.items():
            donor_data.append({'denticity': dent, 'donors': name, 'count': count})
    donor_data = pd.DataFrame(donor_data)
    bottom = np.zeros(len(denticities))
    for idx in range(top_n):
        counts = np.array([list(donor_dict[dent].values())[idx] for dent in denticities])
        ax.bar(denticities, counts, bottom=bottom)
        bottom += counts
    for idx, c in enumerate(ax.containers):
        labels = [list(donor_dict[dent].keys())[idx] for dent in denticities]
        ax.bar_label(c, labels=labels, label_type='center')
    plt.title(f'Denticities & donor types')
    plt.xticks(denticities)
    # plt.legend()
    # plt.xscale('log')
    save_path = Path(save_plots_dir, f'hist_donors.svg')
    plt.savefig(fname=save_path)
    plt.close()

    print('Done!')

    #%% Make histogram of original metals
    plt.figure(figsize=(11.6, 5))
    sns.set(font_scale=0.8)
    metals = data['count_metals'].apply(lambda x: list(x.keys())).explode().value_counts().rename('Number of ligands').to_frame().reset_index(names='Originating metal')
    sns.barplot(data=metals, x='Originating metal', y='Number of ligands', color='b')
    plt.title(f'Original metals in unique ligands')
    # plt.yscale('log')
    # plt.xticks(fontsize=9)
    save_path = Path(save_plots_dir, f'hist_original_metals.svg')
    plt.savefig(fname=save_path)
    plt.close()
    sns.set(font_scale=1)



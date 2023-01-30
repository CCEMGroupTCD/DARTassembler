"""
This script makes statistical plots of the unique ligand database.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from src01.io_custom import load_complex_db
from pathlib import Path
from copy import deepcopy
sns.set_theme()

def unroll_dict_in_df_column(df:pd.DataFrame, column: str, prefix: str):
    df = deepcopy(df)
    df_unrolled = df[column].apply(pd.Series)
    df_unrolled = df_unrolled.rename(columns={col: f'{prefix}{col}' for col in df_unrolled})
    df = df.merge(df_unrolled, left_index=True, right_index=True, validate='1:1', suffixes=('', ''))

    return df

if __name__ == '__main__':

    db_version = '1.2'
    save_plots_dir = f'../data/db_statistics/complex_statistics/v{db_version}'
    db_path = f'../data/tmQMG_Jsons/complex_db_v{db_version}.json'



    save_plots_dir = Path(save_plots_dir)
    save_plots_dir.mkdir(parents=True, exist_ok=True)
    data = pd.DataFrame.from_dict(load_complex_db(path=db_path), orient='index')

    data = unroll_dict_in_df_column(df=data, column='global_props', prefix='gbl_')
    c = data.iloc[0].to_dict()

    #%% plot histograms of global props
    props = {
        'metal_oxi_state': {'discrete': True},
         'total_q': {'discrete': True},
         'Metal_q': {},
         'gbl_charge': {'discrete': True},
        'gbl_molecular_mass': {},
         'gbl_n_atoms': {'discrete': True},
         'gbl_n_electrons': {'discrete': True},
        'gbl_tzvp_lumo_energy': {},
        'gbl_tzvp_homo_energy': {},
        'gbl_tzvp_homo_lumo_gap': {},
        'gbl_homo_lumo_gap_delta': {},
        'gbl_tzvp_electronic_energy': {},
        'gbl_electronic_energy_delta': {},
        'gbl_tzvp_dispersion_energy': {},
        'gbl_dispersion_energy_delta': {},
        'gbl_enthalpy_energy': {},
        'gbl_enthalpy_energy_correction': {},
        'gbl_gibbs_energy': {},
         'gbl_gibbs_energy_correction': {},
         'gbl_zpe_correction': {},
        'gbl_heat_capacity': {},
         'gbl_entropy': {},
         'gbl_tzvp_dipole_moment': {},
        'gbl_dipole_moment_delta': {},
        'gbl_polarisability': {},
        'gbl_lowest_vibrational_frequency': {},
         'gbl_highest_vibrational_frequency': {'ylog': True},
            }
    for prop, configs in props.items():
        plt.figure()
        assert prop in data
        bins = 'auto' if not 'bins' in configs else configs['bins']
        discrete = False if not 'discrete' in configs else configs['discrete']
        log_scale = None if not 'xlog' in configs else configs['xlog']
        sns.histplot(data=data, x=prop, discrete=discrete, log_scale=log_scale, bins=bins)
        if 'ylog' in configs and configs['ylog']:
                plt.yscale('log')
        if 'xlabel_fontsize' in configs:
            plt.xticks(size=configs['xlabel_fontsize'])
        plt.title(f'Distribution of {prop} in complexes')
        save_path = Path(save_plots_dir, f'hist_{prop}.png')
        plt.savefig(fname=save_path, dpi=300)
        plt.close()

    #%% plot histogram of metals
    plt.figure()
    all_metals_hist = data['metal'].value_counts().rename('Count').to_frame().reset_index(names='metal')
    sns.barplot(data=all_metals_hist, x='metal', y='Count', color='b')
    plt.xticks(size=8)
    plt.title(f'Distribution of metal in complexes')
    save_path = Path(save_plots_dir, f'hist_metal.png')
    plt.savefig(fname=save_path, dpi=300)
    plt.close()


    print('Done!')




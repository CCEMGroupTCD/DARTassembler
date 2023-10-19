import os
import shutil
import warnings
import pymatgen.core as pt
import re
import pickle
import numpy as np
from typing import Union, Tuple
import json
from tqdm import tqdm
import random
from sklearn.metrics import r2_score, mean_absolute_error
import pandas as pd
import ase
import networkx as nx
import matplotlib
from pathlib import Path
from DARTassembler.src.constants.Periodic_Table import DART_Element
import seaborn as sns
from DARTassembler.src.ligand_extraction.io_custom import load_json
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import cclib
sns.set_style("ticks")



def hartree2ev(energy: float) -> float:
    return energy * 27.2114

if __name__ == "__main__":
    dft_csv = '../analyze_complexes/data.csv'
    xtb_csv = '231013_relaxations/xtb_relaxations.csv'
    xtb_structure_csv = '231013_relaxations_and_original_dirs/xtb_data.csv'
    outdir = 'plots'
    expected_donors = ['P', 'N', 'Br', 'C']
    homo = 'HOMO'
    lumo = 'LUMO'
    hlgap = 'HL_gap'

    # %% Plot histograms

    df_dft = pd.read_csv(dft_csv)
    df_dft['complex'] = df_dft['dir'].apply(lambda x: Path(x).name)
    df_dft = df_dft.rename(columns={'Metal charge': 'metal_charge'})

    df_xtb = pd.read_csv(xtb_csv)
    df_xtb[homo] = df_xtb[homo].apply(hartree2ev)
    df_xtb[homo] = df_xtb[homo].apply(hartree2ev)
    df_xtb[hlgap] = df_xtb[lumo] - df_xtb[homo]

    df_xtb_struct = pd.read_csv(xtb_structure_csv)
    df_xtb = df_xtb.merge(df_xtb_struct, on='complex')

    df = df_dft.merge(df_xtb, on='complex', suffixes=('_dft', '_xtb'))

    # Remove complexes with NaN values
    df = df.dropna(axis='rows')

    others = ['P-N bite angle', homo, lumo, hlgap, 'Metal charge']
    cols = expected_donors + others

    for prop in ['dist_P', 'dist_N', 'dist_Br', 'dist_C', 'P-N bite angle', homo, lumo, hlgap, 'metal_charge']:
        plt.figure()
        x, y = prop + '_dft', prop + '_xtb'
        sns.scatterplot(data=df, x=x, y=y)

        # r2 = r2_score(df[x], df[y])
        # mae = mean_absolute_error(df[x], df[y])
        # text = f'$r2$: {r2:.2f}\nMAE: {mae:.2f}'
        # plt.annotate(text, xy=(0.85, 0.9), xycoords='axes fraction')

        outpath = Path(outdir, f'comp_{prop}.png')
        plt.savefig(outpath)


    #
    # # Histograms
    # for donor in cols:
    #     plt.figure()
    #     name = f'dist_{donor}' if donor not in others else donor
    #     sns.histplot(data=df, x=name, hue='metal')
    #     plt.xlabel(donor)
    #     plt.savefig(Path(outdir, f'hist_{donor}.svg'))
    #     plt.close()
    #
    # # 2D scatter plots
    # for donor1, donor2 in [('P', 'N'), ('P', 'Br'), ('P', 'C'), ('N', 'Br'), ('N', 'C'), ('Br', 'C')]:
    #     for dist in ['dist', 'dist_norm']:
    #         plt.figure()
    #         sns.scatterplot(data=df, x=f'{dist}_{donor1}', y=f'{dist}_{donor2}', hue='metal', alpha=0.4)
    #         label = 'Distance' if dist == 'dist' else 'Normalized distance'
    #         plt.xlabel(f'{label} to {donor1} in Angstroms')
    #         plt.ylabel(f'{label} to {donor2} in Angstroms')
    #         plt.savefig(Path(outdir, f'scatter_{dist}_{donor1}_{donor2}.svg'))
    #         plt.close()
    #
    # # 2D scatter plots with Mulliken charge as hue
    # for donor1, donor2 in [('P', 'N'), ('P', 'Br'), ('P', 'C'), ('N', 'Br'), ('N', 'C'), ('Br', 'C')]:
    #     for dist in ['dist', 'dist_norm']:
    #         plt.figure()
    #         palette = sns.color_palette("coolwarm", as_cmap=True)#'coolwarm'
    #         markers = ['o', '^']
    #         ax = sns.scatterplot(data=df, x=f'{dist}_{donor1}', y=f'{dist}_{donor2}', hue='Metal charge', style='metal', palette=palette, alpha=1, markers=markers)
    #
    #         norm = plt.Normalize(df['Metal charge'].min(), df['Metal charge'].max())#-0.5, -0.2)#
    #         sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
    #         sm.set_array([])
    #
    #         # Remove the legend and add a colorbar
    #         # ax.get_legend().remove()
    #         ax.figure.colorbar(sm)
    #
    #         h, l = ax.get_legend_handles_labels()
    #         plt.legend(h[6:], l[6:])
    #
    #         label = 'Distance' if dist == 'dist' else 'Normalized distance'
    #         plt.xlabel(f'{label} to {donor1} in Angstroms')
    #         plt.ylabel(f'{label} to {donor2} in Angstroms')
    #         plt.savefig(Path(outdir, f'scatter_metal_charge_{dist}_{donor1}_{donor2}.svg'))
    #         plt.close()

    # Data analytics
    # df_stats = df.groupby('metal').agg(['max', 'min', 'mean', 'std']).T
    # print(df_stats)
    # df_corr = df.corr(method='spearman')
    print('Done!')



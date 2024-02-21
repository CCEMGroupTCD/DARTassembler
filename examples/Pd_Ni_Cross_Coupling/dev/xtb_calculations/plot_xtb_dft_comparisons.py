import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("ticks")



def hartree2ev(energy: float) -> float:
    return energy * 27.2114



if __name__ == "__main__":
    dft_csv = '../analyze_complexes/231019_data_all_extracted.csv'
    xtb_csv = 'data/231013_relaxations/xtb_relaxations.csv'
    xtb_structure_csv = 'data/231013_relaxations_and_original_dirs/xtb_data.csv'
    outdir = 'data/plots'
    expected_donors = ['P', 'N', 'Br', 'C']
    homo = 'homo'
    lumo = 'lumo'
    hlgap = 'hlgap'
    pn_bite_angle = 'P_N_bite_angle'

    # %% Plot histograms

    df_dft = pd.read_csv(dft_csv)
    df_dft['complex'] = df_dft['dir'].apply(lambda x: Path(x).name)
    df_dft = df_dft.rename(columns={'Metal charge': 'metal_charge', 'P-N bite angle': pn_bite_angle, 'HOMO': 'homo', 'LUMO': lumo, 'HL_gap': 'hlgap'})

    df_xtb = pd.read_csv(xtb_csv)
    df_xtb_struct = pd.read_csv(xtb_structure_csv)
    df_xtb = df_xtb.merge(df_xtb_struct, on='complex')
    df_xtb = df_xtb.rename(columns={'P-N bite angle': pn_bite_angle, 'HOMO': 'homo', 'LUMO': lumo, 'HL_gap': 'hlgap'})
    df_xtb[homo] = df_xtb[homo].apply(hartree2ev)
    df_xtb[lumo] = df_xtb[lumo].apply(hartree2ev)
    df_xtb[hlgap] = df_xtb[lumo] - df_xtb[homo]


    df = df_dft.merge(df_xtb, on='complex', suffixes=('_dft', '_xtb'))

    # Remove complexes with NaN values
    df = df.dropna(axis='rows')

    # Make output dir
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    others = [pn_bite_angle, homo, lumo, hlgap, 'Metal charge']
    cols = expected_donors + others
    alpha = 0.6
    df_corr = df.corr(method='pearson')

    for prop in ['dist_P', 'dist_N', 'dist_Br', 'dist_C', pn_bite_angle, homo, lumo, hlgap, 'metal_charge']:
        fig, ax = plt.subplots()
        x, y = prop + '_dft', prop + '_xtb'

        # diagonal line, except for homo-lumo gap because of large offset
        if not prop == hlgap:
            min_value = min(df[x].min(), df[y].min())
            max_value = max(df[x].max(), df[y].max())
            ax.plot([min_value, max_value], [min_value, max_value], color='black', linestyle='dashed', alpha=0.3)

        sns.scatterplot(data=df, x=x, y=y, alpha=alpha, hue='metal', ax=ax)

        corr = df_corr.loc[x, y]
        text = f'Pearson correlation: {corr:.2g}'
        plt.annotate(text, xy=(0.018, 0.81), xycoords='axes fraction')
        plt.legend(loc='upper left')



        outpath = Path(outdir, f'comp_{prop}.pdf')
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



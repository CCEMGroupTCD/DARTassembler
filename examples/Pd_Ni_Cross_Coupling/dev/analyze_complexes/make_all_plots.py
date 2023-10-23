import numpy as np
import pandas as pd
from pathlib import Path
from DARTassembler.src.constants.Periodic_Table import DART_Element
import seaborn as sns
import matplotlib.pyplot as plt


if __name__ == "__main__":
    input_csv = '231019_data_all_extracted.csv'
    outdir = 'plots'
    expected_donors = ['P', 'N', 'Br', 'C']

    # Define names
    pn_bite_angle = 'P_N_bite_angle'
    homo = 'homo'
    lumo = 'lumo'
    hlgap = 'hlgap'
    metal_charge = 'metal_charge'

    # %% Plot histograms
    sns.set_style("ticks")
    df = pd.read_csv(input_csv)

    # Assert no NaN
    assert df.isnull().values.any() == False, 'There are NaN values in the dataframe!'

    for donor in expected_donors:
        df[f'dist_norm_{donor}'] = df[f'dist_{donor}'] - df['metal'].apply(lambda el: DART_Element(el).covalent_radius_angstrom) - DART_Element(donor).covalent_radius_angstrom

    # Population plots
    norm_distances = [f'norm_{donor}' for donor in expected_donors]
    others = [pn_bite_angle, homo, lumo, hlgap, metal_charge]
    cols = expected_donors + norm_distances + others
    for donor in cols:
        name = f'dist_{donor}' if donor not in others else donor
        fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)

        bins = np.linspace(start=df[name].min(), stop=df[name].max(), num=20)
        for ax, metal, invert, color in zip(axes.ravel(), ['Pd', 'Ni'], [False, True], ['C0', 'C1']):
            metal_df = df[df['metal'] == metal]
            sns.histplot(data=metal_df, x=name, bins=bins, color=color, ax=ax, alpha=0.7, legend=True, label=metal)

        max_ytick = max([tick for ax in axes.ravel() for tick in ax.get_yticks()])
        min_ytick = min([tick for ax in axes.ravel() for tick in ax.get_yticks()])
        for ax in axes.ravel():
            ax.set_ylim(min_ytick, max_ytick)

        axes[1].invert_yaxis()
        plt.subplots_adjust(hspace=0)
        plt.legend(loc='lower left')

        plt.savefig(Path(outdir, f'pop_{donor}.svg'))
        plt.close()

    # Histograms
    for donor in cols:
        plt.figure()
        name = f'dist_{donor}' if donor not in others else donor
        sns.histplot(data=df, x=name, hue='metal', bins=25)
        plt.xlabel(donor)
        plt.savefig(Path(outdir, f'hist_{donor}.svg'))
        plt.close()

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

    # 2D scatter plots with Mulliken charge as hue
    for donor1, donor2 in [('P', 'N'), ('P', 'Br'), ('P', 'C'), ('N', 'Br'), ('N', 'C'), ('Br', 'C')]:
        for dist in ['dist']:
            plt.figure()
            palette = sns.color_palette("coolwarm", as_cmap=True)#'coolwarm'
            markers = ['o', '^']
            ax = sns.scatterplot(data=df, x=f'{dist}_{donor1}', y=f'{dist}_{donor2}', hue=metal_charge, style='metal', palette=palette, alpha=1, markers=markers)

            norm = plt.Normalize(df[metal_charge].min(), df[metal_charge].max())#-0.5, -0.2)#
            sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
            sm.set_array([])

            # Remove the legend and add a colorbar
            # ax.get_legend().remove()
            ax.figure.colorbar(sm)

            h, l = ax.get_legend_handles_labels()
            plt.legend(h[6:], l[6:])

            label = 'Distance' if dist == 'dist' else 'Normalized distance'
            plt.xlabel(f'{label} to {donor1} in Angstroms')
            plt.ylabel(f'{label} to {donor2} in Angstroms')
            plt.savefig(Path(outdir, f'scatter_metal_charge_{dist}_{donor1}_{donor2}.svg'))
            plt.close()

    # Data analytics
    df_stats = df.groupby('metal').agg(['max', 'min', 'mean', 'std']).T
    print(df_stats)
    df_corr = df.corr(method='spearman')
    print('Done!')



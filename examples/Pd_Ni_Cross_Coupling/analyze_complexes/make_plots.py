import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
sns.set_style("ticks")
# fix for correct text rendering in some programs
plt.rcParams['svg.fonttype'] = 'none'

# Set input data and output directory paths
input_data = 'data_DFT_relaxed_complexes.csv'
outdir = 'plots'


#======================  Plotting ======================#
if __name__ == '__main__':
    # Define names
    pn_bite_angle = 'P_N_bite_angle'
    homo = 'homo'
    lumo = 'lumo'
    hlgap = 'hlgap'
    metal_charge = 'metal_charge'
    metal = 'metal'

    df = pd.read_csv(input_data)
    # Remove not converged complexes with NaN values
    df = df.dropna(axis='rows')

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Histograms
    for name in [pn_bite_angle, homo, lumo, hlgap, metal_charge]:
        plt.figure()
        sns.histplot(data=df, x=name, hue=metal, bins=25)
        plt.xlabel(name)
        plt.savefig(Path(outdir, f'hist_{name}.svg'))
        plt.close()

    # 2D scatter plot of Metal-P and Metal-N bond distance with Mulliken charge as hue
    donor1 = 'P'
    donor2 = 'N'
    plt.figure()
    palette = sns.color_palette("coolwarm", as_cmap=True)
    markers = ['o', '^']
    ax = sns.scatterplot(data=df, x=f'dist_{donor1}', y=f'dist_{donor2}', hue=metal_charge, style=metal, palette=palette, alpha=1, markers=markers)

    norm = plt.Normalize(df[metal_charge].min(), df[metal_charge].max())
    sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
    sm.set_array([])

    # Add a colorbar
    ax.figure.colorbar(sm)

    h, l = ax.get_legend_handles_labels()
    plt.legend(h[6:], l[6:])

    label = 'Distance'
    plt.xlabel(f'{label} to {donor1} in Angstroms')
    plt.ylabel(f'{label} to {donor2} in Angstroms')
    plt.savefig(Path(outdir, f'scatter_metal_charge_dist_{donor1}_{donor2}.svg'))
    plt.close()

    # 2D scatter plots, three in a row. Left: bite angle vs hlgap, centre: bite angle vs metal charge, right: hlgap vs metal charge
    combinations = [(pn_bite_angle, hlgap), (pn_bite_angle, metal_charge), (hlgap, metal_charge)]
    for x, y in combinations:
        plt.figure()
        sns.scatterplot(data=df, x=x, y=y, hue=metal, alpha=0.4)
        plt.xlabel(x)
        plt.ylabel(y)
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig(Path(outdir, f'scatter_{x}_{y}.svg'))
        plt.close()
    plt.figure(figsize=(5, 5))
    plt.subplot(1, 3, 1)
    sns.scatterplot(data=df, x=pn_bite_angle, y=hlgap, hue=metal)
    plt.xlabel('N-M-P Bite angle (°)')
    plt.ylabel('HOMO-LUMO Gap (eV)')
    plt.legend(loc='upper right')
    plt.subplot(1, 3, 2)
    sns.scatterplot(data=df, x=pn_bite_angle, y=metal_charge, hue=metal)
    plt.xlabel('N-M-P Bite angle (°)')
    plt.ylabel('Metal Mulliken Charge (e)')
    plt.legend(loc='upper right')
    plt.subplot(1, 3, 3)
    sns.scatterplot(data=df, x=hlgap, y=metal_charge, hue=metal)
    plt.xlabel('HOMO-LUMO Gap (eV)')
    plt.ylabel('Metal Mulliken Charge (e)')
    plt.legend(loc='upper right')
    plt.tight_layout()
    # plt.show()
    plt.savefig(Path(outdir, f'pairplot_scatter_metal_charge_bite_angle_hlgap.svg'))
    plt.savefig(Path(outdir, f'pairplot_scatter_metal_charge_bite_angle_hlgap.png'), dpi=400, transparent=True)
    plt.close()


    print('Done!')



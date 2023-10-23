import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
sns.set_style("ticks")


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


    print('Done!')



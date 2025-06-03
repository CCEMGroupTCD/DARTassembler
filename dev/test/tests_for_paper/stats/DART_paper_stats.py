from DARTassembler.src.ligand_extraction.DataBase import LigandDB
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

if __name__ == '__main__':

    n_max = None
    plot_dir = 'plots'


    metalig = LigandDB.load_from_json(n_max=n_max)

    #%% Plot a histogram of the ligand geometries
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(exist_ok=True, parents=True)

    df = metalig.get_ligand_output_df()
    df['haptic'] = df['n_haptic_atoms'].apply(lambda n_haptic: 'haptic' if n_haptic > 0 else 'non-haptic')
    # For each denticity, sort first by name and then by occurrence of the geometries
    order = []
    for dent in range(1, 11):
        df_dent = df[df['n_eff_denticities'] == dent]
        order += df_dent['geometry'].value_counts().index.tolist()
    df = df.sort_values(by='geometry', key=lambda x: x.map(lambda y: order.index(y)))

    fig, ax = plt.subplots()
    sns.histplot(data=df, y='geometry', multiple='stack', legend=False, ax=ax)
    plt.xlabel('Count')
    plt.ylabel('Ligand geometry')
    # Rotate the x-axis labels
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/geometry_histogram.svg')

    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/geometry_histogram_log.svg')
    plt.close()


from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.metalig.mol import Ligand
from tqdm import tqdm
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams['svg.fonttype'] = 'none'

def plot_geometry_histogram_as_svg(ligand_dict: dict, output_file: str = "archetype_histogram.svg"):
    """
    Plot a horizontal bar histogram of ligand geometries using total counts only.
    """
    # Sort geometries by total count descending
    sorted_items = sorted(ligand_dict.items(), key=lambda x: x[1]["num_count"], reverse=True)
    geometries = [k for k, _ in sorted_items]
    total_counts = [v["num_count"] for _, v in sorted_items]

    fig, ax = plt.subplots(figsize=(12, 8))

    # Single horizontal bars (no haptic/non-haptic split)
    bars = ax.barh(geometries, total_counts, color='skyblue')

    # Add total count labels at the right end of bars
    for bar, total in zip(bars, total_counts):
        width = bar.get_width()
        ax.text(width + (0.01 * max(total_counts)),  # small offset from bar end
                bar.get_y() + bar.get_height() / 2,
                f'{total}', va='center', ha='left', fontsize=10)

    # Extend x-axis limits for text padding
    ax.set_xlim(0, max(total_counts) * 1.05)  # 5% extra space on the right

    # Move x-axis to the top
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')

    ax.set_xlabel('Count')
    ax.set_ylabel('Geometry')
    ax.set_title('Ligand Geometry Histogram', pad=20)
    ax.invert_yaxis()  # Largest count at the top
    plt.tight_layout()
    plt.savefig(output_file, format='svg')
    plt.close()

def plot_haptic_pie_chart(ligand_dict: dict, geometry: str = "3-meridional", output_file: str = "3-meridional_haptic_pie.svg"):
    """
    Plot a pie chart showing the percentage of 3-meridional ligands that have haptic coordination.
    """

    if geometry not in ligand_dict:
        print(f"No ligands found with geometry: {geometry}")
        return

    # Get counts for 3-meridional ligands
    total = ligand_dict[geometry]["num_count"]
    haptic = ligand_dict[geometry]["haptic_count"]
    non_haptic = ligand_dict[geometry]["non_haptic_count"]

    labels = ['Haptic', 'Non-haptic']
    sizes = [haptic, non_haptic]
    colors = ['#ff9999', '#66b3ff']

    fig, ax = plt.subplots(figsize=(6, 6))
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        colors=colors,
        autopct='%1.1f%%',
        startangle=90,
        counterclock=False,
        wedgeprops={'edgecolor': 'white'}
    )

    # Draw a circle for a donut shape
    centre_circle = plt.Circle((0, 0), 0.70, fc='white')
    fig.gca().add_artist(centre_circle)

    ax.set_title(f'Haptic Coordination in {geometry} Ligands', pad=20)
    plt.tight_layout()
    plt.savefig(output_file, format='svg')
    plt.close()

if __name__ == "__main__":
    db = LigandDB.from_json(n_max=500).db

    ligand_dict = {}
    for name, data in tqdm(db.items(), desc="Processing ligands"):
        data: Ligand
        geometry = data.geometry
        haptic = int(data.n_haptic_atoms)
        donor_atoms = data.donor_elements # a list for example ["N", "O", "S"] or maybe ["N", "O"]
        donor_atoms.sort()  # Sort donor atoms for consistent representation
        donor_atoms_string = "-".join(donor_atoms)
        xyz = data.get_xyz_string()

        ligand_dict[name] = {
            "archetype": geometry,
            "haptic": haptic,
            "donor_atoms": donor_atoms_string,
            "xyz": xyz}

        if geometry not in ligand_dict:
            # Initialize counts
            ligand_dict[geometry] = {
                "num_count": 1,
                "haptic_count": 1 if haptic != 0 else 0,
                "non_haptic_count": 0 if haptic != 0 else 1
            }
        else:
            ligand_dict[geometry]["num_count"] += 1
            if haptic != 0:
                ligand_dict[geometry]["haptic_count"] += 1
            else:
                ligand_dict[geometry]["non_haptic_count"] += 1

    print(f"Total unique geometries: \n{ligand_dict}")
    plot_geometry_histogram_as_svg(ligand_dict, "archetype_histogram.svg")
    plot_haptic_pie_chart(ligand_dict,
                          geometry="3-meridional",
                          output_file="3-meridional_haptic_pie.svg")


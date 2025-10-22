from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.metalig.mol import Ligand
from tqdm import tqdm
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['svg.fonttype'] = 'none'
import numpy as np


# Goals
# 1. show the variety of number of donor atoms for 1-mono
# 2. show a small chart of the number of unique donor atom configs for 2-cis

def filter_ligands_by_geometry(db: LigandDB, geometry: str):
    """
    Filter ligands by their geometry.
    """
    return {k: v for k, v in db.items() if v.geometry == geometry}


def filter_ligands_by_num_haptic_atoms(db: LigandDB):
    """
    Filter ligands by the number of haptic atoms.
    """
    # Step 1: Get unique set of haptic atom counts
    unique_haptic_counts = set(ligand.n_haptic_atoms for ligand in db.values())

    # Step 2: return a list of ligand databases filtered for each unique count
    filtered_ligands = {}
    for count in unique_haptic_counts:
        filtered_ligands[count] = {k: v for k, v in db.items() if v.n_haptic_atoms == count}
    return filtered_ligands


def plot_mini_log_bar(data: dict, output_file="mini_log_bar.svg"):
    """
    PowerPoint-friendly horizontal bar chart (no scaling, no tight layout).
    Bars ordered 0 donors at top. Exports raw SVG.
    """
    # Sort and prepare
    sorted_data = sorted(data.items(), key=lambda x: x[0])
    categories = [f"{k} donors" for k, _ in sorted_data]
    counts = [v for _, v in sorted_data]

    fig, ax = plt.subplots()  # No figsize specified (use default)
    bars = ax.barh(categories, counts, color='skyblue', edgecolor='k')

    # Log scale
    ax.set_xscale('log')

    # Move x-axis to top
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.set_xlabel('Number of Ligands (log scale)')
    ax.set_ylabel('Donors per Ligand')

    # Annotate counts
    for bar, count in zip(bars, counts):
        ax.text(count, bar.get_y() + bar.get_height()/2, str(count),
                va='center', ha='left')

    # Reverse y-axis
    ax.invert_yaxis()

    # Raw save (no layout adjustments)
    plt.savefig(output_file, format='svg')
    plt.close()
    print(f"Saved PowerPoint-safe mini bar chart to {output_file}")

if __name__ == "__main__":
    db = LigandDB.from_json(n_max=5).db

    # Filter ligands by geometry
    monodentate_ligands = filter_ligands_by_geometry(db, "1-mono")

    # Filter by num haptic atoms
    haptic_ligands = filter_ligands_by_num_haptic_atoms(monodentate_ligands)

    print("done")
    plot_mini_log_bar(data={
        0: 9078,
        2: 461,
        3: 148,
        4: 141,
        5: 51,
        6: 320,
        8: 12
    })

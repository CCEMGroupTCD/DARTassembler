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


def filter_ligands_by_donor_atom_config(db: LigandDB):
    """
    Filter ligands by their donor atom configuration.
    :param db: The ligand database.
    :return: a dictionary with key unique donor atom configurations like 'N-N-C'
             and values as lists of Ligand objects.
    """
    donor_atom_configs = {}
    for ligand in db.values():
        donor_atoms = ligand.donor_elements  # e.g., ["N", "O", "S"]
        donor_atoms.sort()  # Sort for consistent representation
        donor_atoms_string = "-".join(donor_atoms)  # e.g., "C-N-N"

        if donor_atoms_string not in donor_atom_configs:
            donor_atom_configs[donor_atoms_string] = []
        donor_atom_configs[donor_atoms_string].append(ligand)

    return donor_atom_configs

def plot_donor_config_pie(donor_config_dict: dict, output_file="donor_config_pie.svg", group_threshold=None):
    """
    Plot a hollow (donut) pie chart showing the number of ligands per donor atom configuration.

    Args:
        donor_config_dict (dict): {config: [Ligand, Ligand, ...]}
        output_file (str): Output SVG file name
        group_threshold (float): Optional threshold (as fraction of total) to group smaller categories into 'Other'
    """
    # Count number of ligands per config
    counts = {k: len(v) for k, v in donor_config_dict.items()}
    total_ligands = sum(counts.values())

    # Optionally group small categories into 'Other'
    if group_threshold:
        grouped_counts = {}
        other_count = 0
        for config, count in counts.items():
            fraction = count / total_ligands
            if fraction < group_threshold:
                other_count += count
            else:
                grouped_counts[config] = count
        if other_count > 0:
            grouped_counts['Other'] = other_count
        counts = grouped_counts

    # Sort configs by count descending
    sorted_items = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    labels, sizes = zip(*sorted_items)

    # Create pie chart
    fig, ax = plt.subplots(figsize=(6, 6))
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        autopct='%1.1f%%',
        startangle=90,
        counterclock=False,
        wedgeprops=dict(width=0.4, edgecolor='w'),  # Hollow donut
        textprops=dict(color="black")
    )

    ax.set_title("Ligand Donor Configurations", fontsize=14, pad=20)
    plt.savefig(output_file, format='svg')
    plt.close()
    print(f"Saved hollow pie chart to {output_file}")

if __name__ == "__main__":
    db = LigandDB.from_json(n_max=500000).db

    # Filter ligands by geometry
    monodentate_ligands = filter_ligands_by_geometry(db, "3-meridional")

    # Filter ligands by donor atom configuration
    donor_atom_configs = filter_ligands_by_donor_atom_config(monodentate_ligands)

    plot_donor_config_pie(donor_atom_configs, output_file="donor_config_pie.svg", group_threshold=0.02)

    print("done")


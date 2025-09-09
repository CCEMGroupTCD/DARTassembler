import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.metalig.mol import Ligand
import numpy as np

# Path to cached CSV
csv_file = "ligand_metal_occurrence_matrix.csv"


def scatter_metals(df, metal_x, metal_y, log_scale=True, color_by_promiscuity=True):
    """
    Debug function to view ligand occurrence scatterplot between two metals.
    Colors ligands by number of unique metals they bind to (promiscuity).

    Parameters:
        df (pd.DataFrame): Ligand-metal occurrence matrix.
        metal_x (str): Metal symbol for X-axis.
        metal_y (str): Metal symbol for Y-axis.
        log_scale (bool): Apply log10 scaling to counts (default True).
        color_by_promiscuity (bool): Color by number of metals bound (default True).
    """
    if metal_x not in df.columns or metal_y not in df.columns:
        print(f"Error: One or both metals ({metal_x}, {metal_y}) not in DataFrame columns.")
        return

    x = df[metal_x].copy()
    y = df[metal_y].copy()

    # Compute promiscuity: number of metals ligand binds to
    promiscuity = (df > 0).sum(axis=1)

    # Apply log scaling if requested
    if log_scale:
        x = x.replace(0, np.nan).apply(np.log10)
        y = y.replace(0, np.nan).apply(np.log10)

    # Drop ligands where either X or Y is NaN
    mask = ~(x.isna() | y.isna())
    x = x[mask]
    y = y[mask]
    promiscuity = promiscuity[mask]

    plt.figure(figsize=(8, 6))

    if color_by_promiscuity:
        norm = colors.Normalize(vmin=promiscuity.min(), vmax=promiscuity.max())
        cmap = cm.viridis
        scatter = plt.scatter(x, y, c=promiscuity, cmap=cmap, norm=norm, alpha=0.8, edgecolors='k', linewidths=0.3)
        cbar = plt.colorbar(scatter)
        cbar.set_label('Number of Metals Bound (Promiscuity)')
    else:
        plt.scatter(x, y, alpha=0.6, edgecolors='k', linewidths=0.5)

    plt.xlabel(f"{metal_x} occurrence" + (" (log₁₀)" if log_scale else ""))
    plt.ylabel(f"{metal_y} occurrence" + (" (log₁₀)" if log_scale else ""))
    plt.title(f"Ligand Occurrence Scatter: {metal_x} vs {metal_y}")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

if os.path.exists(csv_file):
    print(f"Loading ligand-metal matrix from {csv_file}...")
    df = pd.read_csv(csv_file, index_col=0)
else:
    print("CSV not found. Building matrix from LigandDB (this may take a while)...")
    db = LigandDB.from_json(n_max=5000000).db  # or n_max=None for full dataset

    ligand_metal_dict = {}
    for lig_name, lig_data in db.items():
        lig_data: Ligand
        ligand_metal_dict[lig_name] = lig_data.metal_counts

    df = pd.DataFrame.from_dict(ligand_metal_dict, orient='index').fillna(0).astype(int)

    # Save for future runs
    df.to_csv(csv_file)
    print(f"Matrix saved to {csv_file}")

# Inspect
print(f"Ligand-metal matrix shape: {df.shape}")
print(df.head())

# --- Build Metal-Metal Correlation Matrix ---
# Pearson correlation (default)
correlation_matrix = df.corr(method='spearman')

# Optional: Spearman (rank-based)
# correlation_matrix = df.corr(method='spearman')

# Save correlation matrix
corr_csv = "metal_metal_correlation_matrix.csv"
correlation_matrix.to_csv(corr_csv)
print(f"Metal-metal correlation matrix saved to {corr_csv}")

# --- Plot Heatmap ---
import matplotlib.pyplot as plt

# --- Plot Clustered Heatmap ---
g = sns.clustermap(
    correlation_matrix,
    cmap="coolwarm",
    figsize=(12, 12),
    method="ward",          # linkage method
    metric="euclidean",     # distance metric
    cbar_kws={"label": "Correlation (Pearson)"},
)

# Force all tick labels to display
g.ax_heatmap.set_xticks(range(len(correlation_matrix.columns)))
g.ax_heatmap.set_xticklabels(correlation_matrix.columns, rotation=90, ha='center', fontsize=8)
g.ax_heatmap.set_yticks(range(len(correlation_matrix.index)))
g.ax_heatmap.set_yticklabels(correlation_matrix.index, rotation=0, fontsize=8)

plt.title("Clustered Metal-Metal Correlation")
plt.savefig("metal_metal_correlation_clustermap_full_labels.svg")
plt.show()
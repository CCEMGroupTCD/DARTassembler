import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from ast import literal_eval

from DARTassembler.src.ligand_extraction.utilities_Molecule import stoichiometry2atomslist

sns.set(style="whitegrid")

if __name__ == "__main__":

    dirs = ['../DART_output', '../DART_output_targeted']
    outdir = 'data_output'

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    for dir in dirs:
        dir = Path(dir)
        table_path = Path(dir, 'info_table.csv')
        df = pd.read_csv(table_path)
        df = df[df['success'] == True]
        df['ligand stoichiometries'] = df['ligand stoichiometries'].apply(literal_eval)
        atoms = df['ligand stoichiometries'].apply(lambda stois: [stoichiometry2atomslist(stoi) for stoi in stois])
        atoms = atoms.explode().explode()

        plt.figure()
        sns.histplot(atoms)
        plt.xlabel('Elements')


        plt.savefig(Path(outdir, f'histogram_{dir.name}.png'), dpi=300)



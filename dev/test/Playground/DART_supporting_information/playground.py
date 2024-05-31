import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
from DARTassembler.src.constants.Paths import project_path, default_ligand_db_path
from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db
matplotlib.use('TkAgg')
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none' # for correct text rendering in some programs
plt.rcParams['savefig.facecolor'] = 'white'

from pathlib import Path
sns.set_theme(style='ticks')






if __name__ == '__main__':

    db_version = '1.7'
    db_path = project_path().extend(*f'data/final_db_versions/unique_ligand_db_v{db_version}.json'.split('/'))
    nmax = None


    ligands = LigandDB.load_from_json(default_ligand_db_path, n_max=nmax)
    df_ligands = pd.DataFrame.from_dict({name: lig.get_ligand_output_info(add_confident_charge=True) for name, lig in ligands.db.items()}, orient='index')
    metalig_csv_path = Path(default_ligand_db_path).with_suffix('.csv')
    df_metalig = pd.read_csv(metalig_csv_path, index_col=0)


    # Get all unique ligands with same ligand graph
    df_ligands['ligand_graph_wo_metal'] = [lig.graph_hash for lig in ligands.db.values()]
    donors = df_ligands.groupby('ligand_graph_wo_metal')['Donors'].unique()
    # Reduce to ligands with multiple different donors
    donors = donors[donors.apply(lambda x: len(x) > 1)]
    # Reduce to ligands with 3 different donors
    donors3 = donors[donors.apply(lambda x: len(x) == 3)].to_frame().reset_index()
    interesting_graph_hash_wo_metal = '7795ba07bb55fbb9cf75d0407aa80e03'    # chosen from donors3
    interesting_ligands = df_ligands[df_ligands['ligand_graph_wo_metal'] == interesting_graph_hash_wo_metal]







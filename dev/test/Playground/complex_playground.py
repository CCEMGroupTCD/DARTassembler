"""
Playground which reads in the unique ligand db and lets you play with it.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
from DARTassembler.src.constants.Paths import project_path

from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db

matplotlib.use('TkAgg')
from DARTassembler.src.ligand_extraction.DataBase import LigandDB, ComplexDB
sns.set_theme()


if __name__ == '__main__':

    db_version = '1.7'
    db_path = project_path().extend(f'data/final_db_versions/complex_db_v{db_version}.json')

    # nmax = False
    #
    # complexes = ComplexDB.load_from_json(db_path, n_max=nmax)
    # df = complexes.to_dataframe()
    #
    # csvpath = db_path.with_suffix('.csv')
    # print(df.columns)
    # df.to_csv(csvpath, index=False)

    df_complexes = pd.read_csv(db_path.with_suffix('.csv'))
    n_comp_with_mos = (~df_complexes['Metal OS'].isna()).sum()
    print(f'Metals with OS: {n_comp_with_mos}')

    print('Done')




"""
Playground which reads in the unique ligand db and lets you play with it.
"""
import pandas as pd
import seaborn as sns
import matplotlib
from DARTassembler.src.constants.paths import project_path

try:  # Avoid error when running on server
    matplotlib.use('TkAgg')
except ImportError:
    pass
sns.set_theme()


if __name__ == '__main__':

    db_version = '1.7'
    db_path = project_path().extend(f'data/final_db_versions/complex_db_v{db_version}.json')

    # n_max = False
    #
    # complexes = ComplexDB.from_json(db_path, n_max=n_max)
    # df = complexes.to_dataframe()
    #
    # csvpath = db_path.with_suffix('.csv')
    # print(df.columns)
    # df.to_csv(csvpath, index=False)

    df_complexes = pd.read_csv(db_path.with_suffix('.csv'))
    n_comp_with_mos = (~df_complexes['Metal OS'].isna()).sum()
    print(f'Metals with OS: {n_comp_with_mos}')

    print('Done')




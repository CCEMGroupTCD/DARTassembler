"""
This script is only temporarily necessary. It adds necessary global information from the CSD to the input global_properties.csv so that it can be used without reading in other files.
"""
import pandas as pd
from pathlib import Path

if __name__ == '__main__':
    tmQMg_dir = '../database/tmQMg_fixed_gbl_props'


    original_gbl_props = Path(tmQMg_dir, 'global_mol_properties_original_before_fix.csv')
    csd_info_file = Path(tmQMg_dir, 'CSD.csv')

    df = pd.read_csv(original_gbl_props)

    print('Merge global molecular information from CSD into big dataframe.')
    df_CSD = pd.read_csv(csd_info_file)
    df_CSD = df_CSD.drop(columns=['Unnamed: 0'])
    df_CSD = df_CSD.rename(columns={'metal_nr_if_exists': 'metal_oxi_state', 'date': 'CSD_date', 'name': 'CSD_iupac_name', 'metal_name': 'metal'})

    # drop CSD_date temporarily because this is not used in the messy implementation either, to make it comparable.
    df_CSD.drop(columns=['CSD_date', 'metal_str_if_exists'], inplace=True)

    df = pd.merge(df, df_CSD, on='CSD_code')

    out_gbl_props_file = Path(tmQMg_dir, 'global_mol_properties.csv')
    df.to_csv(out_gbl_props_file, index=False)

    print('Done!')




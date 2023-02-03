"""
This script is for generating a dataset for machine learning of the input complexes.
"""
import numpy as np
import pandas as pd
import rdkit
from src01.io_custom import load_unique_ligand_db, load_complex_db, load_full_ligand_db, save_unique_ligand_db, save_full_ligand_db, save_complex_db
from src01.DataBase import ComplexDB
from src01.utilities import unroll_dict_into_columns
from src11_machine_learning.dataset_preparation.descriptors import SOAP_3D

if __name__ == '__main__':

    db_version = '1.3'
    dataset = f'../data/final_db_versions/complex_db_v{db_version}.json'    # all complexes
    # dataset = '../data/tmQMG_Jsons_fixed_gbl_props_cutoffs_test/complex_db_original_connected.json' # only 1000 complexes for testing
    out_csv = f'../database/machine_learning/descriptors_complex_db_v{db_version}.csv'


    print('Load complex database.')
    db = load_complex_db(dataset, 'class')

    # Reduce to only the complexes where we have the os given.
    db = {key: val for key, val in db.items() if not np.isnan(val.metal_oxi_state)}

    ase_complexes = [c.mol for c in db.values()]
    r_cut = 10.0
    n_max= 2
    l_max= 1
    crossover= False
    average= 'inner'
    soap = SOAP_3D(
                    ase_molecules=ase_complexes,
                    r_cut=r_cut,
                    n_max=n_max,
                    l_max=l_max,
                    crossover=crossover,
                    average=average,
                    )
    print('Start creating SOAP descriptors for each complex.')
    soap_desc = soap.calculate_descriptors(only_from_metal=True)
    print('Created SOAP descriptors successfully.')

    db_dict = {name: mol.write_to_mol_dict() for name, mol in db.items()}
    df = pd.DataFrame.from_dict(db_dict, orient='index')
    df = unroll_dict_into_columns(df, dict_col='global_props', prefix='gbl_')

    db_drop_props = ['atomic_props', 'global_props', 'graph_dict', 'ligands']
    df = df.drop(columns=db_drop_props)

    soap_cols = [f'soap_{i}' for i in range(soap_desc.shape[1])]
    soap_desc = pd.DataFrame(data=soap_desc, columns=soap_cols, index=df.index)
    df = df.join(soap_desc, validate='1:1')

    print(f'Save csv with complexes, properties and descriptors to {out_csv}.')
    df.to_csv(out_csv, index=False)

    print('Done!')









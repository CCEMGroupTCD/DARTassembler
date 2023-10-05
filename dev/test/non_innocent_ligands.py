import numpy as np

from DARTassembler.src.constants.Paths import default_ligand_db_path
from DARTassembler.src.ligand_extraction.io_custom import iterate_unique_ligand_db, iterate_complex_db
from DARTassembler.src.constants.Paths import project_path
import pandas as pd

if __name__ == '__main__':
    complex_db_path = project_path().extend(*'data/final_db_versions/complex_db_v1.7.json'.split('/'))
    nmax = False

    all_ligands = {}
    for c_id, c in iterate_complex_db(complex_db_path, n_max=nmax):
        metal_os = c['metal_oxi_state']
        if not np.isnan(metal_os):
            ligands = c['ligands']
            confident_charges = [lig['pred_charge_is_confident'] for lig in ligands]

            if sum(confident_charges) == len(ligands) -1:   # if only one ligand charge is missing and all others known
                nilig = [lig for conf, lig in zip(confident_charges, ligands) if not conf]
                assert len(nilig) == 1, confident_charges
                nilig = nilig[0]

                if nilig['denticity'] > 0:  # skip unconnected ligands
                    sum_lig_charges = sum(lig['pred_charge'] for lig in ligands if lig['pred_charge_is_confident'])
                    nilig['inno_charge'] = c['charge'] - metal_os - sum_lig_charges
                    assert not np.isnan(nilig['inno_charge']), (c['charge'], metal_os, sum_lig_charges)
                    all_ligands[nilig['name']] = nilig
    df_all_ligands = pd.DataFrame.from_dict(all_ligands, orient='index')
    df_all_ligands['inno_charge'] = df_all_ligands['inno_charge'].astype(int)

    #%% Reduce all ligands to unique ligands with all different charges and other properties
    unique_props = ['n_atoms',  # properties which I know that they are the same for all ligands with same unique_name
       'n_hydrogens', 'n_protons', 'graph_hash', 'n_bonds',
       'heavy_atoms_graph_hash', 'stoichiometry',
       'graph_hash_with_metal',
       'heavy_atoms_graph_hash_with_metal', 'has_betaH',
       'has_neighboring_coordinating_atoms','denticity', 'unique_name',
       'pred_charge', 'pred_charge_is_confident']
    changing_props = ['has_bond_order_attribute', 'has_unknown_bond_orders',
       'has_good_bond_orders', 'bond_order_graph_hash', 'original_complex_id', 'local_elements', 'was_connected_to_metal', 'original_metal',
       'original_metal_position', 'original_metal_symbol', 'original_metal_os',
       'is_centrosymmetric', 'centrosymmetry_ang_dev', 'ligand_to_metal']
    df_ligands = df_all_ligands.groupby('unique_name')[['inno_charge'] + unique_props].agg(list)
    df_ligands['inno_charge'] = df_ligands['inno_charge'].apply(np.unique)
    df_ligands['n_charges'] = df_ligands['inno_charge'].apply(len)
    df_only_ni_ligands = df_ligands.loc[df_ligands['n_charges'] > 1]

    for col in unique_props:
        df_only_ni_ligands[col] = df_only_ni_ligands[col].apply(lambda x: x[0])
    ni_ligands = df_only_ni_ligands.to_dict(orient='index')
    df_only_ni_ligands = df_only_ni_ligands.to_dict(orient='index')
    for uname, row in ni_ligands.items():
        charges = row['inno_charge']
        props = {prop: [] for prop in changing_props}
        for charge in charges:
            ligs = df_all_ligands.loc[df_all_ligands['unique_name'] == uname].to_dict(orient='index')
            lig = [lig for lig in ligs.values() if lig['inno_charge'] == charge][0]
            for prop in changing_props:
                props[prop].append(lig[prop])
        df_only_ni_ligands[uname].update(props)
    df_only_ni_ligands = pd.DataFrame.from_dict(df_only_ni_ligands, orient='index')



    n_ni_ligands = len(df_only_ni_ligands)
    n_all_ni_ligands = df_only_ni_ligands['n_charges'].sum()
    print(f'Num. redox active ligands with diff. charges: {n_all_ni_ligands}\nNum. different redox active ligands: {n_ni_ligands}')
    print('Done!')












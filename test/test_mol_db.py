"""
This script is for testing the molecule databases. Currently it just reads them in.
"""
from src01.DataBase import MoleculeDB, LigandDB, ComplexDB
# from src02_Pre_Assembly_Filtering.FilteringStage import FilterStage

if __name__ == '__main__':
    n_test = False
    db_version = '1.5'
    unique_ligands_path = f'../data/final_db_versions/unique_ligand_db_v{db_version}.json'
    complex_db_path = f'../data/final_db_versions/complex_db_v{db_version}.json'


    # Ligands
    db_ulig = LigandDB.from_json(
                            json_=unique_ligands_path,
                            type_='Ligand',
                            max_number=n_test
                            )
    lig = db_ulig.get_first_entry()

    lig_props = list(lig.__dict__.keys())
    lig_should_have = ['stoichiometry', 'atomic_props', 'global_props', 'graph_dict', 'denticity', 'ligand_to_metal', 'local_elements', 'name', 'CSD_code', 'graph_hash', 'unique_name', 'occurrences', 'count_denticities', 'count_metals', 'n_denticities', 'n_metals', 'chosen_denticity_fraction', 'all_ligand_names', 'pred_charge', 'pred_charge_is_confident']
    diff = set(lig_should_have).symmetric_difference(lig_props)
    print('Difference of properties:', diff)

    # Molecules
    db_mol = MoleculeDB.from_json(
                                    json_=unique_ligands_path,
                                    type_='Molecule',
                                    max_number=n_test
                                    )
    mol = db_mol.get_first_entry()

    # Complexes
    db_complex = ComplexDB.from_json(
                                    json_=complex_db_path,
                                    type_='Complex',
                                    max_number=n_test
                                    )
    complex = db_complex.get_first_entry()
    complexes_should_have = ['atomic_props', 'global_props', 'graph_dict', 'mol_id', 'stoichiometry',
       'metal', 'metal_oxi_state', 'charge', 'ligands.py', 'graph_hash']
    complex_props = list(complex.__dict__.keys())
    diff = set(complexes_should_have).symmetric_difference(complex_props)
    print('Difference of properties:', diff)


    # Filter = FilterStage(db_ulig)
    #
    # print(f"Before Filtering: {len(Filter.database.db)}")
    #
    # lig = db_ulig.get_first_entry()
    #
    # lig_props = list(lig.__dict__.keys())
    # lig_should_have = ['stoichiometry', 'atomic_props', 'global_props', 'graph_dict', 'denticity', 'ligand_to_metal', 'local_elements', 'name', 'CSD_code', 'graph_hash', 'unique_name', 'occurrences', 'count_denticities', 'count_metals', 'n_denticities', 'n_metals', 'chosen_denticity_fraction', 'all_ligand_names', 'pred_charge', 'pred_charge_is_confident']
    # diff = set(lig_should_have).symmetric_difference(lig_props)
    # print('Difference of properties:', diff)
    #
    # # Molecules
    # db_mol = MoleculeDB.from_json(
    #                                 json_=unique_ligands_path,
    #                                 type_='Molecule',
    #                                 max_number=n_test
    #                                 )
    # mol = db_mol.get_first_entry()
    #
    # # Complexes
    # db_complex = ComplexDB.from_json(
    #                                 json_=complex_db_path,
    #                                 type_='Complex',
    #                                 max_number=n_test
    #                                 )
    # complex = db_complex.get_first_entry()
    # complexes_should_have = ['atomic_props', 'global_props', 'graph_dict', 'mol_id', 'stoichiometry',
    #    'metal', 'metal_oxi_state', 'charge', 'ligands', 'graph_hash']
    # complex_props = list(complex.__dict__.keys())
    # diff = set(complexes_should_have).symmetric_difference(complex_props)
    # print('Difference of properties:', diff)



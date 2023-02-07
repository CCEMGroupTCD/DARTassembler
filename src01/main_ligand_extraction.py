"""
This is the main script for the extraction of ligands from a database.
"""
from src01.ligand_extraction import LigandExtraction

def main(database_path: str, data_store_path: str, calculate_charges: bool=True, overwrite_atomic_properties: bool=True, use_existing_input_json: bool=True, exclude_not_fully_connected_complexes: bool=True ,get_only_unique_ligand_db_without_charges: bool=False):

    db = LigandExtraction(
                            database_path=database_path,
                            data_store_path=data_store_path,
                            exclude_not_fully_connected_complexes=exclude_not_fully_connected_complexes
                            )
    db.run_ligand_extraction(
                                calculate_charges=calculate_charges,
                                overwrite_atomic_properties=overwrite_atomic_properties,
                                use_existing_input_json=use_existing_input_json,
                                get_only_unique_ligand_db_without_charges=get_only_unique_ligand_db_without_charges
                            )

    return

if __name__ == '__main__':

    database_path = '../database/tmQMg_fixed_gbl_props_cutoffs'             # in github
    data_store_path = "../data/tmQMG_Jsons_fixed_gbl_props_cutoffs_full"    # directory where we want to store the jsons
    calculate_charges = True                                                # if you want to run charge assignment after ligand extraction, takes ~30 min on tmQMg
    overwrite_atomic_properties = True                                      # if atomic properties json should be overwritten, not really critical
    use_existing_input_json = False                                          # if the existing input json should be used or the process started from the xzy files
    exclude_not_fully_connected_complexes = True                            # script not ready for unconnected graphs yet
    get_only_unique_ligand_db_without_charges = False                       # For graph benchmark useful, reduces runtime because it ignores charge assignment and updating the complex and full ligand db.

    main(
            database_path=database_path,
            data_store_path=data_store_path,
            calculate_charges=calculate_charges,
            overwrite_atomic_properties=overwrite_atomic_properties,
            use_existing_input_json=use_existing_input_json,
            exclude_not_fully_connected_complexes=exclude_not_fully_connected_complexes,
            get_only_unique_ligand_db_without_charges=get_only_unique_ligand_db_without_charges
        )
    
"""
This is the main script for the extraction of ligands from a database.
"""
from src01.ligand_extraction import LigandExtraction


def main(database_path_: str,
         data_store_path_: str,
         calculate_charges_: bool = True,
         overwrite_atomic_properties_: bool = True,
         use_existing_input_json_: bool = True,
         exclude_not_fully_connected_complexes_: bool = True,
         get_only_unique_ligand_db_without_charges_: bool = False
         ):

    db = LigandExtraction(
                            database_path=database_path_,
                            data_store_path=data_store_path_,
                            exclude_not_fully_connected_complexes=exclude_not_fully_connected_complexes_,
                            testing=100
                            )
    db.run_ligand_extraction(
                                calculate_charges=calculate_charges_,
                                overwrite_atomic_properties=overwrite_atomic_properties_,
                                use_existing_input_json=use_existing_input_json_,
                                get_only_unique_ligand_db_without_charges=get_only_unique_ligand_db_without_charges_
                            )

    return


if __name__ == '__main__':

    # some example database paths
    tmqm_path = "../data_input/tmQM"            # directory of the complex database in correct input format with xyz files and global properties for tmQM
    tmqmg_path = "../data_input/tmQMG"          # directory of the complex database in correct input format with xyz files and global properties for tmQMG
    csd_G_MM_path = "../data_input/CSD_MM_G"    # directory of the complex database in correct input format with xyz files and global properties for MM from CSD with Graphs

    data_store_path = "../data/tmQMG_Jsons_fixed_gbl_props_cutoffs_full"    # directory where we want to store the jsons
    calculate_charges = True                                                # if you want to run charge assignment after ligand extraction, takes ~30 min on tmQMg
    overwrite_atomic_properties = True                                      # if atomic properties json should be overwritten, not really critical
    use_existing_input_json = True                                         # if the existing input json should be used or the process started from the xzy files
    exclude_not_fully_connected_complexes = True                            # script not ready for unconnected graphs yet
    get_only_unique_ligand_db_without_charges = False                       # For graph benchmark useful, reduces runtime because it ignores charge assignment and updating the complex and full ligand db.

    main(
            database_path_=tmqm_path,
            data_store_path_=data_store_path,
            calculate_charges_=calculate_charges,
            overwrite_atomic_properties_=overwrite_atomic_properties,
            use_existing_input_json_=use_existing_input_json,
            exclude_not_fully_connected_complexes_=exclude_not_fully_connected_complexes,
            get_only_unique_ligand_db_without_charges_=get_only_unique_ligand_db_without_charges
        )

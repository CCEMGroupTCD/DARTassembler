"""
This is the main script for the extraction of ligands from a database.
"""
from src01.ligand_extraction import LigandExtraction
from typing import Union


def select_example_database(DB: str) -> (str, str):
    """
    some example database paths
    returns:
        directory of the complex database in correct input format with xyz files and global properties
        directory where we want to store the jsons
    """

    if DB.lower() == "tmqm":
        return "../data_input/tmQM", "../data_output/tmQM_Jsons"
    elif DB.lower() == "tmqmg":
        return "../data_input/tmQMG", "../data_output/tmQMG_Jsons"
    elif DB.lower() == "csd_mm_g":
        return "../data_input/CSD_MM_G", "../data_output/CSD_MM_G_Jsons"
    else:
        # unknown DB
        return "", ""


def main(database_path_: str,
         data_store_path_: str,
         calculate_charges_: bool = True,
         overwrite_atomic_properties_: bool = True,
         use_existing_input_json_: bool = True,
         exclude_not_fully_connected_complexes_: bool = True,
         get_only_unique_ligand_db_without_charges_: bool = False,
         testing_: Union[bool, int] = False
         ):
    db = LigandExtraction(
        database_path=database_path_,
        data_store_path=data_store_path_,
        exclude_not_fully_connected_complexes=exclude_not_fully_connected_complexes_,
        testing=testing_
    )

    db.run_ligand_extraction(
        calculate_charges=calculate_charges_,
        overwrite_atomic_properties=overwrite_atomic_properties_,
        use_existing_input_json=use_existing_input_json_,
        get_only_unique_ligand_db_without_charges=get_only_unique_ligand_db_without_charges_
    )

    return


if __name__ == '__main__':

    #
    #
    database_path, data_store_path = select_example_database(DB="CSD_MM_G")

    #
    testing = 500           # if we would like to only do a test run

    calculate_charges = True  # if you want to run charge assignment after ligand extraction, takes ~30 min on tmQMg
    overwrite_atomic_properties = True  # if atomic properties json should be overwritten, not really critical
    use_existing_input_json = True  # if the existing input json should be used or the process started from the xzy files
    exclude_not_fully_connected_complexes = False  # script not ready for unconnected graphs yet
    get_only_unique_ligand_db_without_charges = True  # For graph benchmark useful, reduces runtime because it ignores charge assignment and updating the complex and full ligand db.

    main(
        database_path_=database_path,
        data_store_path_=data_store_path,
        calculate_charges_=calculate_charges,
        overwrite_atomic_properties_=overwrite_atomic_properties,
        use_existing_input_json_=use_existing_input_json,
        exclude_not_fully_connected_complexes_=exclude_not_fully_connected_complexes,
        get_only_unique_ligand_db_without_charges_=get_only_unique_ligand_db_without_charges,
        testing_=testing
    )

"""
This is the main script for the extraction of ligands.py from a database.
"""
from src01.ligand_extraction import LigandExtraction
from typing import Union
from constants.constants import project_path


def select_example_database(DB: str) -> (str, str):
    """
    some example database paths
    returns:
        directory of the complex database in correct input format with xyz files and global properties
        directory where we want to store the jsons
    """

    if DB.lower() == "tmqm":
        return f"{project_path}/data_input/tmQM", f"{project_path}/data_output/tmQM_Jsons"
    elif DB.lower() == "tmqmg":
        return f"{project_path}/data_input/tmQMG", f"{project_path}/data_output/tmQMG_Jsons"
    elif DB.lower() == "csd_mm":
        return f"{project_path}/data_input/CSD_MM", f"{project_path}/data_output/CSD_MM_Jsons"
    elif DB.lower() == "csd_mm_g":
        return f"{project_path}/data_input/CSD_MM_G", f"{project_path}/data_output/CSD_MM_G_Jsons"
    else:
         # unknown DB
        raise ValueError(f'Database not recognized: {DB}')


def main(database_path_: str,
         data_store_path_: str,
         calculate_charges_: bool = True,
         overwrite_atomic_properties_: bool = True,
         use_existing_input_json_: bool = True,
         exclude_not_fully_connected_complexes_: bool = True,
         testing_: Union[bool, int] = False,
         graph_strat_: str = "default",
         exclude_charged_complexes: bool = False,
         max_charge_iterations: Union[int, None] = 10,
         **kwargs
         ):

    db = LigandExtraction(
        database_path=database_path_,
        data_store_path=data_store_path_,
        exclude_not_fully_connected_complexes=exclude_not_fully_connected_complexes_,
        testing=testing_,
        graph_strat=graph_strat_,
        exclude_charged_complexes=exclude_charged_complexes
    )

    db.run_ligand_extraction(
        calculate_charges=calculate_charges_,
        overwrite_atomic_properties=overwrite_atomic_properties_,
        use_existing_input_json=use_existing_input_json_,
        max_charge_iterations=max_charge_iterations,
        **kwargs
    )

    return db


if __name__ == '__main__':

    # example databases, choose between: tmqm, tmqmG, CSD_MM_G
    database = "CSD_MM_G"

    testing = 50_000           # if we would like to only do a test run. Set to False for full run
    graph_strategy = "default"  # the desired graph strategy: default, ase_cutoff, CSD, pymatgen_NN, molsimplifyGraphs

    calculate_charges = True  # if you want to run charge assignment after ligand extraction, takes ~30 min on tmQMg
    overwrite_atomic_properties = True  # if atomic properties json should be overwritten. Only necessary after changing input files.
    use_existing_input_json = False  # if the existing input json should be used or the process started from the xzy files


    # Input complex filters
    exclude_not_fully_connected_complexes = True  # only keep complexes which are fully connected
    exclude_charged_complexes = True   # Keep only input complexes with charge of 0



    database_path, data_store_path = select_example_database(DB=database)
    db = main(
        database_path_=database_path,
        data_store_path_=data_store_path,
        calculate_charges_=calculate_charges,
        overwrite_atomic_properties_=overwrite_atomic_properties,
        use_existing_input_json_=use_existing_input_json,
        exclude_not_fully_connected_complexes_=exclude_not_fully_connected_complexes,
        testing_=testing,
        graph_strat_=graph_strategy,
        exclude_charged_complexes=exclude_charged_complexes,
    )

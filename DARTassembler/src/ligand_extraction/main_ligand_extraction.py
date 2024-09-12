"""
This is the main script for the extraction of ligands from a database.
"""
from DARTassembler.src.ligand_extraction.ligand_extraction import LigandExtraction
from typing import Union
from DARTassembler.src.constants.Paths import project_path

def main(database_path_: str,
         data_store_path_: str,
         overwrite_atomic_properties_: bool = True,
         use_existing_input_json_: bool = True,
         exclude_not_fully_connected_complexes_: bool = True,
         testing_: Union[bool, int] = False,
         graph_strat_: str = "default",
         exclude_charged_complexes: bool = False,
         max_charge_iterations: Union[int, None] = 10,
         store_database_in_memory: bool = False,
         **kwargs
         ):

    db = LigandExtraction(
        database_path=database_path_,
        data_store_path=data_store_path_,
        exclude_not_fully_connected_complexes=exclude_not_fully_connected_complexes_,
        testing=testing_,
        graph_strat=graph_strat_,
        exclude_charged_complexes=exclude_charged_complexes,
        store_database_in_memory=store_database_in_memory
    )

    db.run_ligand_extraction(
        overwrite_atomic_properties=overwrite_atomic_properties_,
        use_existing_input_json=use_existing_input_json_,
        max_charge_iterations=max_charge_iterations,
        **kwargs
    )

    return db


if __name__ == '__main__':

    # example databases, choose between: tmqm, tmqmG, CSD
    database_path = project_path().extend(*'data_input/CSD'.split('/'))
    data_store_path = project_path().extend(*'data_output/CSD'.split('/'))

    testing = False           # if we would like to only do a test run. Set to False for full run
    graph_strategy = "CSD"  # the desired graph strategy: default, ase_cutoff, CSD

    overwrite_atomic_properties = False  # if atomic properties json should be overwritten. Only necessary after changing input files.
    use_existing_input_json = False  # if the existing input json should be used or the process started from the xzy files
    store_database_in_memory = False    # if the database should be stored in memory. Only use if you have enough RAM, but can speed up the pipeline by maybe 30%.

    # Input complex filters
    exclude_not_fully_connected_complexes = False  # only keep complexes which are fully connected
    exclude_charged_complexes = False   # Keep only input complexes with charge of 0


    db = main(
        database_path_=database_path,
        data_store_path_=data_store_path,
        overwrite_atomic_properties_=overwrite_atomic_properties,
        use_existing_input_json_=use_existing_input_json,
        exclude_not_fully_connected_complexes_=exclude_not_fully_connected_complexes,
        testing_=testing,
        graph_strat_=graph_strategy,
        exclude_charged_complexes=exclude_charged_complexes,
        store_database_in_memory=store_database_in_memory
    )

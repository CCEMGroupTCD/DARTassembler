"""
Building up on the new main_ligand_extraction.py
"""

from src01.main_ligand_extraction import main as run
from src01.main_ligand_extraction import select_example_database
from src01.DataBase import LigandDB, MoleculeDB


if __name__ == "__main__":

    #
    # example databases, choose between: tmqm, tmqmG, CSD_MM_G
    database_path, data_store_path = select_example_database(DB="tmQMG")

    #
    # testing = 1000  # if we would like to only do a test run (only works from the second run on)
    graph_strategy = "default"  # the desired graph strategy: default, ase_cutoff, CSD, pymatgen_NN, molsimplifyGraphs

    #
    #
    # calculate_charges = True  # if you want to run charge assignment after ligand extraction, takes ~30 min on tmQMg
    # overwrite_atomic_properties = True  # if atomic properties json should be overwritten, not really critical
    # use_existing_input_json = True  # if the existing input json should be used or the process started from the xzy files
    # exclude_not_fully_connected_complexes = False  # script not ready for unconnected graphs yet
    # get_only_unique_ligand_db_without_charges = True  # For graph benchmark useful, reduces runtime because it ignores charge assignment and updating the complex and full ligand db.

    run(
        database_path_=database_path,
        data_store_path_=data_store_path,
        calculate_charges_=False,
        overwrite_atomic_properties_=False,
        use_existing_input_json_=False,
        exclude_not_fully_connected_complexes_=False,
        get_only_unique_ligand_db_without_charges_=True,
        testing_=False,
        graph_strat_=graph_strategy
    )

    # Read in the DB

    Complex_DB = MoleculeDB.from_json(json_=f"{data_store_path}/complex_db.json",
                                      type_="Molecule",
                                      max_number=100
                                      )

    Ligand_DB = LigandDB.from_json(
        json_=f"{data_store_path}/tmQM_Ligands_full.json",
        type_="Ligand",
        max_number=100
    )

    unique_Ligands = LigandDB.from_json(
        json_=f"{data_store_path}/tmQM_Ligands_unique.json",
        type_="Ligand",
        max_number=100
    )
"""
Building up on the new main_ligand_extraction.py,
is going to become the graph testing thing lateron
"""
import networkx as nx

from src01.main_ligand_extraction import main as run
from src01.main_ligand_extraction import select_example_database
from src01.DataBase import LigandDB, MoleculeDB


if __name__ == "__main__":


    # testing = 1000  # if we would like to only do a test run (only works from the second run on)
    graph_strategy = "CSD"  # the desired graph strategy: default, ase_cutoff, CSD, pymatgen_NN, molsimplifyGraphs
    selected_DB = "CSD_MM_G"

    if selected_DB.lower() in ["tmqm", "tmqmg"] and graph_strategy == "CSD":
        raise ValueError("CSD Graphs can not be used with tmqm atomic properties")

    #
    database_path, data_store_path = select_example_database(DB=selected_DB)
    #
    run(
        database_path_=database_path,
        data_store_path_=data_store_path,
        calculate_charges_=False,
        overwrite_atomic_properties_=False,
        use_existing_input_json_=False,
        exclude_not_fully_connected_complexes_=False,
        get_only_unique_ligand_db_without_charges_=True,
        testing_=100,
        graph_strat_=graph_strategy
    )

    # Read in the DB

    # Complex_DB = MoleculeDB.from_json(json_=f"{data_store_path}/complex_db.json",
    #                                   type_="Molecule",
    #                                   max_number=100
    #                                   )
    #
    # # number of unconnected complexes
    # metric01 = len([mol for mol in Complex_DB.db.values() if nx.is_connected(mol.graph) is False])
    #
    # #
    # Ligand_DB = LigandDB.from_json(
    #     json_=f"{data_store_path}/tmQM_Ligands_full.json",
    #     type_="Ligand",
    #     max_number=100
    # )
    # #
    # metric02 = len(Ligand_DB.db)                                                        # number of total ligands
    # metric03 = len([lig for lig in Ligand_DB.db.values() if lig.denticity == -1])       # number of isolated ligands
    #
    # #
    # unique_Ligands = LigandDB.from_json(
    #     json_=f"{data_store_path}/tmQM_Ligands_unique.json",
    #     type_="Ligand",
    #     max_number=100
    # )
    #
    # metric04 = sum([len(ul.count_denticities) for ul in unique_Ligands.db.values()])        # number of unique ligands wirh denticity

    #
    print("done")

"""
I always needed a file where, which could simply be run in the debugger
mode to get quick acess to the tmQM and its ligands.py
in the RCA_Molecule and RCA_Ligand format respecitvely
and here it is
"""
from DARTassembler.src.ligand_extraction.DataBase import MoleculeDB
from src01.DataLoader import DataLoader


def csd_graphs_not_usable_yet():
    """
    The problem is, that we atm use the tmQM 3D coordinates as an input
    but use the graphs of the uncurated CSD entries to compare them against
    but this means, that the total graphs of the CSD will be bigger,
    as the unconnected counter ions or solution molecules are still flowing around
    and hence, as the graph has more nodes than the list of 3D coordinates, this will reaise
    an error, as this script shows.
    """
    tmQM_DB = MoleculeDB.from_json(json_=DataLoader(database_path_=database_path).data_for_molDB,
                                   type_="Molecule",
                                   max_number=100,
                                   )

    tmQM_DB2 = MoleculeDB.from_json(json_=DataLoader(database_path_=database_path).data_for_molDB,
                                    type_="Molecule",
                                    max_number=100,
                                    graph_strategy="CSD"
                                    )

    b_mol = list(tmQM_DB2.db.values())[1]

    for i, (a, b) in enumerate(zip(tmQM_DB.db.values(), tmQM_DB2.db.values())):

        a_els = [a.graph.nodes[node]['node_label'] for node in a.graph.nodes]
        b_els = [b.graph.nodes[node]['node_label'] for node in b.graph.nodes]

        if a_els != b_els:
            breakpoint()


if __name__ == "__main__":
    """
    data_path = "../data_output/tmQMG_Jsons"

   

    #
    tmQM_DB = MoleculeDB.from_json(json_=f'{data_path}/tmQMG.json', type_="Molecule", n_max=["NIBTAT"])

    # Create the LigandDB from the tmQM
    tmQM_Ligands = LigandDB.from_json(json_=f'{data_path}/tmQMG_Ligands_full.json', type_="Ligand", n_max=["CSD-NIBTAT-03-a"])

    lig = list(tmQM_Ligands.db.values())[0]
    mol = list(tmQM_DB.db.values())[0]

    tmQM_unique_Ligands = LigandDB.from_json(json_='../data/New_DB_jsons/tmQM_ligands_unique.json', type_="Ligand")



    #with open(f"{data_store_path}/tmQMG_Ligands_full.json") as file:
    #    full_ligands_dict = json.load(file)

    print("Playground established")

    print("done")
    """
    database_path = '../data_input/tmQM'

    DB = MoleculeDB.from_json(json_="../../data_output/CSD_MM_G_Jsons/tmQMG.json",
                              type_="Molecule",
                              max_number=100,
                              graph_strategy="CSD"
                              )

    tmQM_DB = MoleculeDB.from_json(json_=DataLoader(database_path_=database_path).data_for_molDB,
                                   type_="Molecule",
                                   max_number=100,
                                   graph_strategy="CSD"
                                   )

    tmQM_DB2 = MoleculeDB.from_json(json_=DataLoader(database_path_=database_path).data_for_molDB,
                                    type_="Molecule",
                                    max_number=100,
                                    graph_strategy="CSD"
                                    )

    """

    a = list(tmQM_DB.db.values())[0]

    tmQM_Ligands = LigandDB.from_MoleculeDB(molDB=tmQM_DB,
                                                 denticity_numbers_of_interest=[-1, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                                                10]
                                                 )
    """
    print("done")

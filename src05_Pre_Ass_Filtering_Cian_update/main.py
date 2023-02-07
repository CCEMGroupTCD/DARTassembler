from src01.DataBase import LigandDB
from src02_Pre_Ass_Filtering_Cian.FilteringStage import FilterStage

if __name__ == "__main__":
    """
    One easy method which includes all the custom filters one could potentially apply to a ligand
    filtered_db = {unique_lig_name: RCA_Ligand} dict-type, which can then be employed for the ligandAssembl
    """

    tmQM_unique_Ligands = LigandDB.from_json(json_='/Users/cianclarke/Documents/PhD/Complex_Assembly/Data/unique_ligand_db_v1.3.json',
                                             type_="Ligand",
                                             max_number=1000
                                             )

    Filter = FilterStage(tmQM_unique_Ligands)

    # print(f"Before Filtering: {len(Filter.database.db)}")

    # Filter.metals_of_interest_filter(metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"])

    # print(f"After metal of interest filter: {len(Filter.database.db)}")

    Filter.denticity_of_interest_filter(denticity_of_interest=[2, 3, 4, 5])

    print(f"After denticity filter: {len(Filter.database.db)}")

    # Filter.filter_functional_group_atoms(atoms_of_interest=["N", "O"])

    # print(f"After NO filter: {len(Filter.database.db)}")

    # Filter.filter_betaHs()

    # print(f"After beta H filter: {len(Filter.database.db)}")

    Filter.box_excluder_filter()

    print(f"After box_excluder_filter: {len(Filter.database.db)}")

    Filter.add_constant_ligands()

    Filter.database.to_json(path="../data/Filtered_Jsons/filteredLigDB_Cian_Ct_ff_graph.json")

    print("Filtering Done")

    #todo: filter out ligands with coordinating atoms --> "iff contains these" "doesnt contain these" etc...
    #todo: filter out ligands that have a denticity occurrence below a certain percentage -> chosen denticity fraction
    #todo: filter by molecular weight,
    #todo: filter by atoms that must be in the molecule
    #todo: filter by ligands that dont contain particular atoms
    #todo: filter out ligands with 5 atoms in a plane
    #todo: subgraph match
    #todo: each filter should be a sub filter of the denticity filters




from src01.DataBase import LigandDB
from src02_Pre_Assembly_Filtering.FilteringStage import FilterStage
from src01.main_ligand_extraction import select_example_database

# ["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"]


if __name__ == "__main__":
    """
    One easy method which includes all the custom filters one could potentially apply to a ligand
    filtered_db = {unique_lig_name: RCA_Ligand} dict-type, which can then be employed for the ligandAssembl
    """

    #
    selected_db = "CSD_MM_G"        # from which DB the ligands should come
    _, datapath = select_example_database(DB=selected_db)

    tmQM_unique_Ligands = LigandDB.from_json(json_=f"{datapath}/tmQM_Ligands_unique.json",
                                             type_="Ligand",
                                             max_number=5000
                                             )

    Filter = FilterStage(tmQM_unique_Ligands)


    print(f"Before Filtering: {len(Filter.database.db)}")

    #
    #
    # Filter: Denticity of Interest

    Filter.denticity_of_interest_filter(denticity_of_interest=[2, 3, 4, 5])

    print(f"After Denticity Filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Metal of Interest

    Filter.metals_of_interest_filter(denticity=2, metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"])

    Filter.metals_of_interest_filter(denticity=3, metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"])

    Filter.metals_of_interest_filter(denticity=4, metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"])

    Filter.metals_of_interest_filter(denticity=5, metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"])

    print(f"After metal of interest filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Coordinating Atoms

    Filter.filter_coordinating_group_atoms(denticity=2, atoms_of_interest=["C", "H", "B"], instruction="must_exclude")

    Filter.filter_coordinating_group_atoms(denticity=3, atoms_of_interest=["C", "H", "B"], instruction="must_exclude")

    Filter.filter_coordinating_group_atoms(denticity=4, atoms_of_interest=["C", "H", "B"], instruction="must_exclude")

    Filter.filter_coordinating_group_atoms(denticity=5, atoms_of_interest=["C", "H", "B"], instruction="must_exclude")

    print(f"After Coordinating Group filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Ligand Atoms

    Filter.filter_ligand_atoms(denticity=2, atoms_of_interest=['C', 'H'], instruction="must_at_least_contain")

    Filter.filter_ligand_atoms(denticity=3, atoms_of_interest=['C', 'H'], instruction="must_at_least_contain")

    Filter.filter_ligand_atoms(denticity=4, atoms_of_interest=['C', 'H'], instruction="must_at_least_contain")

    Filter.filter_ligand_atoms(denticity=5, atoms_of_interest=['C', 'H'], instruction="must_at_least_contain")

    print(f"After Ligand Atom filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Beta Hydrogen

    # Filter.filter_betaHs()

    # print(f"After beta H filter: {len(Filter.database.db)}")


    #
    #
    # Filter: Molecular Weight

    # Filter.filter_molecular_weight(denticity=2, atomic_weight_min=200.1, atomic_weight_max=300.2)

    # Filter.filter_molecular_weight(denticity=3, atomic_weight_min=200.1, atomic_weight_max=300.2)

    # Filter.filter_molecular_weight(denticity=4, atomic_weight_min=200.1, atomic_weight_max=300.2)

    # Filter.filter_molecular_weight(denticity=5, atomic_weight_min=200.1, atomic_weight_max=300.2)

    # print(f"After Molecular_Weight_filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Denticity Fraction
    # 0.9 as a default
    Filter.filter_denticity_fraction(denticity=2, fraction=0.90)

    Filter.filter_denticity_fraction(denticity=3, fraction=0.90)

    Filter.filter_denticity_fraction(denticity=4, fraction=0.90)

    Filter.filter_denticity_fraction(denticity=5, fraction=0.90)

    print(f"After Denticity Fraction: {len(Filter.database.db)}")

    #
    #
    # Filter: Box Excluder

    Filter.box_excluder_filter()

    print(f"After Box Excluder_filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Add Constant Ligands

    Filter.add_constant_ligands()

    # provide an overwrite option?
    Filter.database.to_json(path="../data/Filtered_Jsons/filteredLigDB_refactor_test_030223.json")

    print("Filtering Done")

    # todo: filter out ligands with 5 atoms in a plane
    # todo: subgraph match


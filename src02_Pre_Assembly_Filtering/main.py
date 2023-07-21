from src01.DataBase import LigandDB
from constants.Paths import project_path, default_ligand_db_path
from src02_Pre_Assembly_Filtering.FilteringStage import FilterStage
from test.Integration_Test import IntegrationTest
from pathlib import Path

# ["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"]
lanthanides_actinides = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]



if __name__ == "__main__":

    ligand_db_path = default_ligand_db_path
    outpath = project_path().extend("src14_Assembly_Unit_Test", "ligandfilters", "filtered_ligand_db_v1.7.json")
    """
    One easy method which includes all the custom filters one could potentially apply to a ligand
    filtered_db = {unique_lig_name: RCA_Ligand} dict-type, which can then be employed for the ligandAssembl
    """
    #'/Users/cianclarke/Documents/PhD/Complex_Assembly/Data/unique_ligand_db_v1.6.json'
    #TODO: Timo be careful I have stuck my ligand json somewhere else I dont know if we have
    # a dedicated place for it in the code but we probably should
    # that's why I have no refactored this path
    tmQM_unique_Ligands = LigandDB.from_json(json_=ligand_db_path,
                                             type_="Ligand",
                                             max_number=1000
                                             )

    Filter = FilterStage(tmQM_unique_Ligands)

    print(f"Before Filtering: {len(Filter.database.db)}")

    #
    #
    # Filter: Denticity of Interest

    Filter.denticity_of_interest_filter(denticity_of_interest=[1, 2, 3, 4, 5])

    print(f"After Denticity Filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Neighbouring Coordinating Atoms

    Filter.filter_neighbouring_coordinating_atoms()

    print(f"After Neighbouring Coordinating Atoms: {len(Filter.database.db)}")

    #
    #
    # Filter: Charge Confidence

    Filter.filter_charge_confidence(filter_for="confident")

    print(f"After Charge Confidence Filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Metal of Interest

    # Filter.metals_of_interest_filter(denticity=2, metals_of_interest=["Pd", "Ni"])

    Filter.metals_of_interest_filter(denticity=3, metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"])

    Filter.metals_of_interest_filter(denticity=4, metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"])

    # Filter.metals_of_interest_filter(denticity=5, metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"])

    # print(f"After metal of interest filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Coordinating Atoms

    Filter.filter_coordinating_group_atoms(denticity=3, atoms_of_interest=["N", "N", "N"], instruction="must_contain_and_only_contain")

    # Filter.filter_coordinating_group_atoms(denticity=2, atoms_of_interest=["O", "O"], instruction="must_contain_and_only_contain")

    # Filter.filter_coordinating_group_atoms(denticity=4, atoms_of_interest=["N", "O", "P"], instruction="must_only_contain_in_any_amount")

    # Filter.filter_coordinating_group_atoms(denticity=5, atoms_of_interest=["P", "N", "O", "S"], instruction="must_only_contain_in_any_amount")

    print(f"After Coordinating Group filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Ligand Atoms

    # Filter.filter_ligand_atoms(denticity=2, atoms_of_interest=['C', 'H', 'N', 'O'], instruction="must_only_contain_in_any_amount")

    Filter.filter_ligand_atoms(denticity=2, atoms_of_interest=['C', 'H'], instruction="must_at_least_contain")

    # Filter.filter_ligand_atoms(denticity=2, atoms_of_interest=['Si', 'B'], instruction="must_exclude")

    # Filter.filter_ligand_atoms(denticity=3, atoms_of_interest=['C', 'H'], instruction="must_at_least_contain")

    # Filter.filter_ligand_atoms(denticity=3, atoms_of_interest=['Si', 'B'], instruction="must_exclude")

    # Filter.filter_ligand_atoms(denticity=4, atoms_of_interest=['C', 'H'], instruction="must_at_least_contain")

    # Filter.filter_ligand_atoms(denticity=4, atoms_of_interest=['Si', 'B'], instruction="must_exclude")

    # Filter.filter_ligand_atoms(denticity=5, atoms_of_interest=['C', 'H'], instruction="must_at_least_contain")

    # Filter.filter_ligand_atoms(denticity=5, atoms_of_interest=['Si', 'B'], instruction="must_exclude")

    print(f"After Ligand Atom filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Beta Hydrogen

    Filter.filter_betaHs()

    print(f"After beta H filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Molecular Weight

    # Filter.filter_molecular_weight(denticity=2, atomic_weight_min=1.1, atomic_weight_max=300.2)
    #
    Filter.filter_molecular_weight(denticity=3, atomic_weight_min=30.1, atomic_weight_max=300.2)

    Filter.filter_molecular_weight(denticity=4, atomic_weight_min=80.1, atomic_weight_max=300.2)

    Filter.filter_molecular_weight(denticity=5, atomic_weight_min=80.1, atomic_weight_max=300.2)

    print(f"After Molecular_Weight_filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Denticity Fraction

    # Filter.filter_denticity_fraction(denticity=2, fraction=0.90)

    # Filter.filter_denticity_fraction(denticity=3, fraction=0.90)

    # Filter.filter_denticity_fraction(denticity=4, fraction=0.90)

    # Filter.filter_denticity_fraction(denticity=5, fraction=0.90)

    # print(f"After Denticity Fraction: {len(Filter.database.db)}")

    #
    #
    # Filter: Box Excluder

    # Filter.box_excluder_filter()

    print(f"After Box Excluder_filter: {len(Filter.database.db)}")

    #
    #
    # Filter: even / odd electron filter

    Filter.filter_even_odd_electron(filter_for="even")

    print(f"After even / odd electron filter: {len(Filter.database.db)}")

    #
    #
    # Filter: filter ligand_charges

    Filter.filter_ligand_charges(denticity=3, charge=0)

    Filter.filter_ligand_charges(denticity=2, charge=0)
    #
    # Filter.filter_ligand_charges(denticity=1, charge=-1)

    print(f"After charge filter: {len(Filter.database.db)}")

    #
    #
    # Filter: filter monodentate ligands with lines of symmetry that pass through the coordinating atom

    # Filter.filter_symmetric_monodentate_ligands(instruction="Add", threshold=8)

    print(f"After symmetric filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Atom_Count

    Filter.filter_atom_count(denticity=1, number=10, instruction="less_than")

    print(f"After monodentate size filter: {len(Filter.database.db)}")

    #
    #
    # Filter: Add Constant Ligands

    # Filter.add_constant_ligands()   # not working atm because of CSD_code not in global props

    # Filter.filter_sub_structure_search(denticity=2, SMARTS="c1ncccc1", instruction="must_include") # not implemented yet


    Filter.database.to_json(path=str(outpath), json_lines=True)
    print(f"END: {len(Filter.database.db)}")
    print("Filtering Done")

    # todo: subgraph match


    test = IntegrationTest(new_dir=outpath.parent, old_dir=Path(outpath.parent.parent, 'ligandfilters_benchmark'))
    test.compare_all()
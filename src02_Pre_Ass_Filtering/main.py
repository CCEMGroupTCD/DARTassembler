from src01.DataBase import LigandDB
from src02_Pre_Ass_Filtering.FilteringStage import FilterStage

# only for testing
from src02_Pre_Ass_Filtering.test import id_list


if __name__ == "__main__":
    """
    One easy method which includes all the custom filters one could potentially apply to a ligand
    filtered_db = {unique_lig_name: RCA_Ligand} dict-type, which can then be employed for the ligandAssembl
    """

    tmQM_unique_Ligands = LigandDB.from_json(json_='../data/tmQMG_Jsons/tmQMG_Ligands_unique.json',
                                             type_="Ligand",
                                             # only for testing
                                             # identifier_list=id_list
                                             )

    Filter = FilterStage(tmQM_unique_Ligands)

    Filter.metals_of_interest_filter(metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"])

    Filter.denticity_of_interest_filter(denticity_of_interest=[2, 3, 4, 5])

    Filter.filter_functional_group_atoms(atoms_of_interest=["N", "O"])

    Filter.filter_betaHs()

    Filter.add_constant_ligands()

    Filter.database.to_json(path="../data/Filtered_Jsons/filteredLigDB.json")

    #F = LigandDB.from_json(json_="../src02_Pre_Ass_Filtering/filteredLigDB.json", type_="Ligand")

    print("Filtering Done")

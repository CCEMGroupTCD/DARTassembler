from src.process import LigandDatabase
import pickle
from copy import deepcopy


class FilterHandler:

    def __init__(self, database: LigandDatabase):
        """
        We modify the filtered database dict and the filtered database object we initialized in the ligand databse earlier
        the dict will keep track of all intermediate ligand_dicts
        and the filtered database will be the quick accesss to the ligand_dict after all filters applied
        """
        self.database = database

        # add track of applied filters to db object
        # by adding an attribute to the ligand DB
        self.database.filtered_database_dict = {}        # {1 : [Name of the applied filter, database after the step],
                                                        # 2 : [Name of the applied filter, db after step 2 AND 1, ...}

        # keeps track of the current filtered db, so that we can apply all filter consecutively
        # add this as an attribute of the DB as the filterhandler is exclusively working on that object
        self.database.filtered_database = deepcopy(self.database.full_ligand_dict)

    def filter_duplicates(self):
        """
        In this step we filter duplicates out of the DB
        """

        # generate a new attribute for the database, i.e. expand it
        duplicant_filtered_ligand_dict = {}

        for denticity, ligand_list in self.database.filtered_database.items():
            duplicant_filtered_ligand_dict[denticity] = list(set(ligand_list))

        key_ = len(self.database.filtered_database_dict.keys())
        self.database.filtered_database_dict[key_] = ["Duplicant Filter", duplicant_filtered_ligand_dict]

        # update the filtered ligand db
        self.database.filtered_database = duplicant_filtered_ligand_dict

    def filter_N_and_O_functional_groups(self):
        """
        Only leave in ligands where we have functional atoms equal to N and/or O
        """

        # generate a new attribute for the database, i.e. expand it
        no_filtered_ligand_dict = {}

        for denticity, ligand_list in self.database.filtered_database.items():
            no_filtered_ligand_dict[denticity] = [lig for lig in ligand_list if lig.NO_check() is True]

        key_ = len(self.database.filtered_database_dict.keys())
        self.database.filtered_database_dict[key_] = ["NO Filter", no_filtered_ligand_dict]

        # update the filtered ligand db
        self.database.filtered_database = no_filtered_ligand_dict

    def filter_betaHs(self):
        """
        Filter out all ligands with beta Hydrogen in it
        """

        # generate a new attribute for the database, i.e. expand it
        betaH_filtered_ligand_dict = {}

        for denticity, ligand_list in self.database.filtered_database.items():
            betaH_filtered_ligand_dict[denticity] = [lig for lig in ligand_list if lig.betaH_check() is False]

        key_ = len(self.database.filtered_database_dict.keys())
        self.database.filtered_database_dict[key_] = ["beta H Filter", betaH_filtered_ligand_dict]

        # update the filtered ligand db
        self.database.filtered_database = betaH_filtered_ligand_dict

    def safe(self, safe_path: str):
        """
        print results and safe them to local pickle files
        """

        for key, item in self.database.filtered_database_dict.items():
            print(f"In the {key} step, filter {item[0]} was applied")

        with open(safe_path, "wb") as file:
            pickle.dump(self.database, file)


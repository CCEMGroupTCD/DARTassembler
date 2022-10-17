from src.process import LigandDatabase
import pickle


class FilterHandler:

    def __init__(self, database: LigandDatabase):
        self.full_db = database

        # add track of applied filters to db object
        self.full_db.filters_applied = {}           # In order: {1:First filter,...}

    def filter_duplicates(self, safe_path=None, i=1):
        """
        :param safe_path: Storepath for the new LigandDB object, if safe is desired
        :param i: The position of this filter in the filtering process.
        """

        # generate a new attribute for the database, i.e. expand it
        self.full_db.duplicant_filtered_ligand_dict = {}

        for denticity, ligand_list in self.full_db.full_ligand_dict.items():
            self.full_db.duplicant_filtered_ligand_dict[denticity] = list(set(ligand_list))

        self.full_db.filters_applied[i] = "Duplicant Filter Applied"

        if safe_path is not None:
            with open(safe_path, "wb") as handle:
                pickle.dump(self.full_db, handle)

    def filter_N_and_O_functional_groups(self, safe_path=None, i=1):
        """
        :param safe_path: Storepath for the new LigandDB object, if safe is desired
        :param i: The position of this filter in the filtering process.
        """

        # generate a new attribute for the database, i.e. expand it
        self.full_db.NO_filtered_ligand_dict = {}

        for denticity, ligand_list in self.full_db.full_ligand_dict.items():

            self.full_db.NO_filtered_ligand_dict[denticity] = [lig for lig in ligand_list if lig.NO_check() is True]

        self.full_db.filters_applied[i] = "NO Filter Applied"

        if safe_path is not None:
            with open(safe_path, "wb") as handle:
                pickle.dump(self.full_db, handle)

    def filter_betaHs(self, safe_path=None, i=1):
        """
        :param safe_path: Storepath for the new LigandDB object, if safe is desired
        :param i: The position of this filter in the filtering process.
        """

        # generate a new attribute for the database, i.e. expand it
        self.full_db.NO_filtered_ligand_dict = {}

        for denticity, ligand_list in self.full_db.full_ligand_dict.items():

            self.full_db.NO_filtered_ligand_dict[denticity] = [lig for lig in ligand_list if lig.betaH_check() is False]

        self.full_db.filters_applied[i] = "NO Filter Applied"

        if safe_path is not None:
            with open(safe_path, "wb") as handle:
                pickle.dump(self.full_db, handle)
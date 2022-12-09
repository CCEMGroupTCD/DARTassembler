from pathlib import Path
from copy import deepcopy

from src01.DataBase import MoleculeDB, LigandDB
from src01.utilities_Molecule import original_metal_ligand

from src02_Pre_Ass_Filtering.Box_Excluder_Filter import box_filter
from src02_Pre_Ass_Filtering.constant_Ligands import get_monodentate_list, get_reactant


class FilterStage:

    def __init__(self,
                 database: [MoleculeDB, LigandDB],
                 safe_path: [str, Path] = None):
        """
        :param database:
        :param safe_path:   The path we want to safe the filtered databases to. If None, no saving at all
        """
        # The object we are working on. Is hopefully saved, because the filter stage doesnt make copy itself of it
        self.database = database

        self.safe_path = safe_path
        self.safe = True if self.safe_path is not None else False
        self.filter_tracking = {}

    def safe_and_document_after_filterstep(self, filtering_step_name: str):

        print(f"Filtering step {filtering_step_name} succesfull applied")

        # Create a new attr for database for backtracking
        self.database.filters_applied = self.filter_tracking

        # safe to desired folder
        self.database.to_json(path=f"{self.safe_path}/DB_after_{filtering_step_name}.json")

    def metals_of_interest_filter(self, metals_of_interest: [str, list[str]] = None):
        """

        """
        print("Metal of Interest Filter running")

        if metals_of_interest is None:
            print("No metals of interest selected, no filtering at all")
        elif isinstance(metals_of_interest, str):
            metals_of_interest = [metals_of_interest]

        new_db = deepcopy(self.database.db)

        for identifier, ligand in self.database.db.items():
            om = original_metal_ligand(ligand)
            if not (om in metals_of_interest or om is None):
                del new_db[identifier]

        self.database.db = new_db

        #self.db = {identifier: ligand for identifier, ligand in self.db.items() if om in metals_of_interest or om is None}

        self.filter_tracking[len(self.filter_tracking)] = f"Metals of interest filter with {metals_of_interest}"

        #self.safe_and_document_after_filterstep(filtering_step_name="Metal_of_Interest")

    def denticity_of_interest_filter(self, denticity_of_interest: [int, list[int]] = None):
        """

        """
        print("Denticity Filter running")

        if denticity_of_interest is None:
            print("No denticities of interest selected, no filtering at all")
        elif isinstance(denticity_of_interest, int):
            denticity_of_interest = [denticity_of_interest]

        new_db = deepcopy(self.database.db)

        for identifier, ligand in self.database.db.items():
            if not ligand.denticity in denticity_of_interest:
                del new_db[identifier]

        self.database.db = new_db

        # old: self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items() if ligand.denticity in denticity_of_interest}

        self.filter_tracking[len(self.filter_tracking)] = f"Metals of interest filter with {denticity_of_interest}"

        #self.safe_and_document_after_filterstep(filtering_step_name="Denticity_of_Interest")

    def filter_functional_group_atoms(self, atoms_of_interest: [str, list[str]] = None):
        """
        Only leave in ligands where we have functional atoms equal to N and/or O
        """
        print("FunctionalGroup Filter running")

        self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items()
                   if ligand.functional_atom_check(atoms_of_interest) is True}

        self.filter_tracking[len(self.filter_tracking)] = f"Functional Atom filter with {atoms_of_interest}"

        #self.safe_and_document_after_filterstep(filtering_step_name="FunctionalAtoms_of_Interest")

    def filter_betaHs(self):
        """
        Filter out all ligands with beta Hydrogen in it
        """
        print("betaH Filter running")

        self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items()
                   if ligand.betaH_check() is False}

        self.filter_tracking[len(self.filter_tracking)] = f"betaH Filter"

        #self.safe_and_document_after_filterstep(filtering_step_name="betaH Filter")

    def box_excluder_filter(self):
        """
        Filter out all ligands that violate the box Filter
        """
        print("Box Excluder Filter running")

        self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items()
                   if box_filter(ligand) is True}

        self.filter_tracking[len(self.filter_tracking)] = f"Box Filter"

        #self.safe_and_document_after_filterstep(filtering_step_name="Box Filter")

    def add_constant_ligands(self):
        """
        Now we add the constant Ligands we defined in constant_Ligands
        """
        print("Adding constant Ligands")

        for lig in get_monodentate_list() + get_reactant():
            self.database.db[lig.name] = lig

        #self.db.to_json(path=f"{self.safe_path}/DB_after_Adding_Const_Ligands.json")
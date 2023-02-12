from pathlib import Path
from copy import deepcopy

from src01.DataBase import MoleculeDB, LigandDB

from src02_Pre_Assembly_Filtering.box_filter import box_filter
from src02_Pre_Assembly_Filtering.constant_Ligands import get_monodentate_list, get_reactant



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


    def metals_of_interest_filter(self, denticity: int = None, metals_of_interest: [str, list[str]] = None):
        """
        This function has been updated by Cian to work with the updated .mol file
        """
        print(f"Filtering: Metal of Interest -> Denticity: {denticity}")
        if (metals_of_interest is None) or (denticity is None):
            print("!!!Warning!!! -> No metal of interest or denticity selected -> Proceeding to next filter")
            return
        elif isinstance(metals_of_interest, str):
            metals_of_interest = [metals_of_interest]
        new_db = deepcopy(self.database.db)

        for unq_name, ligand in self.database.db.items():
            # Iterate through each ligand. If the denticity matches
            # the one specified by the user {denticity: int}, but it was never part of
            # complex with a metal specified by the user {metals_of_interest: [str, list[str]]} then it is
            # removed
            if int(ligand.denticity) == int(denticity):
                #print(list(ligand.count_metals))
                metal_of_interest_present = False
                for metal in metals_of_interest:
                    if metal in list(ligand.count_metals):
                        #print("Matching Metal: " + str(metal))
                        metal_of_interest_present = True
                    else:
                        pass
                if metal_of_interest_present:
                    pass
                    #print("Metal of Interest Pass")
                else:
                    del new_db[unq_name]
                    #print("Metal of Interest Fail -> Deleting ligand")
            else:
                pass

        self.database.db = new_db
        # self.db = {identifier: ligand for identifier, ligand in self.db.items() if om in metals_of_interest or om is None}
        self.filter_tracking[len(self.filter_tracking)] = f"Metals of interest filter with {metals_of_interest}"
        # self.safe_and_document_after_filterstep(filtering_step_name="Metal_of_Interest")


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


        # self.safe_and_document_after_filterstep(filtering_step_name="Denticity_of_Interest")

    def filter_coordinating_group_atoms(self, denticity: int = None, atoms_of_interest: [str, list[str]] = None, instruction: str = None):
        # instruction = must_contain_and_only_contain,this means that if the coordinating atoms specified by the user are exactly the same as the coordinating groups of the ligand then the ligand can pass
        # instruction = must_at_least_contain        ,this means that if the coordinating atoms specified by the user are a subset of that of the ligand then the ligand can pass
        # instruction = must_exclude                 ,this means that if the coordinating atoms specified by the user are not contained in any amount in the ligand then the ligand can pass
        # instruction = must_only_contain_in_any_amount   ,this means that if all coordinating atoms specified by the user are contained to some degree in the ligand and no other coordinating atoms are conatined then the ligand can pass
        # This filter only applies to ligands of the specified denticities. ligands with other denticities are allowed to pass
        """
        Only leave in ligands where we have functional atoms equal to specified atoms
        """
        print("Functional Group Filter running")

        print(f"Filtering: Coordinating_atom_type -> Denticity: {denticity}")
        if (atoms_of_interest is None) or (denticity is None) or (instruction is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")
            return
        if ((denticity != len(atoms_of_interest)) and instruction == True) or ((denticity < len(atoms_of_interest)) and instruction == False):
            print("!!!Warning!!! -> The specified denticity is not consistent with specified coordinating atoms -> Proceeding to next filter ")

        elif isinstance(atoms_of_interest, str):
            atoms_of_interest = [atoms_of_interest]

        new_db = deepcopy(self.database.db)

        for unq_name, ligand in self.database.db.items():
            #print("\n")
            #print(unq_name)
            #print(sorted(list(ligand.local_elements)))
            if int(ligand.denticity) == int(denticity):
                coordinating_atoms_present = False
                # If the denticity of the ligand matches that specified by the user
                if ((sorted(list(ligand.local_elements)) == sorted(atoms_of_interest)) and instruction == "must_contain_and_only_contain") or \
                        (all(elem in list(ligand.local_elements) for elem in atoms_of_interest) and instruction == "must_at_least_contain") or \
                        ((any(elem in list(ligand.local_elements) for elem in atoms_of_interest) == False) and instruction == "must_exclude") or \
                        ((all(elem in atoms_of_interest for elem in list(ligand.local_elements))) and instruction == "must_only_contain_in_any_amount"):
                    coordinating_atoms_present = True
                else:
                    pass
                if coordinating_atoms_present:
                    #print("Matching Coordinating groups PASS")
                    pass
                else:
                    #print("Matching Coordinating groups Fail")
                    del new_db[unq_name]
            else:
                # If the denticities don't match then we don't apply the filter and just skip
                pass
        self.database.db = new_db
        self.filter_tracking[len(self.filter_tracking)] = f"Functional Atom filter with {atoms_of_interest}"

        # self.safe_and_document_after_filterstep(filtering_step_name="FunctionalAtoms_of_Interest")

    def filter_ligand_atoms(self, denticity: int = None, atoms_of_interest: [str, list[str]] = None, instruction: str = None):
        # instruction = must_contain_and_only_contain,this means that if the coordinating atoms specified by the user are exactly the same as the coordinating groups of the ligand then the ligand can pass
        # instruction = must_at_least_contain        ,this means that if the coordinating atoms specified by the user are a subset of that of the ligand then the ligand can pass
        # instruction = must_exclude                 ,this means that if the coordinating atoms specified by the user are not contained in any amount in the ligand then the ligand can pass
        # instruction = must_only_contain_in_any_amount   ,this means that if all coordinating atoms specified by the user are contained to some degree in the ligand and no other coordinating atoms are conatined then the ligand can pass
        # This filter only applies to ligands of the specified denticities. ligands with other denticities are allowed to pass
        """
        Only leave in ligands where we have functional atoms equal to specified atoms
        """
        print("FunctionalGroup Filter running")

        print(f"Filtering: ligand_atom_type")
        if (atoms_of_interest is None) or (instruction is None) or (denticity is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")
            return

        elif isinstance(atoms_of_interest, str):
            atoms_of_interest = [atoms_of_interest]

        new_db = deepcopy(self.database.db)

        for unq_name, ligand in self.database.db.items():
            #print("\n")
            #print(unq_name)
            #print(sorted(list(ligand.atomic_props["atoms"])))
            if int(ligand.denticity) == int(denticity):
                coordinating_atoms_present = False
                # If the denticity of the ligand matches that specified by the user
                if ((sorted(list(ligand.atomic_props["atoms"])) == sorted(atoms_of_interest)) and instruction == "must_contain_and_only_contain") or \
                        (all(elem in list(ligand.atomic_props["atoms"]) for elem in atoms_of_interest) and instruction == "must_at_least_contain") or \
                        ((any(elem in list(ligand.atomic_props["atoms"]) for elem in atoms_of_interest) == False) and instruction == "must_exclude") or \
                        ((all(elem in atoms_of_interest for elem in list(ligand.atomic_props["atoms"]))) and instruction == "must_only_contain_in_any_amount"):
                    coordinating_atoms_present = True
                else:
                    pass
                if coordinating_atoms_present:
                    #print("Matching Coordinating groups PASS")
                    pass
                else:
                    #print("Matching Coordinating groups Fail")
                    del new_db[unq_name]
            else:
                # If the denticities don't match then we don't apply the filter and just skip
                pass
        self.database.db = new_db
        self.filter_tracking[len(self.filter_tracking)] = f"Functional Atom filter with {atoms_of_interest}"

    # self.safe_and_document_after_filterstep(filtering_step_name="FunctionalAtoms_of_Interest")


    def filter_betaHs(self):
        """
        Filter out all ligands with beta Hydrogen in it
        """

        #print("betaH Filter running")

        self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items()
                            if ligand.betaH_check() is False}

        self.filter_tracking[len(self.filter_tracking)] = f"betaH Filter"

        # self.safe_and_document_after_filterstep(filtering_step_name="betaH Filter")

    def filter_molecular_weight(self, denticity: int = None, atomic_weight_min: float = None, atomic_weight_max: float = None):

        new_db = deepcopy(self.database.db)

        if (denticity is None) or (atomic_weight_min is None) or (atomic_weight_max is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")

        else:
            for unq_name, ligand in self.database.db.items():
                #print(unq_name)
                molecular_range_condition_satisified = False
                MW = ligand.global_props["molecular_weight"]
                print(MW)
                if (atomic_weight_min <= MW < atomic_weight_max) and (ligand.denticity == denticity):
                    pass
                    #print("Molecular Weight Pass")
                else:
                    del new_db[unq_name]
                    #print("Molecular Weight Fail")
                #print("\n")
        self.database.db = new_db
        self.filter_tracking[len(self.filter_tracking)] = f"Molecular Weight Filter with MW_MIN_{atomic_weight_min} and MW_MAX_{atomic_weight_max}"

    def filter_denticity_fraction(self, denticity: int = None, fraction: float = None):
        # This filter will filter out ligands whose dominating denticity does not account for greater than the proportion
        # specified by the user {fraction: float = None} of the total occurences of the complex
        new_db = deepcopy(self.database.db)
        if (denticity is None) or (fraction is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")

        else:
            for unq_name, ligand in self.database.db.items():
                #print(unq_name)
                occurence = float(ligand.occurrences)
                #print(occurence)
                current_highest_occurrence = 0
                #print(ligand.count_denticities)
                for key, value in ligand.count_denticities.items():
                    if value > current_highest_occurrence:
                        current_highest_occurrence = value
                    else:
                        pass
                occurence_fraction = float(current_highest_occurrence/occurence)
                #print(occurence_fraction)
                if occurence_fraction > fraction:
                    pass
                    #print("PASS")
                else:
                    del new_db[unq_name]
                    #print("FAIL")
                #print("\n")
            self.database.db = new_db
            self.filter_tracking[len(self.filter_tracking)] = f"Denticity Fraction: {fraction}"


    def box_excluder_filter(self):
        """
        Filter out all ligands that violate the box Filter
        """
        print("Box Excluder Filter running")
        self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items() if box_filter(ligand) is True}
        self.filter_tracking[len(self.filter_tracking)] = f"Box Filter"


        # self.safe_and_document_after_filterstep(filtering_step_name="Box Filter")


    def add_constant_ligands(self):
        """
        Now we add the constant Ligands we defined in constant_Ligands
        """
        print("Adding constant Ligands")

        for lig in get_monodentate_list() + get_reactant():
            self.database.db[lig.name] = lig


        # self.db.to_json(path=f"{self.safe_path}/DB_after_Adding_Const_Ligands.json")


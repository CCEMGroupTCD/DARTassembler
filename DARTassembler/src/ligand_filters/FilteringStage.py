import warnings
from pathlib import Path
from copy import deepcopy
import numpy as np

from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.DataBase import MoleculeDB, LigandDB
from typing import Union

from DARTassembler.src.ligand_extraction.Molecule import pseudo_metal
from DARTassembler.src.ligand_filters.constant_Ligands import get_monodentate_list
from DARTassembler.src.ligand_extraction.utilities_Molecule import has_smarts_pattern


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

    def stoichiometry_filter(self, stoichiometry, denticities: Union[list[int], None] = None):
        to_delete = []
        for unq_name, ligand in self.database.db.items():
            if denticities is None or ligand.denticity in denticities:
                if not ligand.if_same_stoichiometry(stoichiometry):
                    to_delete.append(unq_name)
        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}

    def metals_of_interest_filter(self, metals_of_interest: [str, list[str]], denticity: int = None):
        """
        This function has been updated by Cian to work with the updated .mol file
        """
        # print(f"Filtering: Metal of Interest -> Denticity: {denticity}")
        denticity = self.ensure_denticities_is_list(denticity)

        if isinstance(metals_of_interest, str):
            metals_of_interest = [metals_of_interest]

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            # Iterate through each ligand. If the denticity matches
            # the one specified by the user {denticity: int}, but it was never part of
            # complex with a metal specified by the user {metals_of_interest: [str, list[str]]} then it is
            # removed
            if ligand.denticity in denticity:
                # print(list(ligand.count_metals))
                metal_of_interest_present = False
                for metal in metals_of_interest:
                    if metal in list(ligand.count_metals):
                        # print("Matching Metal: " + str(metal))
                        metal_of_interest_present = True
                if not metal_of_interest_present:
                    to_delete.append(unq_name)
                    # print("Metal of Interest Fail -> Deleting ligand")

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}

        # self.db = {identifier: ligand for identifier, ligand in self.db.items() if om in metals_of_interest or om is None}
        self.filter_tracking[len(self.filter_tracking)] = f"Metals of interest filter with {metals_of_interest}"
        # self.safe_and_document_after_filterstep(filtering_step_name="Metal_of_Interest")

    def graph_hash_with_metal_filter(self, graph_hashes_with_metal: list):
        to_delete = []
        for unq_name, ligand in self.database.db.items():
            if ligand.graph_hash_with_metal not in graph_hashes_with_metal:
                to_delete.append(unq_name)

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Graph IDs: {graph_hashes_with_metal}"

        return


    def denticity_of_interest_filter(self, denticity_of_interest: [int, list[int]]):
        """
        """
        # print("Denticity Filter running")
        if isinstance(denticity_of_interest, int):
            denticity_of_interest = [denticity_of_interest]

        to_delete = []
        for identifier, ligand in self.database.db.items():
            if not ligand.denticity in denticity_of_interest:
                to_delete.append(identifier)

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}

        # old: self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items() if ligand.denticity in denticity_of_interest}

        self.filter_tracking[len(self.filter_tracking)] = f"Metals of interest filter with {denticity_of_interest}"

        # self.safe_and_document_after_filterstep(filtering_step_name="Denticity_of_Interest")

    def filter_unconnected_ligands(self):
        to_delete = []
        for identifier, ligand in self.database.db.items():
            if ligand.denticity <= 0:
                to_delete.append(identifier)

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}

        self.filter_tracking[len(self.filter_tracking)] = f"Unconnected Ligand Filter"

        return

    def filter_coordinating_group_atoms(self, atoms_of_interest: [str, list[str]], instruction: str, denticity: int = None, ):
        # instruction = must_contain_and_only_contain,this means that if the coordinating atoms specified by the user are exactly the same as the coordinating groups of the ligand then the ligand can pass
        # instruction = must_at_least_contain        ,this means that if the coordinating atoms specified by the user are a subset of that of the ligand then the ligand can pass
        # instruction = must_exclude                 ,this means that if the coordinating atoms specified by the user are not contained in any amount in the ligand then the ligand can pass
        # instruction = must_only_contain_in_any_amount   ,this means that if all coordinating atoms specified by the user are contained to some degree in the ligand and no other coordinating atoms are conatined then the ligand can pass
        # This filter only applies to ligands of the specified denticities. ligands with other denticities are allowed to pass
        """
        Only leave in ligands where we have functional atoms equal to specified atoms
        """
        # print("Functional Group Filter running")

        # print(f"Filtering: Coordinating_atom_type -> Denticity: {denticity}")
        denticity = self.ensure_denticities_is_list(denticity)

        if isinstance(atoms_of_interest, str):
            atoms_of_interest = [atoms_of_interest]

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            # print("\n")
            # print(unq_name)
            # print(sorted(list(ligand.local_elements)))
            if ligand.denticity in denticity:
                coordinating_atoms_present = False
                # If the denticity of the ligand matches that specified by the user
                if ((sorted(list(ligand.local_elements)) == sorted(atoms_of_interest)) and instruction == "must_contain_and_only_contain") or \
                        (all(elem in list(ligand.local_elements) for elem in atoms_of_interest) and instruction == "must_at_least_contain") or \
                        ((any(elem in list(ligand.local_elements) for elem in atoms_of_interest) == False) and instruction == "must_exclude") or \
                        ((all(elem in atoms_of_interest for elem in list(ligand.local_elements))) and instruction == "must_only_contain_in_any_amount"):
                    coordinating_atoms_present = True

                if not coordinating_atoms_present:
                    # print("Matching Coordinating groups Fail")
                    to_delete.append(unq_name)

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Functional Atom filter with {atoms_of_interest}"

        # self.safe_and_document_after_filterstep(filtering_step_name="FunctionalAtoms_of_Interest")

    def filter_ligand_atoms(self, atoms_of_interest: Union[str, list[str]], instruction: str, denticity: int = None):
        # instruction = must_contain_and_only_contain,this means that if the coordinating atoms specified by the user are exactly the same as the coordinating groups of the ligand then the ligand can pass
        # instruction = must_at_least_contain, this means that if the coordinating atoms specified by the user are a subset of that of the ligand then the ligand can pass
        # instruction = must_exclude, this means that if the coordinating atoms specified by the user are not contained in any amount in the ligand then the ligand can pass
        # instruction = must_only_contain_in_any_amount   ,this means that if all coordinating atoms specified by the user are contained to some degree in the ligand and no other coordinating atoms are conatined then the ligand can pass
        # This filter only applies to ligands of the specified denticities. ligands with other denticities are allowed to pass
        """
        Only leave in ligands where we have functional atoms equal to specified atoms
        """
        denticity = self.ensure_denticities_is_list(denticity)

        # print("FunctionalGroup Filter running")

        # print(f"Filtering: ligand_atom_type")

        if isinstance(atoms_of_interest, str):
            atoms_of_interest = [atoms_of_interest]

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            # print("\n")
            # print(unq_name)
            # print(sorted(list(ligand.atomic_props["atoms"])))
            if ligand.denticity in denticity:
                coordinating_atoms_present = False
                # If the denticity of the ligand matches that specified by the user
                if ((sorted(list(ligand.atomic_props["atoms"])) == sorted(atoms_of_interest)) and instruction == "must_contain_and_only_contain") or \
                        (all(elem in list(ligand.atomic_props["atoms"]) for elem in atoms_of_interest) and instruction == "must_at_least_contain") or \
                        ((any(elem in list(ligand.atomic_props["atoms"]) for elem in atoms_of_interest) == False) and instruction == "must_exclude") or \
                        ((all(elem in atoms_of_interest for elem in list(ligand.atomic_props["atoms"]))) and instruction == "must_only_contain_in_any_amount"):
                    coordinating_atoms_present = True

                if not coordinating_atoms_present:
                    to_delete.append(unq_name)

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}

        self.filter_tracking[len(self.filter_tracking)] = f"Functional Atom filter with {atoms_of_interest}"

    # self.safe_and_document_after_filterstep(filtering_step_name="FunctionalAtoms_of_Interest")

    def filter_betaHs(self, denticities: list = None):
        """
        Filter out all ligands with beta Hydrogen in it
        """
        to_delete = []
        denticities = self.ensure_denticities_is_list(denticities)

        for unq_name, ligand in self.database.db.items():
            if ligand.denticity not in denticities:
                continue

            betaH_present = ligand.betaH_check()
            if betaH_present:
                to_delete.append(unq_name)

        for unq_name in to_delete:
            del self.database.db[unq_name]

        self.filter_tracking[len(self.filter_tracking)] = f"betaH Filter"

    def filter_neighbouring_coordinating_atoms(self, denticities: list = None):
        to_delete = []
        denticities = self.ensure_denticities_is_list(denticities)

        for unq_name, ligand in self.database.db.items():
            if ligand.denticity not in denticities:
                continue

            if ligand.check_for_neighboring_coordinating_atoms():
                to_delete.append(unq_name)

        for unq_name in to_delete:
            del self.database.db[unq_name]

        self.filter_tracking[len(self.filter_tracking)] = f"Neighbouring Atom Filter: {0.2}"

    def filter_min_max_value(self, get_value_from_ligand, min: float = None, max: float = None, denticities: list = None):
        to_delete = []
        denticities = self.ensure_denticities_is_list(denticities)

        # If the user doesn't specify min or max this is set to infinity or -infinity respectively to be ignored
        if min is None or np.isnan(min):
            min = -np.inf
        if max is None or np.isnan(max):
            max = np.inf

        for unq_name, ligand in self.database.db.items():
            if ligand.denticity not in denticities:
                continue

            value = get_value_from_ligand(ligand)
            try:
                if not (min <= value <= max):
                    to_delete.append(unq_name)
            except TypeError:
                if not all(min <= v <= max for v in value):
                    to_delete.append(unq_name)

        for unq_name in to_delete:
            del self.database.db[unq_name]

    def filter_occurrences(self, min: int = None, max: int = None, denticities: list = None):
        self.filter_min_max_value(lambda ligand: ligand.occurrences, min, max, denticities)

    def filter_atom_count(self, min: int = None, max: int = None, denticities: list = None):
        self.filter_min_max_value(lambda ligand: ligand.n_atoms, min, max, denticities)

    def filter_molecular_weight(self, min: float = None, max: float = None, denticities: list = None):
        self.filter_min_max_value(lambda ligand: ligand.global_props['molecular_weight'], min, max, denticities)

    def filter_planarity(self, min: float = None, max: float = None, denticities: list = None):
        self.filter_min_max_value(lambda ligand: ligand.calculate_planarity(), min, max, denticities)

    def filter_interatomic_distances(self, min: float = None, max: float = None, denticities: list = None):
        self.filter_min_max_value(lambda ligand: ligand.get_all_inter_atomic_distances_as_list(), min, max, denticities)



    def filter_charge_confidence(self, filter_for: str):
        """
        Filter out all ligands with a charge assignment that is not confident or confident.
        :param filter_for: "confident" or "not_confident"
        """

        to_delete = []
        if (filter_for is None) or (filter_for != ("confident" or "not_confident")):
            print("!!!Warning!!! -> Arguments specified incorrectly  -> Proceeding to next filter")

        else:
            confident = 0
            not_confident = 0
            for unq_name, ligand in self.database.db.items():
                confidence = ligand.pred_charge_is_confident

                if confidence:
                    confident += 1
                elif not confidence:
                    not_confident += 1
                    to_delete.append(unq_name)
                else:
                    print("!!!Fatal Error!!! -> Charge Confidence Incorrectly Specified  -> Aborting Program")
            # print("Charge Assignment Confident:" + str(confident))
            # print("Charge Assignment Not Confident:" + str(not_confident))
            self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
            self.filter_tracking[len(self.filter_tracking)] = f"Charge Filter: {filter_for}"



    def add_constant_ligands(self):
        """
        Now we add the constant Ligands we defined in constant_Ligands
        """
        print("Adding constant Ligands")

        # for lig in get_monodentate_list() + get_reactant():
        for lig in get_monodentate_list():
            self.database.db[lig.name] = lig

        # self.db.to_json(path=f"{self.safe_path}/DB_after_Adding_Const_Ligands.json")

    def filter_even_odd_electron(self, filter_for: str):
        # filter_for = even --> This will extract all ligands with an even number of electrons
        # filter_for = odd  --> This will extract all ligands with an odd number of electrons
        to_delete = []
        if (filter_for != "even") and (filter_for != "odd"):
            print("!!!Warning!!! -> Arguments specified incorrectly  -> Proceeding to next filter")

        else:
            for unq_name, ligand in self.database.db.items():
                # print(ligand)
                electrons = 0
                for atom in ligand.atomic_props['atoms']:
                    Z = DART_Element(atom).atomic_number
                    electrons += Z
                electrons = electrons * (-1)
                num_electrons = electrons + ligand.pred_charge

                if ((num_electrons % 2 == 0) and (filter_for == "even")) or ((num_electrons % 2 != 0) and (filter_for == "odd")):
                    pass
                    # print("PASS")
                    # RCA_Ligand.view_3d(ligand)
                    # print("")
                else:
                    # print("Fail")
                    # RCA_Ligand.view_3d(ligand)
                    # print("")
                    to_delete.append(unq_name)

            self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
            self.filter_tracking[len(self.filter_tracking)] = f"even_odd_electron_filter: {filter_for}"

    def filter_ligand_charges(self,  charge: Union[list, int], denticity: int = None):
        denticity = self.ensure_denticities_is_list(denticity)

        if not charge is None:
            if not isinstance(charge, (list,tuple)):
                charge = [charge]

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            ligand_charge = ligand.pred_charge
            if ligand.denticity in denticity:
                if ligand_charge not in charge:
                    to_delete.append(unq_name)

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Ligand Charge Filter: [{denticity}] [{charge}]"

    def filter_atomic_neighbors(self, atom: str, neighbors: list, denticities: int = None):
        """
        This function will filter out all ligands in which the specified atom has the specified neighbors in the specified amount. For example, if `atom` is 'C' and `neighbors` is ['H', 'H'], then all ligands in which a 'C' atom has two 'H' neighbors will be filtered out.

        Parameters:
        atom: str -> Chemical element of the central element.
        neighbors: list -> List of chemical elements of the neighbors. Can include repeating elements.
        denticity: int -> The denticity of the ligands to filter. If None, all denticities are considered.
        """
        denticities = self.ensure_denticities_is_list(denticities)

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            if ligand.denticity in denticities:
                if ligand.has_specified_atomic_neighbors(atom, neighbors):
                    to_delete.append(unq_name)

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Atomic Neighbors Filter: {atom} {neighbors}"

    def filter_smarts_substructure_search(self, smarts: str, should_be_present: bool, include_metal: bool = False, denticities: list = None):
        """
        This function will filter out all ligands in which the specified SMARTS pattern is present or not present. If the ligand has no valid SMILES string, it is passed through.
        @param smarts: str -> SMARTS pattern to search for.
        @param should_be_present: bool -> If True, the ligands with the SMARTS pattern present will be kept. If False, the ligands with the SMARTS pattern present will be removed.
        @param include_metal: bool -> If True, the ligands will be filtered based on the SMILES string where the coordinating atoms of the ligand are connected to a metal center.
        @param denticities: list -> The denticity of the ligands to filter. If None, all denticities are considered.
        @return: None
        """
        denticities = self.ensure_denticities_is_list(denticities)

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            if ligand.denticity in denticities:
                with_metal = pseudo_metal if include_metal else None        # Optionally include the metal center in the SMILES string.
                smiles = ligand.get_smiles(with_metal=with_metal)
                if smiles is None:  # If the ligand has no valid SMILES string, pass it through.
                    continue

                has_pattern = has_smarts_pattern(smarts=smarts, smiles=smiles)
                if has_pattern != should_be_present:
                    to_delete.append(unq_name)

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"SMARTS Filter: {smarts} {should_be_present}"

    def filter_symmetric_monodentate_ligands(self, instruction: str = None, threshold: float = None):
        to_delete = []
        if (instruction != "Add") and (instruction != "Remove"):
            print("!!!Warning!!! -> Arguments specified incorrectly  -> Proceeding to next filter")

        else:
            for unq_name, ligand in self.database.db.items():
                if ligand.denticity == 1:
                    x_centroid_list = []
                    x_COM_list = []
                    y_centroid_list = []
                    y_COM_list = []
                    z_centroid_list = []
                    z_COM_list = []
                    all_atomic_masses = []
                    for atom_index in ligand.coordinates:
                        atom_pos = ligand.coordinates[atom_index]
                        atom_type = atom_pos[0]

                        x = atom_pos[1][0]
                        x_centroid_list.append(x)
                        x_COM_list.append(x * DART_Element(atom_type).atomic_mass)

                        y = atom_pos[1][1]
                        y_centroid_list.append(y)
                        y_COM_list.append(y * DART_Element(atom_type).atomic_mass)

                        z = atom_pos[1][2]
                        z_centroid_list.append(z)
                        z_COM_list.append(z * DART_Element(atom_type).atomic_mass)

                        all_atomic_masses.append(DART_Element(atom_type).atomic_mass)

                    # [coord_atom_pos,  centre_of_points, centre_of_mass]
                    coord_atom_pos = np.array(ligand.coordinates[ligand.ligand_to_metal[0]][1])
                    centre_of_points = np.array([sum(x_centroid_list) / len(x_centroid_list), sum(y_centroid_list) / len(y_centroid_list), sum(z_centroid_list) / len(z_centroid_list)])
                    centre_of_mass = np.array([sum(x_COM_list) / sum(all_atomic_masses), sum(y_COM_list) / sum(all_atomic_masses), sum(z_COM_list) / sum(all_atomic_masses)])

                    v1 = coord_atom_pos
                    v2 = centre_of_points
                    v3 = v2 - v1
                    if all(v1 != v2):
                        cosine = np.dot(v1 * (-1), v3) / (np.linalg.norm(v1 * (-1)) * np.linalg.norm(v3))
                        angle = np.arccos(cosine)
                        angle = np.degrees(angle)
                    else:
                        angle = 180

                    if (abs(180 - angle) < threshold) or (abs(0 - angle) < threshold):
                        # print("sucess")
                        # ligand.view_3d()
                        pass
                    else:
                        # print("failure")
                        to_delete.append(unq_name)
                        # ligand.view_3d()

                else:
                    pass
        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"monodentate_filter: {threshold}"

    def ensure_denticities_is_list(self, denticities):
        if isinstance(denticities, int):
            denticities = [denticities]
        elif denticities is None:
            denticities = list(range(-10, 100))

        return denticities

    def filter_missing_bond_orders(self, denticities):
        """
        Filter out all ligands in which not all bond orders are present.
        """
        denticities = self.ensure_denticities_is_list(denticities)

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            if ligand.denticity in denticities:
                if not ligand.has_good_bond_orders:
                    to_delete.append(unq_name)

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Missing Bond Order Filter"

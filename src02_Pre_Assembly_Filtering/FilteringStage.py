import warnings
from pathlib import Path
from copy import deepcopy
import numpy as np
from src01.DataBase import MoleculeDB, LigandDB
from typing import List, Dict, Union
from src02_Pre_Assembly_Filtering.Box_Excluder_Filter import box_filter
from src02_Pre_Assembly_Filtering.constant_Ligands import get_monodentate_list, get_reactant
from pymatgen.core.periodic_table import Element
from src01.Molecule import RCA_Ligand
from rdkit.Chem import rdmolfiles
from src05_Assembly_Refactor.building_block_utility import ligand_to_mol
from rdkit.Chem.Draw import rdMolDraw2D


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

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            # Iterate through each ligand. If the denticity matches
            # the one specified by the user {denticity: int}, but it was never part of
            # complex with a metal specified by the user {metals_of_interest: [str, list[str]]} then it is
            # removed
            if int(ligand.denticity) == int(denticity):
                # print(list(ligand.count_metals))
                metal_of_interest_present = False
                for metal in metals_of_interest:
                    if metal in list(ligand.count_metals):
                        # print("Matching Metal: " + str(metal))
                        metal_of_interest_present = True
                    else:
                        pass
                if metal_of_interest_present:
                    pass
                    # print("Metal of Interest Pass")
                else:
                    to_delete.append(unq_name)
                    # print("Metal of Interest Fail -> Deleting ligand")
            else:
                pass

        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}

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

        to_delete = []
        for identifier, ligand in self.database.db.items():
            if not ligand.denticity in denticity_of_interest:
                to_delete.append(identifier)

        self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items() if identifier not in to_delete}

        # old: self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items() if ligand.denticity in denticity_of_interest}

        self.filter_tracking[len(self.filter_tracking)] = f"Metals of interest filter with {denticity_of_interest}"

        # self.safe_and_document_after_filterstep(filtering_step_name="Denticity_of_Interest")

    def filter_coordinating_group_atoms(self, denticity: int = None, atoms_of_interest: [str, list[str]] = None, instruction: str = None):
        # instruction = must_contain_and_only_contain,this means that if the coordinating atoms specified by the user are exactly the same as the coordinating groups of the ligand then the ligand can pass
        # instruction = must_at_least_contain        ,this means that if the coordinating atoms specified by the user are a subset of that of the ligand then the ligand can pass
        # instruction = must_exclude                 ,this means that if the coordinating atoms specified by the user are not contained in any amount in the ligand then the ligand can pass
        # instruction = must_only_contain_in_any_amount   ,this means that if all coordinating atoms specified by the user are contained to some degree in the ligand and no other coordinating atoms are conatined then the ligand can pass
        # This filter only applies to ligands.py of the specified denticities. ligands.py with other denticities are allowed to pass
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


        to_delete = []
        for unq_name, ligand in self.database.db.items():
            # print("\n")
            # print(unq_name)
            # print(sorted(list(ligand.local_elements)))
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
                    # print("Matching Coordinating groups PASS")
                    pass
                else:
                    # print("Matching Coordinating groups Fail")
                    to_delete.append(unq_name)
            else:
                # If the denticities don't match then we don't apply the filter and just skip
                pass
        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Functional Atom filter with {atoms_of_interest}"

        # self.safe_and_document_after_filterstep(filtering_step_name="FunctionalAtoms_of_Interest")

    def filter_ligand_atoms(self, denticity: int = None, atoms_of_interest: [str, list[str]] = None, instruction: str = None):
        # instruction = must_contain_and_only_contain,this means that if the coordinating atoms specified by the user are exactly the same as the coordinating groups of the ligand then the ligand can pass
        # instruction = must_at_least_contain, this means that if the coordinating atoms specified by the user are a subset of that of the ligand then the ligand can pass
        # instruction = must_exclude, this means that if the coordinating atoms specified by the user are not contained in any amount in the ligand then the ligand can pass
        # instruction = must_only_contain_in_any_amount   ,this means that if all coordinating atoms specified by the user are contained to some degree in the ligand and no other coordinating atoms are conatined then the ligand can pass
        # This filter only applies to ligands.py of the specified denticities. ligands.py with other denticities are allowed to pass
        """
        Only leave in ligands.py where we have functional atoms equal to specified atoms
        """
        print("FunctionalGroup Filter running")

        print(f"Filtering: ligand_atom_type")
        if (atoms_of_interest is None) or (instruction is None) or (denticity is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")
            return

        elif isinstance(atoms_of_interest, str):
            atoms_of_interest = [atoms_of_interest]

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            # print("\n")
            # print(unq_name)
            # print(sorted(list(ligand.atomic_props["atoms"])))
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
                    # print("Matching Coordinating groups PASS")
                    pass
                else:
                    # print("Matching Coordinating groups Fail")
                    to_delete.append(unq_name)
            else:
                # If the denticities don't match then we don't apply the filter and just skip
                pass
        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}

        self.filter_tracking[len(self.filter_tracking)] = f"Functional Atom filter with {atoms_of_interest}"

    # self.safe_and_document_after_filterstep(filtering_step_name="FunctionalAtoms_of_Interest")

    def filter_betaHs(self):
        """
        Filter out all ligands.py with beta Hydrogen in it
        """

        # print("betaH Filter running")

        self.database.db = {identifier: ligand for identifier, ligand in self.database.db.items()
                            if ligand.betaH_check() is False}

        self.filter_tracking[len(self.filter_tracking)] = f"betaH Filter"

        # self.safe_and_document_after_filterstep(filtering_step_name="betaH Filter")

    def filter_neighbouring_coordinating_atoms(self):
        # The goal of this filter is to remove ligands.py that have coordinating atoms close to each other

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            # print("######")
            # print(ligand.atomic_props["atoms"])
            # print(ligand.denticity)
            break_condition = False
            if ligand.denticity != 1:
                if not break_condition:
                    for index_1 in ligand.ligand_to_metal:
                        for index_2 in ligand.ligand_to_metal:
                            if index_1 == index_2:
                                pass
                            elif (index_1 != index_2) and not break_condition:
                                position_1 = np.array(ligand.coordinates[index_1][1])
                                position_2 = np.array(ligand.coordinates[index_2][1])
                                cov_1 = Element(ligand.coordinates[index_1][0]).atomic_radius
                                cov_2 = Element(ligand.coordinates[index_2][0]).atomic_radius
                                distance = np.linalg.norm(position_1 - position_2)
                                if distance < (cov_1 + cov_2 + 0.2):
                                    break_condition = True
                                    # RCA_Ligand.view_3d(ligand)
                                    # print("FAIL")
                                    to_delete.append(unq_name)
                                    break
                                else:
                                    pass

            if not break_condition:
                pass
                # print("PASS")
        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Neighbouring Atom Filter: {0.2}"
        print("finished")

    def filter_molecular_weight(self, denticity: int = None, atomic_weight_min: float = None, atomic_weight_max: float = None):

        to_delete = []

        # If the user doesn't specify min or max this is set to infinity or -infinity respectively to be ignored
        if atomic_weight_min is None:
            atomic_weight_min = -np.nan
        if atomic_weight_max is None:
            atomic_weight_max = np.nan

        else:
            for unq_name, ligand in self.database.db.items():
                if denticity is not None and ligand.denticity != denticity:
                    continue
                # print(unq_name)
                molecular_range_condition_satisified = False
                MW = ligand.global_props["molecular_weight"]
                # print(MW)
                if not (atomic_weight_min <= MW < atomic_weight_max):
                    to_delete.append(unq_name)
                    # print("Molecular Weight Fail")
                # print("\n")
        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Molecular Weight Filter with MW_MIN_{atomic_weight_min} and MW_MAX_{atomic_weight_max}"

    def filter_denticity_fraction(self, denticity: int = None, fraction: float = None):
        # This filter will filter out ligands.py whose dominating denticity does not account for greater than the proportion
        # specified by the user {fraction: float = None} of the total occurences of the complex
        to_delete = []
        if (denticity is None) or (fraction is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")

        else:
            for unq_name, ligand in self.database.db.items():
                # print(unq_name)
                occurence = float(ligand.occurrences)
                # print(occurence)
                current_highest_occurrence = 0
                # print(ligand.count_denticities)
                for key, value in ligand.count_denticities.items():
                    if value > current_highest_occurrence:
                        current_highest_occurrence = value
                    else:
                        pass
                occurence_fraction = float(current_highest_occurrence / occurence)
                # print(occurence_fraction)
                if occurence_fraction > fraction:
                    pass
                    # print("PASS")
                else:
                    to_delete.append(unq_name)
                    # print("FAIL")
                # print("\n")
            self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
            self.filter_tracking[len(self.filter_tracking)] = f"Denticity Fraction: {fraction}"

    def filter_charge_confidence(self, filter_for: str = None):
        # filter_for = "confident"
        # filter_for = "not_confident"

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
            print("Charge Assignment Confident:" + str(confident))
            print("Charge Assignment Not Confident:" + str(not_confident))
            self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
            self.filter_tracking[len(self.filter_tracking)] = f"Charge Filter: {filter_for}"

    def box_excluder_filter(self):
        """
        Filter out all ligands.py that violate the box Filter
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

        # for lig in get_monodentate_list() + get_reactant():
        for lig in get_monodentate_list():
            pass
            self.database.db[lig.name] = lig
            pass

        # self.db.to_json(path=f"{self.safe_path}/DB_after_Adding_Const_Ligands.json")

    def filter_even_odd_electron(self, filter_for: str = None):
        # filter_for = even --> This will extract all ligands.py with an even number of electrons
        # filter_for = odd  --> This will extract all ligands.py with an odd number of electrons
        to_delete = []
        if (filter_for != "even") and (filter_for != "odd"):
            print("!!!Warning!!! -> Arguments specified incorrectly  -> Proceeding to next filter")

        else:
            for unq_name, ligand in self.database.db.items():
                # print(ligand)
                electrons = 0
                for atom in ligand.atomic_props['atoms']:
                    Z = Element(atom).data["Atomic no"]
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
                    RCA_Ligand.view_3d(ligand)
                    print("")
                    to_delete.append(unq_name)

            self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
            self.filter_tracking[len(self.filter_tracking)] = f"even_odd_electron_filter: {filter_for}"

    def filter_ligand_charges(self, denticity: int = None, charge: Union[list,None,int]=None):
        if not charge is None:
            if not isinstance(charge, (list,tuple)):
                charge = [charge]

        to_delete = []
        for unq_name, ligand in self.database.db.items():
            ligand_charge = ligand.pred_charge
            if ligand.denticity == denticity:
                if ligand_charge not in charge:
                    to_delete.append(unq_name)
                else:
                    pass
            else:
                pass
        self.database.db = {unq_name: ligand for unq_name, ligand in self.database.db.items() if unq_name not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Ligand Charge Filter: [{denticity}] [{charge}]"

    def filter_sub_structure_search(self, denticity: int = None, SMARTS: str = None, instruction: str = None):
        warnings.warn("This Filter is still under active developement please do not use")
        raise NotImplementedError
        # filter_for = even --> This will extract all ligands.py with an even number of electrons
        # filter_for = odd  --> This will extract all ligands.py with an odd number of electrons
        new_db = deepcopy(self.database.db)
        if ((instruction != "must_include") and (instruction != "must_exclude")) or (denticity is None) or (SMARTS is None):
            print("!!!Warning!!! -> Arguments specified incorrectly  -> Proceeding to next filter")

        else:
            mol_sub_struct = rdmolfiles.MolFromSmarts(SMARTS)
            for unq_name, ligand in self.database.db.items():
                ligand_mol_string = ligand_to_mol(ligand)
                ligand_mol = rdmolfiles.MolFromMolBlock(ligand_mol_string, removeHs=False, sanitize=False, strictParsing=False)
                test = list(ligand_mol.GetSubstructMatch(mol_sub_struct))
                test_bonds = []
                if test == []:
                    print("No fragment detected")
                else:
                    print("fragment detected")
                    for bond in mol_sub_struct.GetBonds():
                        aid1 = test[bond.GetBeginAtomIdx()]
                        aid2 = test[bond.GetEndAtomIdx()]
                        test_bonds.append(ligand_mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
                    d = rdMolDraw2D.MolDraw2DCairo(500, 500)
                    rdMolDraw2D.PrepareAndDrawMolecule(d, ligand_mol, highlightAtoms=test, highlightBonds=test_bonds)
                    d.FinishDrawing()
                    p = d.GetDrawingText()
                    with open('tmp.png', 'wb') as f:
                        f.write(p)
                    print("done")

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
                        x_COM_list.append(x * Element(atom_type).atomic_mass)

                        y = atom_pos[1][1]
                        y_centroid_list.append(y)
                        y_COM_list.append(y * Element(atom_type).atomic_mass)

                        z = atom_pos[1][2]
                        z_centroid_list.append(z)
                        z_COM_list.append(z * Element(atom_type).atomic_mass)

                        all_atomic_masses.append(Element(atom_type).atomic_mass)

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
        self.database.db = {uname: ligand for uname, ligand in self.database.db.items() if uname not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"monodentate_filter: {threshold}"

    def filter_atom_count(self, denticity: int = None, number: int = None, instruction: str = None):
        to_delete = []
        if (instruction != "greater_than") and (instruction != "less_than"):
            print("!!!Warning!!! -> Arguments specified incorrectly  -> Proceeding to next filter")

        else:
            for unq_name, ligand in self.database.db.items():
                num_atoms = ligand.global_props['n_atoms']
                if ligand.denticity == denticity:
                    if (num_atoms < number) and (instruction == "greater_than"):
                        to_delete.append(unq_name)
                    elif (num_atoms >= number) and (instruction == "less_than"):
                        to_delete.append(unq_name)
                    else:
                        #ligand.view_3d()
                        pass
                else:
                    pass
        self.database.db = {uname: ligand for uname, ligand in self.database.db.items() if uname not in to_delete}
        self.filter_tracking[len(self.filter_tracking)] = f"Atom Number Filter: [{number}] [{instruction}]"

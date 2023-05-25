import yaml
from yaml import SafeLoader
from pathlib import Path
from copy import deepcopy
import warnings
from typing import Union
import json
from constants.constants import transition_metal_symbols
import numpy as np
from pymatgen.core.periodic_table import Element
import networkx as nx
import networkx.algorithms.isomorphism as iso
from src01.utilities_graph import node_check, smiles2nx

from src01.DataBase import LigandDB
from src02_Pre_Assembly_Filtering.box_filter import box_filter
from src02_Pre_Assembly_Filtering.constant_Ligands import get_monodentate_list, get_reactant


class Filter:

    def __init__(self,
                 datapath_,
                 ligand_db: LigandDB = None,
                 max_number: int = None,
                 filter_yaml_path: str = "Filter_setup.yml"
                 ):
        """

        :param ligand_db: If we want to give a custom Ligand DB, rather than reading in from the default document
        """
        self.datapath = datapath_

        with open(filter_yaml_path, "r") as file:
            self.setup = yaml.load(file, SafeLoader)

        self.check_setup_yml()

        if ligand_db is None:
            self.unique_Ligands = LigandDB.from_json(
                json_=str(Path(datapath_, "tmQM_Ligands_unique.json")),
                type_="Ligand",
                max_number=max_number
            )
        else:
            self.unique_Ligands = ligand_db

        self.filtered_db = deepcopy(self.unique_Ligands)

    def check_setup_yml(self):
        """
        this only checks if the setup.yml is in the correct format
        """

        if not "denticities" in self.setup:
            self.setup["denticities"] = [1, 2, 3, 4, 5]
        else:
            if not set(self.setup["denticities"]).issubset({1, 2, 3, 4, 5}):
                raise ValueError("Selected denticities out of range for the assembly")

        if not "metals" in self.setup:
            self.setup["metals"] = {}
        else:
            if not set(self.setup["metals"]).issubset(set(self.setup["denticities"])):
                raise ValueError("Metal filters for ligand denticities selected, which are filtered out already")

            for metal_list in self.setup["metals"].values():
                if not set(metal_list).issubset(set(transition_metal_symbols)):
                    raise ValueError("Non Transiton Metal for Filtering selected")

        if not "coordinating_group_atoms" in self.setup:
            self.setup["coordinating_group_atoms"] = {}
        else:
            if not set(self.setup["coordinating_group_atoms"]).issubset(set(self.setup["denticities"])):
                raise ValueError("Coordinating Group atom Filter for ligand denticities selected, which are filtered out already")

            for dent, specifications in self.setup["coordinating_group_atoms"].items():
                if specifications["instruction"] not in ["must_contain_and_only_contain", "must_at_least_contain",
                                                         "must_exclude", "must_only_contain_in_any_amount"]:
                    warnings.warn("Wrong instruction selected, no filtering of coordinating group atoms")
                    self.setup["coordinating_group_atoms"][dent]["instruction"] = None

        if not "ligand_atoms" in self.setup:
            self.setup["ligand_atoms"] = {}
        else:
            if not set(self.setup["ligand_atoms"]).issubset(set(self.setup["denticities"])):
                raise ValueError("Ligand Atom Filter for ligand denticities selected, which are filtered out already")

            for dent, specifications in self.setup["ligand_atoms"].items():
                if specifications["instruction"] not in ["must_contain_and_only_contain", "must_at_least_contain",
                                                         "must_exclude", "must_only_contain_in_any_amount"]:
                    warnings.warn("Wrong instruction selected, no filtering of Ligand atoms")
                    self.setup["ligand_atoms"][dent]["instruction"] = None

        if not "beta_hydrogen_filter" in self.setup:
            self.setup["beta_hydrogen_filter"] = "no"

        if not "denticity_fraction_filter" in self.setup:
            self.setup["denticity_fraction_filter"] = []
        else:
            if not set(self.setup["denticity_fraction_filter"]).issubset(set(self.setup["denticities"])):
                raise ValueError("Denticity Fraction Filtering for ligand denticities selected, which are filtered out already")

        if not "molecular_weight" in self.setup:
            self.setup["molecular_weight"] = {}
        elif self.setup["molecular_weight"] is None:
            self.setup["molecular_weight"] = {}
        else:
            if not set(self.setup["molecular_weight"]).issubset(set(self.setup["denticities"])):
                raise ValueError(
                    "Molecular weight Filter for ligand denticities selected, which are filtered out already")

            for dent, specifications in self.setup["molecular_weight"].items():
                if specifications["min"] > specifications["max"]:
                    raise ValueError("Selected molecular weight minimum is bigger than molecular weight maximum")

        if not 'submolecule_filter' in self.setup:
            self.setup['submolecule_filter'] = []
        elif self.setup['submolecule_filter'] is None:
            self.setup['submolecule_filter'] = []

        return

    def filter_charge_confidence(self,
                                 filter_for: str = None
                                 ):
        """
        # filter_for = "confident"
        # filter_for = "not_confident"
        """

        new_db = deepcopy(self.filtered_db.db)
        if (filter_for is None) or (filter_for != ("confident" or "not_confident")):
            print("!!!Warning!!! -> Arguments specified incorrectly  -> Proceeding to next filter")

        else:
            confident = 0
            not_confident = 0
            for unq_name, ligand in self.filtered_db.db.items():
                confidence = ligand.pred_charge_is_confident

                if confidence:
                    confident += 1
                elif not confidence:
                    not_confident += 1
                    del new_db[unq_name]
                else:
                    print("!!!Fatal Error!!! -> Charge Confidence Incorrectly Specified  -> Aborting Program")
            print("Charge Assignment Confident:" + str(confident))
            print("Charge Assignment Not Confident:" + str(not_confident))
            self.filtered_db.db = new_db

    def denticity_of_interest_filter(self,
                                     denticity_of_interest: [int, list[int]] = None
                                     ):
        """

        """
        print("Denticity Filter running ...")

        if denticity_of_interest is None:
            warnings.warn("No denticities of interest selected, no filtering at all")
        elif isinstance(denticity_of_interest, int):
            denticity_of_interest = [denticity_of_interest]

        new_db = deepcopy(self.filtered_db.db)

        for identifier, ligand in self.filtered_db.db.items():
            if not ligand.denticity in denticity_of_interest:
                del new_db[identifier]

        self.filtered_db.db = new_db

    def metals_of_interest_filter(self,
                                  denticity: int = None,
                                  metals_of_interest: [str, list[str]] = None
                                  ):
        """
        This function has been updated by Cian to work with the updated .mol file
        """

        print(f"Filtering: Metal of Interest -> Denticity: {denticity} ...")

        if (metals_of_interest is None) or (denticity is None):
            warnings.warn("!!!Warning!!! -> No metal of interest or denticity selected -> Proceeding to next filter")
            return
        elif isinstance(metals_of_interest, str):
            metals_of_interest = [metals_of_interest]

        new_db = deepcopy(self.filtered_db.db)

        for unq_name, ligand in self.filtered_db.db.items():

            if int(ligand.denticity) == int(denticity):
                # print(list(ligand.count_metals))
                metal_of_interest_present = False
                for metal in metals_of_interest:
                    if metal in list(ligand.count_metals):
                        # print("Matching Metal: " + str(metal))
                        metal_of_interest_present = True
                    else:
                        pass
                if metal_of_interest_present is True:
                    pass
                    # print("Metal of Interest Pass")
                else:
                    del new_db[unq_name]
                    # print("Metal of Interest Fail -> Deleting ligand")
            else:
                pass

        self.filtered_db.db = new_db

    def filter_coordinating_group_atoms(self,
                                        denticity: int = None,
                                        atoms_of_interest: [str, list[str]] = None,
                                        instruction: str = "must_contain_and_only_contain"
                                        ):
        """
        :param instruction:
            must_contain_and_only_contain,this means that if the coordinating atoms specified by the user are exactly the same as the coordinating groups of the ligand then the ligand can pass
            must_at_least_contain        ,this means that if the coordinating atoms specified by the user are a subset of that of the ligand then the ligand can pass
            must_exclude                 ,this means that if the coordinating atoms specified by the user are not contained in any amount in the ligand then the ligand can pass
            must_only_contain_in_any_amount   ,this means that if all coordinating atoms specified by the user are contained to some degree in the ligand and no other coordinating atoms are conatined then the ligand can pass

        """
        print("Functional Group Filter running ...\n"
              f"Filtering: Coordinating_atom_type -> Denticity: {denticity}")

        if (atoms_of_interest is None) or (denticity is None) or (instruction is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")
            return
        if ((denticity != len(atoms_of_interest)) and instruction == True) or (
                (denticity < len(atoms_of_interest)) and instruction == False):
            print(
                "!!!Warning!!! -> The specified denticity is not consistent with specified coordinating atoms -> Proceeding to next filter ")

        elif isinstance(atoms_of_interest, str):
            atoms_of_interest = [atoms_of_interest]

        new_db = deepcopy(self.filtered_db.db)

        for unq_name, ligand in self.filtered_db.db.items():

            if int(ligand.denticity) == int(denticity):
                coordinating_atoms_present = False
                # If the denticity of the ligand matches that specified by the user
                if ((sorted(list(ligand.local_elements)) == sorted(
                        atoms_of_interest)) and instruction == "must_contain_and_only_contain") or \
                        (all(elem in list(ligand.local_elements) for elem in
                             atoms_of_interest) and instruction == "must_at_least_contain") or \
                        ((any(elem in list(ligand.local_elements) for elem in
                              atoms_of_interest) == False) and instruction == "must_exclude") or \
                        ((all(elem in atoms_of_interest for elem in
                              list(ligand.local_elements))) and instruction == "must_only_contain_in_any_amount"):
                    coordinating_atoms_present = True
                else:
                    pass
                if coordinating_atoms_present:
                    pass
                else:
                    del new_db[unq_name]
            else:
                # If the denticities don't match then we don't apply the filter and just skip
                pass

        self.filtered_db.db = new_db

    def filter_ligand_atoms(self,
                            denticity: int = None,
                            atoms_of_interest: [str, list[str]] = None,
                            instruction: str = "must_contain_and_only_contain"
                            ):
        """
        instruction similar to filter_coordinating group atoms
        """

        print("FunctionalGroup Filter running ... ")

        print(f"Filtering: ligand_atom_type")
        if (atoms_of_interest is None) or (instruction is None) or (denticity is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")
            return

        elif isinstance(atoms_of_interest, str):
            atoms_of_interest = [atoms_of_interest]

        new_db = deepcopy(self.filtered_db.db)

        for unq_name, ligand in self.filtered_db.db.items():

            if int(ligand.denticity) == int(denticity):
                coordinating_atoms_present = False
                # If the denticity of the ligand matches that specified by the user
                if ((sorted(list(ligand.atomic_props["atoms"])) == sorted(
                        atoms_of_interest)) and instruction == "must_contain_and_only_contain") or \
                        (all(elem in list(ligand.atomic_props["atoms"]) for elem in
                             atoms_of_interest) and instruction == "must_at_least_contain") or \
                        ((any(elem in list(ligand.atomic_props["atoms"]) for elem in
                              atoms_of_interest) == False) and instruction == "must_exclude") or \
                        ((all(elem in atoms_of_interest for elem in list(
                            ligand.atomic_props["atoms"]))) and instruction == "must_only_contain_in_any_amount"):
                    coordinating_atoms_present = True
                else:
                    pass
                if coordinating_atoms_present:
                    pass
                else:
                    del new_db[unq_name]
            else:
                # If the denticities don't match then we don't apply the filter and just skip
                pass
        self.filtered_db.db = new_db

    def filter_betaHs(self):
        """
        Filter out all ligands with beta Hydrogen in it
        """

        print("betaH Filter running ...")
        new_db = deepcopy(self.filtered_db.db)

        self.filtered_db.db = {identifier: ligand for identifier, ligand in new_db.items()
                               if ligand.betaH_check() is False}

        return

    def filter_molecular_weight(self,
                                denticity: int = None,
                                atomic_weight_min: float = None,
                                atomic_weight_max: float = None
                                ):

        new_db = deepcopy(self.filtered_db.db)

        if (denticity is None) or (atomic_weight_min is None) or (atomic_weight_max is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")

        else:
            for unq_name, ligand in self.filtered_db.db.items():

                if (atomic_weight_min <= ligand.global_props["molecular_weight"] < atomic_weight_max) and (
                        ligand.denticity == denticity):
                    pass
                else:
                    del new_db[unq_name]

        self.filtered_db.db = new_db

    def filter_denticity_fraction(self,
                                  denticity: int = None,
                                  fraction: float = 0.9
                                  ):
        """
        This filter will filter out ligands whose dominating denticity does not account for greater than the proportion
        specified by the user {fraction: float = None} of the total occurences of the complex
        """

        new_db = deepcopy(self.filtered_db.db)

        if (denticity is None) or (fraction is None):
            print("!!!Warning!!! -> All arguments not specified  -> Proceeding to next filter")

        else:
            for unq_name, ligand in self.filtered_db.db.items():

                occurence = float(ligand.occurrences)

                current_highest_occurrence = 0

                for key, value in ligand.count_denticities.items():
                    if value > current_highest_occurrence:
                        current_highest_occurrence = value
                    else:
                        pass
                occurence_fraction = float(current_highest_occurrence / occurence)

                if occurence_fraction > fraction:
                    pass
                else:
                    del new_db[unq_name]
            self.filtered_db.db = new_db

    def box_excluder_filter(self):
        """
        Filter out all ligands that violate the box Filter
        """
        print("Mandatory Box Filter running ...")
        self.filtered_db.db = {identifier: ligand for identifier, ligand in self.filtered_db.db.items() if
                               box_filter(ligand) is True}

    def filter_neighbouring_coordinating_atoms(self,
                                               skin_: float = 0.2
                                               ):
        # The goal of this filter is to remove ligands that have coordinating atoms close to each other
        new_db = deepcopy(self.filtered_db.db)
        for unq_name, ligand in self.filtered_db.db.items():

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
                                if distance < (cov_1+cov_2+skin_):
                                    break_condition = True
                                    #RCA_Ligand.view_3d(ligand)
                                    #print("FAIL")
                                    del new_db[unq_name]
                                    break
                                else:
                                    pass

            if not break_condition:
                pass

        self.filtered_db.db = new_db
        print("finished")

    def subgraph_filter(self):
        """
        We perform a subgraph search,
        where we want to see if the inserted molecule is part of the ligand
        input-format: smiles
        # todo: At the moment the filter is also positive if the ligand is a subgraph of the input molecule
        is okay for now as we shall assume that the input molecules are small

        We further assume that we are only interested in the core skeleton and all hydrogens could
        in the ligand be replaced by bonds, thus explicit_H = false

        If the molecule too small, we will skip this filter, i.e. len(subG.nodes) < 2
        as this is the same as a single atom filter and we recommend using another filter type
        """

        new_db = deepcopy(self.filtered_db.db)

        for smiles in self.setup["submolecule_filter"]:

            subG = smiles2nx(smiles_str=smiles, explicit_H=False)

            if len(subG.nodes) < 2:
                warnings.warn(f"The skeleton of the selected molecule for the subgraph search consists of a"
                              f"single molecule. We recommend using the Ligand Atom type filter")
                continue

            for unq_name, ligand in self.filtered_db.db.items():
                G = ligand.graph

                subG_m = nx.MultiGraph(subG)
                G_m = nx.MultiGraph(G)

                if iso.MultiGraphMatcher(G_m, subG_m, node_match=node_check).subgraph_is_isomorphic() is True:
                    try:
                        del new_db[unq_name]
                    except KeyError:
                        # already deleted
                        pass
                else:
                    pass

        self.filtered_db.db = new_db
        return

    def run_filters(self):

        print(f"Number of Ligands before Filtering: {len(self.filtered_db.db)}")

        #
        # mandatory charge_filter
        self.filter_charge_confidence()
        print(f"Number of Ligands after Charge Confidence Filter: {len(self.filtered_db.db)}")

        #
        #
        self.denticity_of_interest_filter(
            denticity_of_interest=self.setup["denticities"]
        )

        print(f"Number of Ligands after denticity Filter: {len(self.filtered_db.db)}")

        #
        #
        if self.setup["metals"] == {}:
            print(f"No metal of interest filtering, number of ligands stays the same")
        else:
            for dent, metal_list in self.setup["metals"].items():
                self.metals_of_interest_filter(denticity=dent,
                                               metals_of_interest=metal_list
                                               )
            print(f"Number of Ligands after metals of interest Filtering: {len(self.filtered_db.db)}")

        #
        #
        if self.setup["coordinating_group_atoms"] == {}:
            print(f"No coordinating group atoms filtering, number of ligands stays the same")
        else:
            for dent, specifications in self.setup["coordinating_group_atoms"].items():
                self.filter_coordinating_group_atoms(denticity=dent,
                                                     atoms_of_interest=specifications["atoms"],
                                                     instruction=specifications["instruction"]
                                                     )
            print(f"Number of Ligands after Coordination-Atoms Filtering: {len(self.filtered_db.db)}")

        #
        #
        if self.setup["ligand_atoms"] == {}:
            print(f"No ligand atoms Filtering, number of ligands stays the same")
        else:
            for dent, specifications in self.setup["ligand_atoms"].items():
                self.filter_ligand_atoms(denticity=dent,
                                         atoms_of_interest=specifications["atoms"],
                                         instruction=specifications["instruction"]
                                         )
            print(f"Number of Ligands after ligand atoms Filtering: {len(self.filtered_db.db)}")

        #
        #
        if self.setup["beta_hydrogen_filter"] in ["y", "yes"]:
            self.filter_betaHs()
            print(f"Number of Ligands after beta Hydrogen Filter: {len(self.filtered_db.db)}")
        else:
            print(f"No beta Hydrogen Filter, number of ligands stays the same")

        #
        #
        if self.setup["molecular_weight"] == {}:
            print(f"No molecular weight number of ligands stays the same")
        else:
            for dent, specifications in self.setup["molecular_weight"].items():
                self.filter_molecular_weight(denticity=dent,
                                             atomic_weight_min=specifications["atoms"],
                                             atomic_weight_max=specifications["instruction"]
                                             )
            print(f"Number of Ligands after Molecular weight filter: {len(self.filtered_db.db)}")

        #
        #
        if not self.setup["denticity_fraction_filter"]:
            print(f"No Ligand Fraction Filter, number of ligands stays the same")
        else:
            for dent in self.setup["denticity_fraction_filter"]:
                self.filter_denticity_fraction(denticity=dent)
            print(f"Number of Ligands after Denticity fraction filter: {len(self.filtered_db.db)}")

        #
        # now we run the subgraph filter
        print(f"Subgraph filter for {len(self.setup['submolecule_filter'])} Molecules running")
        self.subgraph_filter()
        print(f"Number of Ligands after Subgraph filter: {len(self.filtered_db.db)}")

        #
        # mandatory filters
        self.filter_neighbouring_coordinating_atoms()
        self.box_excluder_filter()
        print("Filtering done")

    def add_constant_ligands(self):
        """
        Now we add the constant Ligands we defined in constant_Ligands
        """
        print("Adding constant Ligands")
        for lig in get_monodentate_list():
            self.filtered_db.db[lig.name] = lig

    def add_reactant(self):
        print("Add reactant to active site")
        self.filtered_db.db[get_reactant().name] = get_reactant()

    def safe(self,
             other_safe_path: Union[str, Path] = None,
             overwriting: bool = True
             ):
        """
        Just a method to carefully store the filtering results
        :param other_safe_path: If we don't want to safe it at the default location
        """

        if other_safe_path is not None:
            filtering_safe_path = other_safe_path
        else:
            filtering_safe_path = self.datapath

        db_json_path = Path(filtering_safe_path, "Filtering", "filteredLigDB.json")
        setup_json_path = Path(filtering_safe_path, "Filtering", "filterSetup.json")
        db_json_path.parent.mkdir(parents=True,
                                  exist_ok=True
                                  )

        if db_json_path.exists() and overwriting is False:
            print("File exists and no overwritng selected, no safe")
            return
        else:
            # If we store the filtered DB, we will always overwrite the setup according to the saved db

            self.filtered_db.to_json(path=db_json_path)
            with open(setup_json_path, "w") as file:
                json.dump(self.setup, file)
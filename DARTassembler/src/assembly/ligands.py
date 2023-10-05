import warnings
from copy import deepcopy
import random
import itertools
from DARTassembler.src.assembly.Assemble import PlacementRotation
from typing import Union

from DARTassembler.src.assembly.Assembly_Input import LigandCombinationError


class ChooseRandomLigands:
    def __init__(self, database, topology, instruction, max_attempts, metal_oxidation_state: int, total_complex_charge: int, ligand_choice, max_num_assembled_complexes: int):

        self.ligand_dict = self.make_database_list_consistent_with_topology_and_similarity_list(database=database, topology=topology, instruction=instruction)
        self.topology = topology
        self.instruction = instruction
        self.max_loop = max_attempts
        self.metal_ox = metal_oxidation_state
        self.total_charge = total_complex_charge
        self.ligand_choice = ligand_choice
        self.max_num_assembled_complexes = max_num_assembled_complexes

        self.ligand_dic_list = self.get_ligand_dic(deepcopy_ligands=False)

        self.stop_assembly = False


    @staticmethod
    def format_similarity_lists(input_list: list = None, instruction_list: list = None):
        # This function makes the similarity of one list look like that of another
        # i.e.   [3, 3, 5, 1, 7, 3, 2]       -->     [1, 2, 3, 3, 4, 5, 5]
        #       [3, 3, 5, 5, 7, 3, 3]       -->     [1, 2, 3, 3, 4, 5, 5]
        master_index_list = []
        for identity_code in set(instruction_list):
            # This function here returns the indexes for the of the ligands which should be identical
            index_list = [i for i, val in enumerate(instruction_list) if val == identity_code]
            master_index_list.append(index_list)
        for index_list in master_index_list:
            for i in range(len(index_list) - 1):
                input_list[index_list[i + 1]] = input_list[index_list[0]]
        return input_list

    def make_database_list_consistent_with_topology_and_similarity_list(self, database: Union[dict, list[dict]], topology: list, instruction: list) -> list[dict]:
        """
        This ugly function makes sure that the database list is consistent with the topology and similarity list. The output is a list of databases with the same length as the topology and so that the databases are the same if the similarity is the same.
        """
        if not isinstance(database, list):
            database = [database]*len(topology)
        elif len(database) < len(topology):
            # Repeating similarities, therefore we also have to repeat the database
            db_idx = 0
            new_database = []
            for top_idx in range(len(topology)):
                if top_idx == 0:
                    new_database.append(database[db_idx])
                    db_idx += 1
                else:
                    same_as_last = instruction[top_idx] == instruction[top_idx-1]
                    if same_as_last:
                        new_database.append(new_database[-1])
                    else:
                        new_database.append(database[db_idx])
                        db_idx += 1
            database = new_database

        # Double check
        for idx in range(len(instruction)):
            if idx != 0:
                same_similarities = instruction[idx] == instruction[idx-1]
                if same_similarities:
                    same_databases = id(database[idx]) == id(database[idx-1])
                    assert same_databases, f"Similarities and databases are not consistent at index {idx}!"
        assert len(database) == len(topology), f"The number of topologies and the number of ligand databases are not consistent: {len(topology)} != {len(database)}"

        return database

    def get_charge_dic(self, deepcopy_ligands: bool = True):
        dic1 = []
        dic2 = []
        assert len(self.topology) == len(self.ligand_dict), f"The number of topologies and the number of ligand databases are not consistent: {len(self.topology)} != {len(self.ligand_dict)}"
        for dent, ligands in zip(self.topology, self.ligand_dict):
            tmp_dic_1 = {}  # dictionary with keys (denticty) and values(list of all charges)
            tmp_dic_2 = {}  # dictionary with keys (denticty) and values(list of all ligands in the same order as their charges in tmp_dic_1)
            for denticity in ligands.keys():
                tmp_charge_list = []
                tmp_ligand_list = []
                for ligand in ligands[denticity]:
                    try:
                        charge = ligand.pred_charge
                    except:
                        charge = ligand.global_props["LCS_pred_charge"]
                    tmp_charge_list.append(charge)
                    tmp_ligand_list.append(ligand)
                tmp_dic_1.update({f"{denticity}": deepcopy(tmp_charge_list)})
                tmp_dic_2.update({f"{denticity}": deepcopy(tmp_ligand_list) if deepcopy_ligands else tmp_ligand_list})
            dic1.append(tmp_dic_1)
            dic2.append(tmp_dic_2)
        return dic1, dic2

    def get_ligand_dic(self, deepcopy_ligands: bool = True):
        dic_1_list, dic_2_list = self.get_charge_dic(deepcopy_ligands=False)
        ligands = []

        assert len(self.topology) == len(dic_1_list) == len(dic_2_list), f"The number of topologies and the number of ligand databases are not consistent: {len(self.topology)} != {len(dic_1_list)} != {len(dic_2_list)}"
        for dent, dic_1, dic_2 in zip(self.topology, dic_1_list, dic_2_list):
            tmp_dic_3 = {}  # tmp_dic_3 is a dictionary with keys (denticity) and  value (dictionary). This dictionary has keys (unique charge) and values(ligand building blocks))
            for denticity, charge_list, ligand_list in zip(dic_1.keys(), dic_1.values(), dic_2.values()):
                tmp_dic_charge = {}
                for unq_charge in set(charge_list):
                    tmp_list = []
                    for charge, ligand in zip(charge_list, ligand_list):
                        if str(unq_charge) == str(charge):
                            tmp_list.append(ligand)
                        else:
                            pass
                    tmp_dic_charge.update({f"{unq_charge}": tmp_list})
                tmp_dic_3.update({f"{denticity}": deepcopy(tmp_dic_charge) if deepcopy_ligands else tmp_dic_charge})
            ligands.append(tmp_dic_3)
        return ligands

    def charge_list_process(self):
        print("\nStarting Charge Loop")
        m = 0
        charge_dic, _ = self.get_charge_dic(deepcopy_ligands=False)     # No deepcopy for speedup since we don't need the ligands
        while m < self.max_loop:
            charge_list = []
            for dent, charges in zip(self.topology, charge_dic):
                charge_list.append(random.choice(charges[str(dent)]))

            charge_list_out = self.format_similarity_lists(charge_list, self.instruction)
            if sum(charge_list_out) == self.total_charge - self.metal_ox:
                print(f"Charge Resolved After [{m}] Iterations\n")
                return charge_list_out
            else:
                pass
            m += 1
        warnings.warn(
            f"!!!Fatal Error!!! -> The total charge condition [{self.total_charge}] and metal oxidation state [{self.metal_ox}] assigned to the complex [{self.topology} -- {self.instruction}] is not solvable in a realistic time frame -> Exiting Program")
        return None

    def choose_ligands(self) -> Union[dict,None]:
        if self.ligand_choice == "random":
            return self.choose_ligands_randomly()
        elif self.ligand_choice == "all":
            return self.choose_ligands_iteratively()
        else:
            raise ValueError(f"Unknown ligand choice method: {self.ligand_choice}")

    def if_make_more_complexes(self, num_assembled_complexes: int) -> bool:
        max_complexes_reached = num_assembled_complexes >= self.max_num_assembled_complexes
        if self.ligand_choice == "random":
            return not max_complexes_reached
        elif self.ligand_choice == "all":
            return not self.stop_assembly and not max_complexes_reached
        else:
            raise ValueError(f"Unknown ligand choice method: {self.ligand_choice}")

    def choose_ligands_randomly(self) -> Union[dict,None]:
        while True:
            charge_list = self.charge_list_process()
            if charge_list is None:
                raise LigandCombinationError(f'Could not find a valid combination of ligands! Most likely there was no valid combination of ligand charges which fulfills the total charge requirement for the specified metal oxidation state MOS={self.metal_ox_state} and total charge Q={self.total_charge}. Please check your ligand database and/or your assembly input file.')

            # Choose ligands randomly
            ligands = {}
            for i, (denticity, charge, ligand_dic) in enumerate(zip(self.topology, charge_list, self.ligand_dic_list)):
                ligands.update({i: random.choice(ligand_dic[str(denticity)][str(charge)])})

            ligands_out = self.format_similarity_lists(ligands, self.instruction) # replace ligands with the same similarity with the same ligand

            # Deepcopy all chosen ligands. Only doing this here at the end instead of at each intermediate step makes for a huge speedup
            ligands_out = {idx: deepcopy(ligand) for idx, ligand in ligands_out.items()}

            yield ligands_out

    def choose_ligands_iteratively(self) -> dict:
        concat_ligand_dic_list = []
        for dent, ligand_dic in zip(self.topology, self.ligand_dic_list):
            dent = str(dent)
            all_ligands = []
            ligand_dic = ligand_dic[dent]
            for charge in ligand_dic.keys():
                all_ligands.extend(ligand_dic[charge])
            concat_ligand_dic_list.append(all_ligands)
        all_combs = list(itertools.product(*concat_ligand_dic_list))

        i = -1
        chosen_ligand_combinations = set()
        n_combinations = len(all_combs)
        while i < n_combinations:
            i += 1      # infinite for loop
            try:
                ligands = all_combs[i]
            except IndexError:
                self.stop_assembly = True
                yield None

            assert [lig.denticity for lig in ligands] == self.topology, f"Topology of ligands {ligands} does not match the desired topology {self.topology}! This should not happen!"

            ligands = {idx: lig for idx, lig in enumerate(ligands)}
            ligands = self.format_similarity_lists(ligands, self.instruction)  # replace ligands with the same similarity with the same ligand
            ligand_names = tuple(sorted([ligand.unique_name for ligand in ligands.values()]))

            if ligand_names in chosen_ligand_combinations:
                # This combination has already been chosen
                continue
            else:
                chosen_ligand_combinations.add(ligand_names)

            # Check total charge condition
            sum_ligand_charges = sum([ligand.pred_charge for ligand in ligands.values()])
            if sum_ligand_charges == self.total_charge - self.metal_ox:
                # Total charge condition is satisfied
                ligands = {idx: deepcopy(ligand) for idx, ligand in ligands.items()}
                yield ligands

# This class is deprecated
# class ChooseIterativeLigands:
#
#     def __init__(self, database, top_list, charge, metal, random_seed):
#         self.top_list = top_list
#         self.topology = None
#         self.similarity_list = None
#         self.ligand_dict = database
#         self.metal_list = metal
#         self.Total_Charge = charge
#         self.Random_Seed = random_seed
#         self.charge_dic, self.ligand_dic = self.get_charge_dic()
#         self.perform_pre_checks()
#         self.Combinations = None
#
#     def get_charge_dic(self):
#         tmp_dic_1 = {}  # dictionary with keys (denticity) and values(list of all charges)
#         tmp_dic_2 = {}  # dictionary with keys (denticity) and values(list of all ligands in the same order as their charges in tmp_dic_1)
#         for denticity in self.ligand_dict.keys():
#             tmp_charge_list = []
#             tmp_ligand_list = []
#             for ligand in self.ligand_dict[denticity]:
#                 try:
#                     charge = ligand.pred_charge
#                 except:
#                     charge = ligand.global_props["LCS_pred_charge"]
#                 tmp_charge_list.append(charge)
#                 tmp_ligand_list.append(ligand)
#             tmp_dic_1.update({f"{denticity}": deepcopy(tmp_charge_list)})
#             tmp_dic_2.update({f"{denticity}": deepcopy(tmp_ligand_list)})
#         return tmp_dic_1, tmp_dic_2
#
#     def perform_pre_checks(self):
#         if len(self.top_list) != 1:
#             # Ensures only one topology provided
#             print("!!!Fatal Error!!! -> In order to Assemble complexes deterministically only one topology can be provided -> Exiting Program")
#             raise ValueError
#
#         elif len(self.metal_list) != 1:
#             # Ensures only one metal provided
#             print("!!!Fatal Error!!! -> In order to Assemble complexes deterministically only one metal can be provided -> Exiting Program")
#             raise ValueError
#
#         self.topology, self.similarity_list = PlacementRotation.format_topologies(self.top_list[0])
#         if sorted(set(self.topology)) != sorted(set([int(x) for x in self.charge_dic.keys()])):
#             print("!!!Fatal Error!!! -> The ligands provided in the database are not consistent with the desired denticities in the input-> Exiting Program")
#             raise ValueError
#
#         for denticity, charge_list in self.charge_dic.items():
#             if len(set(charge_list)) != 1:
#                 print(f"!!!Fatal Error!!! -> [{len(set(charge_list))}] different charges detected for denticity [{denticity}], This value should be 1 -> Exiting Program")
#                 raise ValueError
#             else:
#                 pass
#
#         if self.Random_Seed != 0:
#             print(f"!!!Fatal Error!!! -> The random seed must be set to zero not [{self.Random_Seed}]-> Exiting Program")
#             raise ValueError
#
#         print("!!!Success!!! -> ChooseIterativeLigands pre-checks have been PASSED -> Progressing to Assembly")
#
#     def obtain_all_combinations(self):
#         tmp_list = []
#         for denticity in self.topology:
#             tmp_list.append(self.ligand_dict[denticity])
#         combo_list = list(itertools.product(*tmp_list))
#         self.Combinations = [list(ele) for ele in combo_list]
#         self.eliminate_duplicates()
#         ans = input(f"DART has determined that there are [{len(self.Combinations)}] possible complexes that could be assembled\n"
#                     f"Would you like to proceed (y/n)")
#         if (ans == "y") or (ans == "yes") or (ans == "Yes") or (ans == "Y") or (ans == "YES"):
#             pass
#         else:
#             exit()
#
#     def eliminate_duplicates(self):
#         # The purpose of this function is to eliminate duplicates from the self.Combinations list
#         # i.e. lists of ligands that are the same apart from order
#         # It is the job of the isomer handler to build isomers not the job of the ligand decider
#         tmp_list = []
#         for ligand_combo_1 in self.Combinations:
#             detected = False
#             for ligand_combo_2 in self.Combinations:
#                 if set(ligand_combo_1) == set(ligand_combo_2) and ((self.Combinations.index(ligand_combo_1)) != (self.Combinations.index(ligand_combo_2))):
#                     print("duplicate detected")
#                     detected = True
#                 else:
#                     pass
#             if not detected:
#                 tmp_list.append(ligand_combo_1)
#             else:
#                 pass
#         self.Combinations = tmp_list
#
#     def choose_ligands(self, iteration: int = None):
#         tmp_dic = {}
#         i = 0
#         for ligand in self.Combinations[iteration]:
#             tmp_dic.update({i: ligand})
#             i = i + 1
#         return tmp_dic


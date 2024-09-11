import warnings
from copy import deepcopy
import random
import itertools
from typing import Union
import logging
from DARTassembler.src.assembly.Assembly_Input import LigandCombinationError

class LigandChoice(object):

    def __init__(self, database, topology, instruction, metal_oxidation_state: int, total_complex_charge: int, max_num_assembled_complexes: Union[int,str]):
        """
        This class is used to choose ligands for the assembly of complexes. It supports both random and iterative ligand choice methods.
        """
        self.ligand_lists = self._get_relevant_ligand_db(database=database, topology=topology, instruction=instruction)
        self.topology = topology
        self.instruction = instruction
        self.metal_ox = metal_oxidation_state
        self.total_charge = total_complex_charge
        self.max_num_assembled_complexes = max_num_assembled_complexes  # int or "all"
        self.ligand_choice = 'all' if max_num_assembled_complexes == "all" else 'random'

        self.continue_assembly = True   # If set to False, the assembler will stop
        self.max_rejected_ligand_combinations = 1_000  # If this many ligand combinations are rejected in a row, the random choice of ligands is exhausted and will be switched to iterative mode
        self.switched_to_iterative = False  # Flag to output a warning if the ligand choice method was switched to "all"

    def _check_good_ligand_charges(self, ligand_combination: list):
        """
        Check if these ligands have the correct sum of charges.
        """
        sum_ligand_charges = sum([ligand.pred_charge for ligand in ligand_combination])
        correct_charges = sum_ligand_charges == self.total_charge - self.metal_ox
        return correct_charges

    def if_make_more_complexes(self, num_assembled_complexes: int) -> bool:
        """
        Returns True if more complexes can be assembled, False otherwise.
        """
        if self.max_num_assembled_complexes == "all":
            return self.continue_assembly
        else: # self.max_num_assembled_complexes is an integer
            return num_assembled_complexes < self.max_num_assembled_complexes

    def _choose_random_ligand_combination_from_db(self) -> list:
        """
        Choose ligands randomly from the ligand databases.
        """
        ligand_combination = []
        for ligand_list in self.ligand_lists:
            # Choose ligands randomly and respect the "same_as_previous" instruction
            chosen_ligand = random.choice(ligand_list) if ligand_list != 'same_as_previous' else ligand_combination[-1]
            ligand_combination.append(chosen_ligand)

        return ligand_combination

    def _choose_iterative_ligand_combination_from_db(self, all_combinations) -> list:
        """
        Choose ligands iteratively from the ligand databases.
        """
        # Choose ligands iteratively, but if the last entry in the topology is "same_as_previous", all_combinations only includes the ligand lists before that. Therefore we have to add it later.
        try:
            prel_ligand_combination = next(all_combinations) # preliminary ligand combination without respecting "same_as_previous"
        except StopIteration:   # No more valid ligand combinations
            return None

        # Add the last ligand to the list of ligands if the last entry in the topology is "same_as_previous"
        ligand_combination = []
        for idx, ligand_list in enumerate(self.ligand_lists):
            if ligand_list == 'same_as_previous':
                chosen_ligand = prel_ligand_combination[-1]
            else:
                chosen_ligand = prel_ligand_combination[idx]
            ligand_combination.append(chosen_ligand)

        return ligand_combination

    def _final_assertions_for_ligand_combination(self, ligand_combination: list) -> None:
        """
        Doublechecks the ligand combination for consistency with all constraints.
        """
        # Charges
        sum_of_charges = sum([ligand.pred_charge for ligand in ligand_combination])
        assert sum_of_charges == self.total_charge - self.metal_ox, f"The sum of charges of the ligand combination {ligand_combination} is not equal to the total charge {self.total_charge} - the metal oxidation state {self.metal_ox} = {self.total_charge - self.metal_ox}! This should not happen!"

        for idx, ligand in enumerate(ligand_combination):
            # Correct denticities
            assert ligand.denticity == self.topology[idx], f"The denticity of ligand {ligand} at index {idx} is not equal to the desired denticity {self.topology[idx]}! This should not happen!"

            # Correct same_as_previous
            if self.ligand_lists[idx] == 'same_as_previous':
                assert ligand.unique_name == ligand_combination[idx-1].unique_name, f"The ligand {ligand.unique_name} at index {idx} is not the same as the ligand {ligand_combination[idx-1].unique_name} at index {idx-1}! This should not happen!"
            if idx > 0 and (self.instruction[idx] == self.instruction[idx-1]):
                assert ligand.unique_name == ligand_combination[idx-1].unique_name, f"The ligand {ligand.unique_name} at index {idx} is not the same as the ligand {ligand_combination[idx-1].unique_name} at index {idx-1}! This should not happen!"

        # Another doublecheck for the similarity lists.
        # Build the similarity list from the ligand combination and check if it is consistent with the input similarity list.
        has_similarity = []
        counter = 1
        for idx, ligand in enumerate(ligand_combination):
            if idx == 0:    # First ligand always has similarity 1
                has_similarity.append(counter)
                counter += 1
            else:
                if self.ligand_lists[idx] == 'same_as_previous':
                    has_similarity.append(has_similarity[-1])
                else:
                    has_similarity.append(counter)
                    counter += 1
        assert has_similarity == self.instruction, f"The similarity lists are not consistent with the ligand combination {ligand_combination} and the topology {self.topology}! This should not happen!"

        return

    def choose_ligands(self) -> Union[dict,None]:
        """
        Choose ligands for the assembly of complexes. This function is a generator and yields ligand combinations. There are two modes, random and iterative, which are chosen by the ligand_choice attribute.
        - ligand_choice = 'all': If the mode is iterative, it's very simple: All structures will be made, though the order is non-random.
        - ligand_choice = 'random': If the mode is random, the function will yield a random ligand combination each time it is called. The function will stop yielding ligand combinations if the maximum number of complexes has been reached or if no more valid ligand combinations can be found. In the latter case, the random mode will also switch to iterative mode to make sure all possible complexes are made.
        """

        # Setup all ligand combinations as iterable. Needed for the iterative ligand choice method.
        assert self.ligand_lists[-1] == 'same_as_previous' if 'same_as_previous' in self.ligand_lists else True, "The 'same_as_previous' instruction must always come last in the list of ligand lists!" # HARDCODED: If the 'same_as_previous' instruction is used, it always comes last in the list of ligand lists
        all_valid_lists = [ligands for ligands in self.ligand_lists if ligands != 'same_as_previous']
        all_ligand_combinations = itertools.product(*all_valid_lists)

        chosen_ligand_combinations = set()  # Store all chosen ligand combinations to avoid duplicates
        count_rejected_ligand_combinations_in_a_row = 0  # Count how many times the same ligand combination has been chosen in a row. If this number gets too high, the random choice of ligands is probably exhausted and the ligand choice method will be switched to "all".
        while True:     # infinite loop, will be broken by function if_make_more_complexes().

            # Choose ligands for a complex
            if self.ligand_choice == 'random':
                ligand_combination = self._choose_random_ligand_combination_from_db()
            elif self.ligand_choice == 'all':
                ligand_combination = self._choose_iterative_ligand_combination_from_db(all_ligand_combinations)
                if ligand_combination is None: # No more valid ligand combinations
                    self.continue_assembly = False
                    break
            else:
                raise ValueError(f"Unknown ligand choice method: {self.ligand_choice}")

            # Check if ligand charges sum up to correct total charge
            if not self._check_good_ligand_charges(ligand_combination):
                continue

            # If a lot of ligand combinations are rejected in a row, this might indicate that the random choice of ligands is exhausted. In this case, switch to iterative mode to make sure all possible complexes are made.
            if self.ligand_choice == 'random' and count_rejected_ligand_combinations_in_a_row > self.max_rejected_ligand_combinations:
                logging.debug(f"DART warning: Random choice of ligands seems exhausted. Switch to iterative ligand choice.")
                self.ligand_choice = 'all'
                self.switched_to_iterative = True   # set flag for printing an info statement later

            # Check if this combination has already been chosen before. Also important if ligand_choice == "all", because the same combination can be chosen multiple times at different ligand sites because of how itertools.product works.
            ligand_names = tuple(sorted(ligand.unique_name for ligand in ligand_combination))   # Sort ligand names to make sure the same combination is always represented by the same tuple
            if ligand_names in chosen_ligand_combinations:  # This combination has already been chosen
                count_rejected_ligand_combinations_in_a_row += 1
                continue
            else:
                count_rejected_ligand_combinations_in_a_row = 0
                chosen_ligand_combinations.add(ligand_names)

            # Final check if the ligand combination fulfills all constraints
            self._final_assertions_for_ligand_combination(ligand_combination)

            # Deepcopy all chosen ligands for safety. Only doing this here at the end instead of at each intermediate step makes for a huge speedup.
            # Note: Huge bottleneck. Leaving out the deepcopy() speeds up the code by a factor of 3. Potentially dangerous though, so if we fix it we should make sure in the code that it is not needed. Preliminary tests with small numbers of complexes show no difference in the results though.
            # ligands_out = {idx: deepcopy(ligand) for idx, ligand in enumerate(ligand_combination)}
            ligands_out = {idx: ligand for idx, ligand in enumerate(ligand_combination)}    # todo

            yield ligands_out

        if len(chosen_ligand_combinations) == 0: # Output error because no valid ligand combinations found
            raise LigandCombinationError(
                f'No valid ligand combinations found which fulfill the metal oxidation state MOS={self.metal_ox} and total charge Q={self.total_charge} requirement! This can happen when the provided metal oxidation state or total charge are too high/low. Please check your ligand database and/or your assembly input file.')

        if self.switched_to_iterative:
            logging.info(f'DART info: This batch was interrupted early because all possible complexes have already been made.')

        self.continue_assembly = False # Stop assembly after the maximum number of complexes has been reached

    def _get_only_relevant_denticities_in_db(self, database: list, topology: list) -> list:
        """
        This function takes a database and a topology and returns a database with only the relevant denticities.
        """
        relevant_db = []
        for idx, dent in enumerate(topology):
            entry = database[idx]
            if entry != 'same_as_previous':
                try:
                    correct_dent_ligands = database[idx][dent]
                except KeyError:
                    raise LigandCombinationError(
                        f'The provided ligand database doesn\'t contain the denticity {dent} in the topology {self.topology} at the {idx + 1}th site! Please check your ligand database and/or your assembly input file')
                relevant_db.append(correct_dent_ligands)
            else:
                relevant_db.append('same_as_previous')

        return relevant_db

    def _get_relevant_ligand_db(self, database: Union[dict, list[dict]],
                                topology: list, instruction: list) -> list[dict]:
        """
        This ugly function makes sure that the database list is consistent with the topology and similarity list. The output is a list of databases with the same length as the topology and so that the databases are the same if the similarity is the same.
        """
        if not isinstance(database, list):
            database = [database] * len(topology)
        elif len(database) < len(topology):
            # Repeating similarities, therefore we also have to repeat the database
            db_idx = 0
            new_database = []
            for top_idx in range(len(topology)):
                if top_idx == 0:
                    new_database.append(database[db_idx])
                    db_idx += 1
                else:
                    same_as_last = instruction[top_idx] == instruction[top_idx - 1]
                    if same_as_last:
                        new_database.append('same_as_previous')
                    else:
                        new_database.append(database[db_idx])
                        db_idx += 1
            database = new_database

        # Double check
        assert len(database) == len(
            topology), f"The number of topologies and the number of ligand databases are not consistent: {len(topology)} != {len(database)}"

        database = self._get_only_relevant_denticities_in_db(database=database, topology=topology)
        return database



# ===========   OUTDATED CODE FOR LIGAND CHOICE   ==========
#
# class ChooseRandomLigands:
#     def __init__(self, database, topology, instruction, max_attempts, metal_oxidation_state: int, total_complex_charge: int, ligand_choice, max_num_assembled_complexes: int):
#
#         self.ligand_dict = self.make_database_list_consistent_with_topology_and_similarity_list(database=database, topology=topology, instruction=instruction)
#         self.topology = topology
#         self.instruction = instruction
#         self.max_loop = max_attempts
#         self.metal_ox = metal_oxidation_state
#         self.total_charge = total_complex_charge
#         self.ligand_choice = ligand_choice
#         self.max_num_assembled_complexes = max_num_assembled_complexes
#
#         self.ligand_dic_list = self.get_ligand_dic(deepcopy_ligands=False)
#
#         self.stop_assembly = False
#
#
#     @staticmethod
#     def format_similarity_lists(input_list: list = None, instruction_list: list = None):
#         # This function makes the similarity of one list look like that of another
#         # i.e.   [3, 3, 5, 1, 7, 3, 2]       -->     [1, 2, 3, 3, 4, 5, 5]
#         #       [3, 3, 5, 5, 7, 3, 3]       -->     [1, 2, 3, 3, 4, 5, 5]
#         master_index_list = []
#         for identity_code in set(instruction_list):
#             # This function here returns the indexes for the of the ligands which should be identical
#             index_list = [i for i, val in enumerate(instruction_list) if val == identity_code]
#             master_index_list.append(index_list)
#         for index_list in master_index_list:
#             for i in range(len(index_list) - 1):
#                 input_list[index_list[i + 1]] = input_list[index_list[0]]
#         return input_list
#
#     def make_database_list_consistent_with_topology_and_similarity_list(self, database: Union[dict, list[dict]], topology: list, instruction: list) -> list[dict]:
#         """
#         This ugly function makes sure that the database list is consistent with the topology and similarity list. The output is a list of databases with the same length as the topology and so that the databases are the same if the similarity is the same.
#         """
#         if not isinstance(database, list):
#             database = [database]*len(topology)
#         elif len(database) < len(topology):
#             # Repeating similarities, therefore we also have to repeat the database
#             db_idx = 0
#             new_database = []
#             for top_idx in range(len(topology)):
#                 if top_idx == 0:
#                     new_database.append(database[db_idx])
#                     db_idx += 1
#                 else:
#                     same_as_last = instruction[top_idx] == instruction[top_idx-1]
#                     if same_as_last:
#                         new_database.append(new_database[-1])
#                     else:
#                         new_database.append(database[db_idx])
#                         db_idx += 1
#             database = new_database
#
#         # Double check
#         for idx in range(len(instruction)):
#             if idx != 0:
#                 same_similarities = instruction[idx] == instruction[idx-1]
#                 if same_similarities:
#                     same_databases = id(database[idx]) == id(database[idx-1])
#                     assert same_databases, f"Similarities and databases are not consistent at index {idx}!"
#         assert len(database) == len(topology), f"The number of topologies and the number of ligand databases are not consistent: {len(topology)} != {len(database)}"
#
#         return database
#
#     def get_charge_dic(self, deepcopy_ligands: bool = True):
#         dic1 = []
#         dic2 = []
#         assert len(self.topology) == len(self.ligand_dict), f"The number of topologies and the number of ligand databases are not consistent: {len(self.topology)} != {len(self.ligand_dict)}"
#         for dent, ligands in zip(self.topology, self.ligand_dict):
#             tmp_dic_1 = {}  # dictionary with keys (denticty) and values(list of all charges)
#             tmp_dic_2 = {}  # dictionary with keys (denticty) and values(list of all ligands in the same order as their charges in tmp_dic_1)
#             for denticity in ligands.keys():
#                 tmp_charge_list = []
#                 tmp_ligand_list = []
#                 for ligand in ligands[denticity]:
#                     try:
#                         charge = ligand.pred_charge
#                     except:
#                         charge = ligand.global_props["LCS_pred_charge"]
#                     tmp_charge_list.append(charge)
#                     tmp_ligand_list.append(ligand)
#                 tmp_dic_1.update({f"{denticity}": deepcopy(tmp_charge_list)})
#                 tmp_dic_2.update({f"{denticity}": deepcopy(tmp_ligand_list) if deepcopy_ligands else tmp_ligand_list})
#             dic1.append(tmp_dic_1)
#             dic2.append(tmp_dic_2)
#         return dic1, dic2
#
#     def get_ligand_dic(self, deepcopy_ligands: bool = True):
#         dic_1_list, dic_2_list = self.get_charge_dic(deepcopy_ligands=False)
#         ligands = []
#
#         assert len(self.topology) == len(dic_1_list) == len(dic_2_list), f"The number of topologies and the number of ligand databases are not consistent: {len(self.topology)} != {len(dic_1_list)} != {len(dic_2_list)}"
#         for dent, dic_1, dic_2 in zip(self.topology, dic_1_list, dic_2_list):
#             tmp_dic_3 = {}  # tmp_dic_3 is a dictionary with keys (denticity) and  value (dictionary). This dictionary has keys (unique charge) and values(ligand building blocks))
#             for denticity, charge_list, ligand_list in zip(dic_1.keys(), dic_1.values(), dic_2.values()):
#                 tmp_dic_charge = {}
#                 for unq_charge in set(charge_list):
#                     tmp_list = []
#                     for charge, ligand in zip(charge_list, ligand_list):
#                         if str(unq_charge) == str(charge):
#                             tmp_list.append(ligand)
#                         else:
#                             pass
#                     tmp_dic_charge.update({f"{unq_charge}": tmp_list})
#                 tmp_dic_3.update({f"{denticity}": deepcopy(tmp_dic_charge) if deepcopy_ligands else tmp_dic_charge})
#             ligands.append(tmp_dic_3)
#         return ligands
#
#     def charge_list_process(self):
#         logging.debug("\nStarting Charge Loop")
#         m = 0
#         charge_dic, _ = self.get_charge_dic(deepcopy_ligands=False)     # No deepcopy for speedup since we don't need the ligands
#         while m < self.max_loop:
#             charge_list = []
#             for dent, charges in zip(self.topology, charge_dic):
#                 charge_list.append(random.choice(charges[str(dent)]))
#
#             charge_list_out = self.format_similarity_lists(charge_list, self.instruction)
#             if sum(charge_list_out) == self.total_charge - self.metal_ox:
#                 logging.debug(f"Charge Resolved After [{m}] Iterations\n")
#                 return charge_list_out
#             else:
#                 pass
#             m += 1
#         logging.warning(
#             f"!!!Fatal Error!!! -> The total charge condition [{self.total_charge}] and metal oxidation state [{self.metal_ox}] assigned to the complex [{self.topology} -- {self.instruction}] is not solvable in a realistic time frame -> Exiting Program")
#         return None
#
#     def choose_ligands(self) -> Union[dict,None]:
#         if self.ligand_choice == "random":
#             return self.choose_ligands_randomly()
#         elif self.ligand_choice == "all":
#             return self.choose_ligands_iteratively()
#         else:
#             raise ValueError(f"Unknown ligand choice method: {self.ligand_choice}")
#
#     def if_make_more_complexes(self, num_assembled_complexes: int) -> bool:
#         max_complexes_reached = num_assembled_complexes >= self.max_num_assembled_complexes
#         if self.ligand_choice == "random":
#             return not max_complexes_reached
#         elif self.ligand_choice == "all":
#             return not self.stop_assembly and not max_complexes_reached
#         else:
#             raise ValueError(f"Unknown ligand choice method: {self.ligand_choice}")
#
#     def choose_ligands_randomly(self) -> Union[dict,None]:
#         chosen_ligand_combinations = set()
#         # for _ in range(self.max_num_assembled_complexes):
#         # while len(chosen_ligand_combinations) < self.max_num_assembled_complexes:
#         while True:
#             # state_before_choice = random.getstate()  # Save the state of the random number generator for debugging. todo: remove
#             charge_list = self.charge_list_process()
#             # random.setstate(state_before_choice)  # Reset the state of the random number generator for debugging. todo: remove
#             if charge_list is None:
#                 raise LigandCombinationError(f'No valid ligand combinations found which fulfills the metal oxidation state MOS={self.metal_ox} and total charge Q={self.total_charge} requirement! This can happen when the provided metal oxidation state or total charge are too high/low. Please check your ligand database and/or your assembly input file.')
#
#             # Choose ligands randomly
#             ligands = {}
#             for i, (denticity, charge, ligand_dic) in enumerate(zip(self.topology, charge_list, self.ligand_dic_list)):
#                 ligands.update({i: random.choice(ligand_dic[str(denticity)][str(charge)])})
#
#             ligands_out = self.format_similarity_lists(ligands, self.instruction) # replace ligands with the same similarity with the same ligand
#
#             # Check if this combination has already been chosen before
#             ligand_names = tuple(sorted(ligand.unique_name for ligand in ligands_out.values()))
#             if ligand_names in chosen_ligand_combinations:
#                 # This combination has already been chosen
#                 continue
#             else:
#                 chosen_ligand_combinations.add(ligand_names)
#
#             # Deepcopy all chosen ligands for safety. Only doing this here at the end instead of at each intermediate step makes for a huge speedup
#             ligands_out = {idx: deepcopy(ligand) for idx, ligand in ligands_out.items()}
#
#             yield ligands_out
#
#     def choose_ligands_iteratively(self) -> dict:
#         concat_ligand_dic_list = []
#         for idx, (dent, ligand_dic) in enumerate(zip(self.topology, self.ligand_dic_list)):
#             dent = str(dent)
#             all_ligands = []
#             try:
#                 ligand_dic = ligand_dic[dent]
#                 for charge in ligand_dic.keys():
#                     all_ligands.extend(ligand_dic[charge])
#             except KeyError:
#                 raise LigandCombinationError(f'The provided ligand database doesn\'t contain the denticity {dent} in the topology {self.topology} at the {idx+1}th site! Please check your ligand database and/or your assembly input file')
#             concat_ligand_dic_list.append(all_ligands)
#         all_combs = list(itertools.product(*concat_ligand_dic_list))
#
#         if len(all_combs) == 0:
#             raise LigandCombinationError(f'No valid ligand combinations found which fulfills the metal oxidation state MOS={self.metal_ox} and total charge Q={self.total_charge} requirement! This can happen when the provided metal oxidation state or total charge are too high/low. Please check your ligand database and/or your assembly input file.')
#
#         i = -1
#         chosen_ligand_combinations = set()
#         n_combinations = len(all_combs)
#         while i < n_combinations:
#             i += 1      # infinite for loop
#             try:
#                 ligands = all_combs[i]
#             except IndexError:
#                 self.stop_assembly = True
#                 yield None
#
#             assert [lig.denticity for lig in ligands] == self.topology, f"Topology of ligands {ligands} does not match the desired topology {self.topology}! This should not happen!"
#
#             ligands = {idx: lig for idx, lig in enumerate(ligands)}
#             ligands = self.format_similarity_lists(ligands, self.instruction)  # replace ligands with the same similarity with the same ligand
#             ligand_names = tuple(sorted([ligand.unique_name for ligand in ligands.values()]))
#
#             if ligand_names in chosen_ligand_combinations:
#                 # This combination has already been chosen
#                 continue
#             else:
#                 chosen_ligand_combinations.add(ligand_names)
#
#             # Check total charge condition
#             sum_ligand_charges = sum([ligand.pred_charge for ligand in ligands.values()])
#             if sum_ligand_charges == self.total_charge - self.metal_ox:
#                 # Total charge condition is satisfied
#                 ligands = {idx: deepcopy(ligand) for idx, ligand in ligands.items()}
#                 yield ligands

import random
import itertools
from pathlib import Path
from typing import Union
import logging
from DARTassembler.src.metalig.mol import Ligand


class LigandChoice(object):

    def __init__(self, ligand_dbs, total_ligand_charges: int, n_max_complexes: Union[int,str]):
        """
        This class is used to choose ligands for the assembly of complexes. It supports both random and iterative ligand choice methods.
        :param ligand_dbs: List of ligand databases. Each database can be a LigandDB object or 'same_as_previous' to indicate that the same ligand as in the previous position should be used.
        :param total_ligand_charges: Total charge of the ligands in the complex. This is used to filter out ligand combinations which do not fulfill this constraint.
        :param n_max_complexes: Maximum number of complexes to be assembled. If 'all', all possible complexes will be assembled.
        """
        self.ligand_dbs = ligand_dbs
        self.total_ligand_charges = total_ligand_charges
        self.n_max_complexes = n_max_complexes  # int or 'all'
        self.ligand_choice = 'all' if n_max_complexes == 'all' else 'random'
        self.continue_assembly = True   # If set to False, the assembler will stop
        self.max_rejected_ligand_combinations = 1_000  # If this many ligand combinations are rejected in a row, the random choice of ligands is exhausted and will be switched to iterative mode
        self.switched_to_iterative = False  # Flag to output a warning if the ligand choice method was switched to "all"

    def _check_good_ligand_charges(self, ligand_combination: list):
        """
        Check if these ligands have the correct sum of charges.
        """
        if self.total_ligand_charges is None:
            correct_charges = True
        else:
            total_ligand_charges = sum([ligand.charge for ligand in ligand_combination])
            correct_charges = total_ligand_charges == self.total_ligand_charges
        return correct_charges

    def if_make_more_complexes(self, num_assembled_complexes: int) -> bool:
        """
        Returns True if more complexes can be assembled, False otherwise.
        """
        if self.n_max_complexes == "all":
            return self.continue_assembly
        else: # self.n_max_complexes is an integer
            return num_assembled_complexes < self.n_max_complexes

    def _choose_random_ligand_combination_from_db(self) -> list:
        """
        Choose ligands randomly from the ligand databases.
        """
        ligand_combination = []
        for ligand_list in self.ligand_dbs:
            # Choose ligands randomly and respect the "same_as_previous" instruction
            if ligand_list == 'same_as_previous':
                chosen_ligand = ligand_combination[-1]
            else:
                ligand_list = list(ligand_list.db.values())
                chosen_ligand = random.choice(ligand_list)
            ligand_combination.append(chosen_ligand)

        return ligand_combination

    def _choose_iterative_ligand_combination_from_db(self, all_combinations: list[tuple[str]]) -> list:
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
        for idx, ligand_list in enumerate(self.ligand_dbs):
            if ligand_list == 'same_as_previous':
                chosen_ligand = prel_ligand_combination[-1]
            else:
                chosen_ligand = prel_ligand_combination[idx]
            # Fixed a bug where this function would return the name of the ligand rather than the ligand object itself
            ligand_combination.append(ligand_list.db[chosen_ligand])

        return ligand_combination

    def _final_assertions_for_ligand_combination(self, ligand_combination: list[Ligand]) -> None:
        """
        Doublechecks the ligand combination for consistency with all constraints.
        """
        # Charges
        if self.total_ligand_charges is not None:
            sum_of_charges = sum([ligand.charge for ligand in ligand_combination])
            assert sum_of_charges == self.total_ligand_charges, f"The sum of charges {sum_of_charges} of the ligand combination {ligand_combination} is not equal to the total_ligand_charges {self.total_ligand_charges}."

        for idx, ligand in enumerate(ligand_combination):
            # Correct `same_as_previous`
            if self.ligand_dbs[idx] == 'same_as_previous':
                assert ligand.unique_name == ligand_combination[idx-1].unique_name, f"The ligand {ligand.unique_name} at index {idx} is not the same as the ligand {ligand_combination[idx-1].unique_name} at index {idx-1}!"

        return

    def choose_ligands(self) -> Union[dict,None]:
        """
        Choose ligands for the assembly of complexes. This function is a generator and yields ligand combinations. There are two modes, random and iterative, which are chosen by the ligand_choice attribute.
        - ligand_choice = 'all': If the mode is iterative, it's very simple: All structures will be made, though the order is non-random.
        - ligand_choice = 'random': If the mode is random, the function will yield a random ligand combination each time it is called. The function will stop yielding ligand combinations if the maximum number of complexes has been reached or if no more valid ligand combinations can be found. In the latter case, the random mode will also switch to iterative mode to make sure all possible complexes are made.
        """
        # Setup all ligand combinations as iterable. Needed for the iterative ligand choice method.
        assert self.ligand_dbs[-1] == 'same_as_previous' if 'same_as_previous' in self.ligand_dbs else True, "The 'same_as_previous' instruction must always come last in the list of ligand lists!" # HARDCODED: If the 'same_as_previous' instruction is used, it always comes last in the list of ligand lists
        all_valid_lists = [ligands.db for ligands in self.ligand_dbs if ligands != 'same_as_previous']
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
                count_rejected_ligand_combinations_in_a_row = 0     # novel ligand comb, reset counter
                chosen_ligand_combinations.add(ligand_names)

            # Final check if the ligand combination fulfills all constraints
            self._final_assertions_for_ligand_combination(ligand_combination)

            yield ligand_combination

        if len(chosen_ligand_combinations) == 0: # Output error because no valid ligand combinations found
            raise LigandCombinationError(
                f'No valid ligand combinations found which fulfill the total_ligand_charges {self.total_ligand_charges} requirement. Please check your ligand database and/or your assembly input file.')

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
                        f'The provided ligand database doesn\'t contain the denticity {dent} in the topology {topology} at the {idx + 1}th site! Please check your ligand database and/or your assembly input file')
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
        assert len(database) == len(topology), f"The number of topologies and the number of ligand databases are not consistent: {len(topology)} != {len(database)}"

        database = self._get_only_relevant_denticities_in_db(database=database, topology=topology)
        return database


class LigandCombinationError(Exception):
    """
    Exception raised for errors in the choice of ligands.
    """
    def __init__(self, message: str, file: str='', batch_name:str =''):
        file = Path(file).name
        if file != '':
            file = f" in input file '{file}'"
        total_message = f"\n\t--> Ligand Combination Error{batch_name}{file}:\n\t\t{message}"
        super().__init__(total_message)

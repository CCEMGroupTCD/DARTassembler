import warnings
from copy import deepcopy
import random
import itertools
from src05_Assembly_Refactor.Assemble import PlacementRotation


class ChooseRandomLigands:
    def __init__(self, database, topology, instruction, max_attempts, metal_oxidation_state: int = None, total_complex_charge: int = None):
        self.ligand_dict = database
        self.topology = topology
        self.instruction = instruction
        self.max_loop = max_attempts
        self.metal_ox = metal_oxidation_state
        self.total_charge = total_complex_charge

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

    def get_charge_dic(self):
        tmp_dic_1 = {}  # dictionary with keys (denticty) and values(list of all charges)
        tmp_dic_2 = {}  # dictionary with keys (denticty) and values(list of all ligands in the same order as their charges in tmp_dic_1)
        for denticity in self.ligand_dict.keys():
            tmp_charge_list = []
            tmp_ligand_list = []
            for ligand in self.ligand_dict[denticity]:
                try:
                    charge = ligand.pred_charge
                except:
                    charge = ligand.global_props["LCS_pred_charge"]
                tmp_charge_list.append(charge)
                tmp_ligand_list.append(ligand)
            tmp_dic_1.update({f"{denticity}": deepcopy(tmp_charge_list)})
            tmp_dic_2.update({f"{denticity}": deepcopy(tmp_ligand_list)})
        return tmp_dic_1, tmp_dic_2

    def get_ligand_dic(self, dic_1, dic_2):
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
            tmp_dic_3.update({f"{denticity}": deepcopy(tmp_dic_charge)})
        return tmp_dic_3

    def charge_list_process(self):
        print("\nStarting Charge Loop")
        m = 0
        charge_dic, ligand_dic = self.get_charge_dic()
        while m < self.max_loop:
            charge_list = []
            for dent in self.topology:
                charge_list.append(random.choice(charge_dic[str(dent)]))

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

    def choose_ligands(self):
        charge_list = self.charge_list_process()
        dic_1, dic_2 = self.get_charge_dic()
        ligand_dic = self.get_ligand_dic(dic_1=dic_1, dic_2=dic_2)
        ligands = {}
        i = 0
        if charge_list is None:
            return None
        else:
            pass
        for denticity, charge in zip(self.topology, charge_list):
            ligands.update({i: random.choice(ligand_dic[str(denticity)][str(charge)])})
            i = i + 1
        ligands_out = self.format_similarity_lists(ligands, self.instruction)
        return ligands_out


class ChooseIterativeLigands:

    def __init__(self, database, top_list, charge, metal, random_seed):
        self.top_list = top_list
        self.topology = None
        self.similarity_list = None
        self.ligand_dict = database
        self.metal_list = metal
        self.Total_Charge = charge
        self.Random_Seed = random_seed
        self.charge_dic, self.ligand_dic = self.get_charge_dic()
        self.perform_pre_checks()
        self.Combinations = None

    def get_charge_dic(self):
        tmp_dic_1 = {}  # dictionary with keys (denticity) and values(list of all charges)
        tmp_dic_2 = {}  # dictionary with keys (denticity) and values(list of all ligands in the same order as their charges in tmp_dic_1)
        for denticity in self.ligand_dict.keys():
            tmp_charge_list = []
            tmp_ligand_list = []
            for ligand in self.ligand_dict[denticity]:
                try:
                    charge = ligand.pred_charge
                except:
                    charge = ligand.global_props["LCS_pred_charge"]
                tmp_charge_list.append(charge)
                tmp_ligand_list.append(ligand)
            tmp_dic_1.update({f"{denticity}": deepcopy(tmp_charge_list)})
            tmp_dic_2.update({f"{denticity}": deepcopy(tmp_ligand_list)})
        return tmp_dic_1, tmp_dic_2

    def perform_pre_checks(self):
        if len(self.top_list) != 1:
            # Ensures only one topology provided
            print("!!!Fatal Error!!! -> In order to Assemble complexes deterministically only one topology can be provided -> Exiting Program")
            raise ValueError

        elif len(self.metal_list) != 1:
            # Ensures only one metal provided
            print("!!!Fatal Error!!! -> In order to Assemble complexes deterministically only one metal can be provided -> Exiting Program")
            raise ValueError

        self.topology, self.similarity_list = PlacementRotation.format_topologies(self.top_list[0])
        if sorted(set(self.topology)) != sorted(set([int(x) for x in self.charge_dic.keys()])):
            print("!!!Fatal Error!!! -> The ligands provided in the database are not consistent with the desired denticities in the input-> Exiting Program")
            raise ValueError

        for denticity, charge_list in self.charge_dic.items():
            if len(set(charge_list)) != 1:
                print(f"!!!Fatal Error!!! -> [{len(set(charge_list))}] different charges detected for denticity [{denticity}], This value should be 1 -> Exiting Program")
                raise ValueError
            else:
                pass

        if self.Random_Seed != 0:
            print(f"!!!Fatal Error!!! -> The random seed must be set to zero not [{self.Random_Seed}]-> Exiting Program")
            raise ValueError

        print("!!!Success!!! -> ChooseIterativeLigands pre-checks have been PASSED -> Progressing to Assembly")

    def obtain_all_combinations(self):
        tmp_list = []
        for denticity in self.topology:
            tmp_list.append(self.ligand_dict[denticity])
        combo_list = list(itertools.product(*tmp_list))
        self.Combinations = [list(ele) for ele in combo_list]
        self.eliminate_duplicates()
        ans = input(f"DART has determined that there are [{len(self.Combinations)}] possible complexes that could be assembled\n"
                    f"Would you like to proceed (y/n)")
        if (ans == "y") or (ans == "yes") or (ans == "Yes") or (ans == "Y") or (ans == "YES"):
            pass
        else:
            exit()

    def eliminate_duplicates(self):
        # The purpose of this function is to eliminate duplicates from the self.Combinations list
        # i.e. lists of ligands that are the same apart from order
        # It is the job of the isomer handler to build isomers not the job of the ligand decider
        tmp_list = []
        for ligand_combo_1 in self.Combinations:
            detected = False
            for ligand_combo_2 in self.Combinations:
                if set(ligand_combo_1) == set(ligand_combo_2) and ((self.Combinations.index(ligand_combo_1)) != (self.Combinations.index(ligand_combo_2))):
                    print("duplicate detected")
                    detected = True
                else:
                    pass
            if not detected:
                tmp_list.append(ligand_combo_1)
            else:
                pass
        self.Combinations = tmp_list

    def choose_ligands(self, iteration: int = None):
        tmp_dic = {}
        i = 0
        for ligand in self.Combinations[iteration]:
            tmp_dic.update({i: ligand})
            i = i + 1
        return tmp_dic

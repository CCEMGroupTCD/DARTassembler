import json
import os
import ast
import networkx as nx
from src01.graph_utility import graph_from_graph_dict, view_graph
from tqdm import tqdm
import matplotlib.pyplot as plt

sub_graph_iodine = {'graph': {'1': {}}, 'node_attributes': {'1': {'metal_neighbor': True, 'node_label': 'C', 'orig_idx': 1}}}
sub_graph_tri_pyridne = {'graph': {'1': {'2': {}, '10': {}}, '2': {'1': {}, '3': {}, '4': {}}, '3': {'2': {}}, '4': {'2': {}, '5': {}, '6': {}}, '5': {'4': {}},
                                   '6': {'4': {}, '7': {}, '8': {}}, '7': {'6': {}}, '8': {'6': {}, '9': {}, '10': {}}, '9': {'8': {}}, '10': {'1': {}, '8': {}, '13': {}},
                                   '11': {'13': {}, '16': {}}, '12': {'13': {}, '14': {}}, '13': {'10': {}, '11': {}, '12': {}}, '14': {'12': {}, '15': {}, '39': {}},
                                   '15': {'14': {}, '16': {}}, '16': {'11': {}, '15': {}, '18': {}}, '17': {'18': {}, '25': {}}, '18': {'16': {}, '17': {}, '19': {}},
                                   '19': {'18': {}, '20': {}, '21': {}}, '20': {'19': {}}, '21': {'19': {}, '22': {}, '23': {}}, '22': {'21': {}}, '23': {'21': {}, '24': {}, '25': {}},
                                   '24': {'23': {}}, '25': {'17': {}, '23': {}, '26': {}}, '26': {'25': {}}, '39': {'14': {}, '40': {}, '41': {}}, '40': {'39': {}}, '41': {'39': {}}},
                         'node_attributes': {'1': {'metal_neighbor': True, 'node_label': 'N', 'orig_idx': 1}, '2': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 2},
                                             '3': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 3}, '4': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 4},
                                             '5': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 5}, '6': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 6},
                                             '7': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 7}, '8': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 8},
                                             '9': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 9}, '10': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 10},
                                             '11': {'metal_neighbor': True, 'node_label': 'N', 'orig_idx': 11}, '12': {'metal_neighbor': False, 'node_label': 'N', 'orig_idx': 12},
                                             '13': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 13}, '14': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 14},
                                             '15': {'metal_neighbor': False, 'node_label': 'N', 'orig_idx': 15}, '16': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 16},
                                             '17': {'metal_neighbor': True, 'node_label': 'N', 'orig_idx': 17}, '18': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 18},
                                             '19': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 19}, '20': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 20},
                                             '21': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 21}, '22': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 22},
                                             '23': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 23}, '24': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 24},
                                             '25': {'metal_neighbor': False, 'node_label': 'C', 'orig_idx': 25}, '26': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 26},
                                             '39': {'metal_neighbor': False, 'node_label': 'N', 'orig_idx': 39}, '40': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 40},
                                             '41': {'metal_neighbor': False, 'node_label': 'H', 'orig_idx': 41}}}

test_input_dict = {"lig_atoms_must_include": str(["C", "H", "N"]),
                   "lig_atoms_must_not_include": str(["S", "H", "Cl"]),
                   "coord_atoms_must_include": str(["N"]),
                   "coord_atoms_must_not_include": str(["C"]),
                   "sub_graph": str(sub_graph_iodine)}


class JsonFilter:
    def __init__(self, input_dict):
        self.input_dict = input_dict
        self.input_file = "/Users/cianclarke/Documents/PhD/Complex_Assembly/Data/tmQM_Ligands_unique_v1.2.json"
        self.output_file = "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/output_test.json"
        ###
        # The various inputs from the .json file
        self.unique_name = None
        self.ligand_info = None
        self.stoichiometry = None
        self.atomic_props = None
        self.global_props = None
        self.graph_dict = None
        self.denticity = None
        self.ligand_to_metal = None
        self.local_elements = None
        self.name = None
        self.CSD_code = None
        self.graph_hash = None
        self.occurrences = None
        self.count_denticities = None
        self.count_metals = None
        self.n_denticities = None
        self.n_metals = None
        self.all_ligand_names = None
        self.pred_charge = None
        self.pred_charge_is_confident = None

        if self.input_dict["lig_atoms_must_include"] is not None:
            print("filter_1 active")
            pass

        if self.input_dict["lig_atoms_must_not_include"] is not None:
            print("filter_2 active")
            pass

        if self.input_dict["coord_atoms_must_include"] is not None:
            print("filter_3 active")
            pass

        if self.input_dict["coord_atoms_must_not_include"] is not None:
            print("filter_4 active")
            pass

        if self.input_dict["sub_graph"] is not None:
            print("filter_5 active")
            pass

    def Iterator(self):
        # This function iterates through the entire json file
        with open(self.input_file) as data_file:
            data = json.load(data_file)
            for unique_name, ligand_info in tqdm(zip(data.keys(), data.values())):
                self.unique_name = unique_name
                self.ligand_info = ligand_info
                print("The unique name is: " + str(unique_name))
                for ligand_info_key, ligand_info_value in zip(self.ligand_info.keys(), self.ligand_info.values()):
                    ###
                    print("ligand_info_key: " + str(ligand_info_key))
                    print("ligand_info_value: " + str(ligand_info_value))
                    self.stoichiometry = self.ligand_info["stoichiometry"]
                    self.atomic_props = self.ligand_info["atomic_props"]
                    self.global_props = self.ligand_info["global_props"]
                    self.graph_dict = self.ligand_info["graph_dict"]
                    self.denticity = self.ligand_info["denticity"]
                    self.ligand_to_metal = self.ligand_info["ligand_to_metal"]
                    self.local_elements = self.ligand_info["local_elements"]
                    self.name = self.ligand_info["name"]
                    self.CSD_code = self.ligand_info["CSD_code"]
                    self.graph_hash = self.ligand_info["graph_hash"]
                    self.occurrences = self.ligand_info["occurrences"]
                    self.count_denticities = self.ligand_info["count_denticities"]
                    self.count_metals = self.ligand_info["count_metals"]
                    self.n_denticities = self.ligand_info["n_denticities"]
                    self.n_metals = self.ligand_info["n_metals"]
                    self.all_ligand_names = self.ligand_info["all_ligand_names"]
                    self.pred_charge = self.ligand_info["pred_charge"]
                    self.pred_charge_is_confident = self.ligand_info["pred_charge_is_confident"]
                    ###

                if self.input_dict["lig_atoms_must_include"] is not None:
                    print(self.filter_including_atoms())
                    pass

                if self.input_dict["lig_atoms_must_not_include"] is not None:
                    print(self.filter_excluding_atoms())

                if self.input_dict["coord_atoms_must_include"] is not None:
                    print(self.filter_including_atoms_coord())
                    pass

                if self.input_dict["coord_atoms_must_not_include"] is not None:
                    print(self.filter_excluding_atoms_coord())
                    pass

                if self.input_dict["sub_graph"] is not None:
                    print(self.filter_sub_graph_search())
                    if self.filter_sub_graph_search():
                        self.view_ligand()
                    pass

    def view_ligand(self):
        # This utility essentially just allows us to view the various ligands
        view_file = "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/test/filter_view.xyz"
        os.system(f'touch {view_file}')
        with open(view_file, 'w') as f:
            f.write(str(len(self.atomic_props["atoms"])))
            f.write("\n\n")
            for atom, x, y, z in zip(self.atomic_props["atoms"], self.atomic_props["x"], self.atomic_props["y"], self.atomic_props["z"]):
                f.write(f'{atom}   {x}    {y}    {z} \n')
        os.system(f'ase-gui {view_file}')
        next_ = input("press enter to continue")

    def filter_including_atoms(self):
        # This filter will only allow ligands past if and only if they contain all the atoms selected by the user
        input_list = ast.literal_eval(self.input_dict["lig_atoms_must_include"])
        matched_atoms = []
        for atom in self.atomic_props["atoms"]:
            if atom in input_list:
                matched_atoms.append(atom)
        if set(matched_atoms) == set(input_list):
            return "pass"
        else:
            return "fail"

    def filter_excluding_atoms(self):
        # This filter will only allow ligands past if and only if none of the atoms selected by the user are contained in the ligand
        input_list = ast.literal_eval(self.input_dict["lig_atoms_must_not_include"])
        for atom in self.atomic_props["atoms"]:
            if atom in input_list:
                return "fail"
        return "pass"

    def filter_including_atoms_coord(self):
        # This filter will only allow ligands past if and only if the coordinating groups contain at least one the atoms selected by the user
        input_list = ast.literal_eval(self.input_dict["coord_atoms_must_include"])
        matched_atoms = []
        for atom in self.local_elements:
            if atom in input_list:
                matched_atoms.append(atom)
        if set(matched_atoms) == set(input_list):
            return "pass"
        else:
            return "fail"
        pass

    def filter_excluding_atoms_coord(self):
        # This filter will only allow ligands past if and only if none of the atoms selected by the user are contained in the functional groups of the ligand
        input_list = ast.literal_eval(self.input_dict["coord_atoms_must_not_include"])
        for atom in self.local_elements:
            if atom in input_list:
                return "fail"
        return "pass"

    @staticmethod
    def node_check(dict1, dict2):
        return dict1["node_label"] == dict2["node_label"]

    def filter_sub_graph_search(self):
        ligand_graph = graph_from_graph_dict(self.graph_dict)
        sub_graph = graph_from_graph_dict(ast.literal_eval(self.input_dict["sub_graph"]))
        a = nx.algorithms.isomorphism.ISMAGS(ligand_graph,sub_graph,node_match=self.node_check).is_isomorphic()
        return a

instance = JsonFilter(test_input_dict)

instance.Iterator()

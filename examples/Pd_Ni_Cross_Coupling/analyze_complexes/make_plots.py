import os
import shutil
import warnings
import pymatgen.core as pt
import re
import pickle
import numpy as np
from typing import Union, Tuple
import json
from tqdm import tqdm
import random
import pandas as pd
import ase
import networkx as nx
import matplotlib
from pathlib import Path
from DARTassembler.src.constants.Periodic_Table import DART_Element
import seaborn as sns
from DARTassembler.src.ligand_extraction.io_custom import load_json
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import cclib


class GaussianCalculation:
    def __init__(self, calc_path, expected_donors):
        self.calc_name = Path(calc_path).name
        self.calc_path = calc_path
        self.expected_donors = expected_donors
        self.file_path_list = self.get_files()
        ###Keeps track of how many times a calculation has been restarted
        self.restart_status = self.calculation_restarted()  # Dictionary That contains info on how many restarts
        self.restart_tag = self.get_continuation_stamp()  # i.e. _cont or _cont_cont or _cont_cont_cont
        self.com, self.log, self.chk, self.pkl, = self.get_latest_files()
        self.fchk = self.chk.replace(".chk", ".fchk") if self.chk is not None else None
        self.fchk_parser = cclib.io.ccread(self.fchk)
        self.log_atoms, self.log_x, self.log_y, self.log_z = self.log_to_xyz()
        self.metal_idx, self.metal_symbol = self.get_metal_atom()
        self.metal_position = [self.log_x[self.metal_idx], self.log_y[self.metal_idx], self.log_z[self.metal_idx]]
        self.ase_mol = self.get_ase_mol()
        self.coordinates = self.ase_mol.positions
        self.n_donors = len(self.expected_donors)
        self.donor_idc, self.donor_elements, self.donor_distances = self.get_coordinating_atoms()
        self.homo, self.lumo, self.hl_gap = self.get_homo_lumo()
        self.metal_mulliken_charge = self.get_metal_mulliken_charge()
        # self.Au_C = None
        # self.Au_N = None
        # self.C_aromatic = None
        # self.N_aromatic = None
        # self.gold_binding_precedence = None
        # self.C_ring_hetero_atoms = []
        # self.N_ring_hetero_atoms = []
        # self.bidentate_bite_angle = None
        # self.get_NBO_data()
        # self.Au_C_bond_distance = None
        # self.Au_N_bond_distance = None
        # self.shortest_distance_between_donor_atoms = None
        # self.ring_info = None
        # self.donors_in_rings = None
        # self.total_Lewis = None
        # self.cpu_time_in_hours = float(self.get_cpu_hours())
        self.num_atoms = len(self.log_atoms)
        # self.N_ring_size = None
        # self.C_ring_size = None
        # self.Num_atoms = len(self.log_atoms)
        # self.CSD_list = None
        # self.load_pkl()
        self.molecular_weight = self.get_Molecular_weight()
        self.longest_distance = None
        self.get_longest_distance_between_atoms()
        assert len(self.log_atoms) == len(self.log_x) == len(self.log_y) == len(self.log_z)
        # assert self.Au_C is not None
        # assert self.Au_N is not None
        # assert self.molecular_weight is not None
        # assert self.longest_distance is not None
        # assert self.Au_C_bond_distance is not None
        # assert self.Au_N_bond_distance is not None
        # assert self.bidentate_bite_angle is not None
        # assert self.N_ring_size is not None
        # assert self.C_ring_size is not None
        # assert self.C_aromatic is not None
        # assert self.N_aromatic is not None
        # assert self.CSD_list is not None

    def get_bite_angle(self, atoms: list) -> float:
        """
        Calculates the bite angle between two coordinating atoms with the metal atom as the vertex.
        @atoms: list of two elements that are the coordinating atoms, e.g. ["P", "N"]
        """
        if not len(atoms) == 2:
            raise ValueError(f"Expected two atoms, but got {len(atoms)}")
        if not all([atom in self.donor_elements for atom in atoms]):
            raise ValueError(f"Expected two donor atoms, but got {atoms}")

        donor_idc = [idx for idx, el in zip(self.donor_idc, self.donor_elements) if el in atoms]
        positions = [self.coordinates[idx] for idx in donor_idc] + [self.metal_position]

        return self.angle_between_three_points(*positions)



    def get_ase_mol(self):
        atomic_props = {"x": self.log_x, "y": self.log_y, "z": self.log_z, "atoms": self.log_atoms}
        coord_list_3D = [[atomic_props[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
                         enumerate(atomic_props["x"])]
        atom_list = atomic_props["atoms"]
        mol = ase.Atoms(atom_list, positions=coord_list_3D)
        return mol

    def get_metal_atom(self) -> tuple[int, str]:
        """
        This function returns the metal atom as a tuple of the index and the element symbol.
        """
        metals = [(i, el) for i, el in enumerate(self.log_atoms) if DART_Element(el).is_metal]
        assert len(metals) == 1
        idx, el = metals[0]
        return idx, el

    def get_coordinating_atoms(self) -> tuple[list[int], list[str], list[float]]:
        """
        This function returns the indices, element symbols and distances of the coordinating atoms.
        @return: tuple of indices, element symbols, distances of the coordinating atoms
        """
        distances = self.ase_mol.get_all_distances()[self.metal_idx]
        # donor_dists = sorted(distances)[1:self.n_donors+1]
        # assert len(donor_dists) == self.n_donors, f"Expected {self.n_donors} donor atoms, but got {len(donor_dists)}"
        # max_donor_distances_to_metal = max(donor_dists)# exclude metal which will be the first element after sorting
        atoms = [(idx, el, dist) for idx, el, dist in zip(range(len(self.log_atoms)), self.log_atoms, distances)]
        atoms = sorted(atoms, key=lambda x: x[2])

        donors_not_yet_found = self.expected_donors.copy()
        donors = []
        for idx, el, dist in atoms:
            if el in donors_not_yet_found:
                donors.append((idx, el, dist))
                donors_not_yet_found.remove(el)

        #
        # idc = []
        # elements = []
        # donor_distances = []
        # for donor in self.expected_donors:
        #     for i, (el, distance) in enumerate(zip(self.log_atoms, distances)):
        #         if el == donor and distance <= max_donor_distances_to_metal and i != self.metal_idx:
        #             idc.append(i)
        #             elements.append(el)
        #             donor_distances.append(distance)
        # a = 1
        idc, elements, donor_distances = zip(*donors)
        a = 1
        assert sorted(elements) == sorted(self.expected_donors), f"Expected donors: {self.expected_donors}, but got: {elements}"

        return idc, elements, donor_distances




    def get_files(self):
        file_list = os.listdir(self.calc_path)
        file_path_list = []
        for file in file_list:
            file_path_list.append(self.calc_path + "/" + file)
        return file_path_list

    def get_latest_files(self):
        tmp_com = None
        tmp_log = None
        tmp_pkl = None
        tmp_chk = None
        for file_path in self.file_path_list:
            if (str(file_path).endswith(".com")) and (str(file_path).count("_cont") == self.restart_status[".com"]):
                tmp_com = file_path
            elif (str(file_path).endswith(".log")) and (str(file_path).count("_cont") == self.restart_status[".log"]):
                tmp_log = file_path
            elif (str(file_path).endswith(".json")) and (str(file_path).count("_cont") == self.restart_status[".json"]):
                tmp_pkl = file_path
            elif (str(file_path).endswith(".chk")) and (str(file_path).count("_cont") == self.restart_status[".chk"]):
                tmp_chk = file_path
            else:
                pass
        assert (tmp_com is not None) and (tmp_log is not None) and (tmp_pkl is not None) and (tmp_chk is not None)
        return tmp_com, tmp_log, tmp_chk, tmp_pkl

    def calculation_restarted(self):
        # This also serves as a check
        # to make sure all the correct
        # files are present
        tmp = {}
        com_detected = []
        log_detected = []
        chk_detected = []
        pkl_detected = []
        for file_path in self.file_path_list:
            if str(file_path).endswith(".com"):
                com_detected.append(str(file_path).count("_cont"))
            elif str(file_path).endswith(".log"):
                log_detected.append(str(file_path).count("_cont"))
            elif str(file_path).endswith(".json"):
                pkl_detected.append(str(file_path).count("_cont"))
            elif str(file_path).endswith(".chk"):
                chk_detected.append(str(file_path).count("_cont"))
        tmp.update({".com": max(com_detected),
                    ".log": max(log_detected),
                    ".json": max(pkl_detected),
                    ".chk": max(chk_detected)})
        assert (len(com_detected) != 0) and (len(chk_detected) != 0) and (len(log_detected) != 0) and (len(pkl_detected) != 0)
        # assert tmp[".com"] == tmp[".chk"]
        return tmp

    def get_continuation_stamp(self):
        stamp = "_cont"
        final_stamp = ""
        for i in range(self.restart_status[".com"] + 1):
            final_stamp = final_stamp + stamp
        return final_stamp

    def log_to_xyz(self):
        with open(self.log, 'r') as log_file:
            rline = log_file.readlines()

        detected_start_coord = False
        already_got_coords = False
        symbol_list = []
        x_coord_list = []
        y_coord_list = []
        z_coord_list = []
        i = 0
        Standard_orientation_count_total = 0
        Standard_orientation_count = 0
        for line in rline:
            if "                        Standard orientation:                        " in line:
                Standard_orientation_count_total += 1
            else:
                pass

        for line in rline:
            if "                        Standard orientation:                        " in line:
                Standard_orientation_count += 1
            else:
                pass

            if "---------------------------------------------------------------------" in line:
                detected_start_coord = False

            if detected_start_coord:
                # This regex finds all numbers in a string and puts them in a list
                atom_list = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", str(line))
                atomic_number = int(atom_list[1])
                x_coord = float(atom_list[3])
                y_coord = float(atom_list[4])
                z_coord = float(atom_list[5].strip("\n"))
                symbol_list.append(pt.Element.from_Z(atomic_number).symbol)
                x_coord_list.append(x_coord)
                y_coord_list.append(y_coord)
                z_coord_list.append(z_coord)
                # essentially this line will pull the last set of coordinates from the .log file
            if ("---------------------------------------------------------------------" in line) and \
                    ("Number     Number       Type             X           Y           Z" in rline[i - 1]) and \
                    ("Center     Atomic      Atomic             Coordinates (Angstroms)" in rline[i - 2]) and \
                    ("---------------------------------------------------------------------" in rline[i - 3]) and \
                    ("                        Standard orientation:                        " in rline[i - 4]) and not already_got_coords and (
                    (Standard_orientation_count_total - Standard_orientation_count == 1) or Standard_orientation_count_total == 1):
                already_got_coords = True
                detected_start_coord = True
            else:
                pass
            i = i + 1
        return symbol_list, x_coord_list, y_coord_list, z_coord_list

    def get_NBO_data(self):
        with open(self.log, "r") as f:
            lines = f.readlines()

        start = "SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS"
        end = "NATURAL BOND ORBITALS (Summary):"
        start_found = False
        end_found = False
        result_N_found = False
        result_C_found = False

        carbon_interaction_sum = 0
        nitrogen_interaction_sum = None
        for idx, line in enumerate(lines):

            if start in line:
                start_found = True
            if end in line:
                end_found = True

            # NOTE: to save me having to loop over the file again I am going to extract the lewis acid data here too
            if ('Total Lewis' in line) and (line.endswith(")")) and ('%' in line):
                print("haven't completed this bit yet")
                raise ValueError

            if start_found and not end_found:
                line_list = line.split()

                # NITROGEN
                # Here we make sure that LP comes before N which comes before Au which should come before Cl
                if (")Au" in line) and ("-Cl" in line) and ("LP" in line) and (len(line_list) > 3) and ("N" in line) and (line.find("LP") < line.find("N")) and (line.find("N") < line.find("Au")) and (
                        line.find("N") < line.find("Cl")):
                    interaction_energy_N = line_list[-3]
                    if nitrogen_interaction_sum is None:
                        nitrogen_interaction_sum = float(interaction_energy_N)
                    if nitrogen_interaction_sum is not None:
                        if float(interaction_energy_N) > float(nitrogen_interaction_sum):
                            nitrogen_interaction_sum = float(interaction_energy_N)
                        else:
                            pass
                    result_N_found = True
                else:
                    pass

                """
                # NITROGEN
                # Here we make sure that LP comes before N which comes before Au which should come before Cl
                if (")Au" in line) and ("-Cl" in line) and ("LP" in line) and (len(line_list) > 3) and ("N" in line) and (line.find("LP") < line.find("N")) and (line.find("N") < line.find("Au")) and (
                        line.find("N") < line.find("Cl")):
                    interaction_energy_N = line_list[-3]
                    nitrogen_interaction_sum = interaction_energy_N
                    # print(line.strip() + f"         Interaction_N: [{interaction_energy_N}]")
                    result_N_found = True
                else:
                    pass"""

                # CARBON
                if (")Au" in line) and ("BD" in line) and ("RY" in line) and (len(line_list) > 3) and ("- C" in line) and (") C" not in line) and ("Cl" not in line) and ("N" not in line) and (
                        line.find("BD") < line.find("- C")):
                    interaction_energy_C = float(line_list[-3])

                    if interaction_energy_C >= 5.0:
                        carbon_interaction_sum += interaction_energy_C
                    # print(line.strip() + f"         Interaction_C: [{interaction_energy_C}]")
                    result_C_found = True
                else:
                    pass

        if (not result_N_found) or (not result_C_found):
            self.Au_N = 1000 + round(random.uniform(0, 100), 6)
            self.Au_C = 1000 + round(random.uniform(0, 100), 6)
            print(f"!!!Warning!!!-> No value found for Au-N or Au_C -> [{self.calc_name}]")
            # raise ValueError

        else:
            self.Au_N = nitrogen_interaction_sum
            self.Au_C = carbon_interaction_sum

    def load_pkl(self):
        file = open(self.pkl, 'rb')
        ligand = pickle.load(file)
        elem = ligand.local_elements
        indexes = ligand.ligand_to_metal
        metal_coords = np.array([self.log_x[0], self.log_y[0], self.log_z[0]])
        carbon_atom_coords = np.array([self.log_x[indexes[elem.index('C')] + 1], self.log_y[indexes[elem.index('C')] + 1], self.log_z[indexes[elem.index('C')] + 1]])
        nitrogen_atom_coords = np.array([self.log_x[indexes[elem.index('N')] + 1], self.log_y[indexes[elem.index('N')] + 1], self.log_z[indexes[elem.index('N')] + 1]])
        Au_N_bond_lenght = np.linalg.norm(metal_coords - nitrogen_atom_coords)
        Au_C_bond_lenght = np.linalg.norm(metal_coords - carbon_atom_coords)

        self.Au_N_bond_distance = Au_N_bond_lenght
        self.Au_C_bond_distance = Au_C_bond_lenght

        self.bidentate_bite_angle = self.angle_between_three_points(carbon_atom_coords, metal_coords, nitrogen_atom_coords)

        # source donor atom
        source_atomic_index = ligand.ligand_to_metal[0]
        source_atomic_symbol = ligand.local_elements[0]
        # target donor atom
        target_atomic_index = ligand.ligand_to_metal[1]
        target_atomic_symbol = ligand.local_elements[1]
        shortest_path_lenght = nx.shortest_path_length(ligand.graph,
                                                       source=ligand.graph_index_to_atomic_index[source_atomic_index],
                                                       target=ligand.graph_index_to_atomic_index[target_atomic_index])
        self.shortest_distance_between_donor_atoms = shortest_path_lenght

        # These lists contain int values of the size of each ring that each donor atom is part of
        C_ring_sizes = []
        N_ring_sizes = []

        #In this section of code I want to get a list of all the csd codes that this ligand has been seen in.
        all_csd_codes = []
        for ligand_name in ligand.all_ligand_names:
            CSD_code = str(ligand_name).split("-")[1]
            all_csd_codes.append(CSD_code)
        self.CSD_list = list(set(all_csd_codes)) # this is a unique list of all the complexes that contain this ligand

        # Here we loop through all the rings
        for cycle in nx.minimum_cycle_basis(ligand.graph):
            # for cycle in nx.cycle_basis(ligand.graph):
            # Then we loop through all the nodes
            for node in cycle:
                # If the source node is contained in the ring
                if int(ligand.graph_index_to_atomic_index[source_atomic_index]) == int(node):
                    if source_atomic_symbol == "N":
                        N_ring_sizes.append(cycle)
                        # self.N_ring_size = len(cycle)
                    elif source_atomic_symbol == "C":
                        C_ring_sizes.append(cycle)
                        # self.C_ring_size = len(cycle)
                    else:
                        warnings.warn("!!!Warning!!! --> the source or the target are not the donor atoms --> Exiting Program")
                        raise ValueError

                # If the target node is contained in the ring
                if int(ligand.graph_index_to_atomic_index[target_atomic_index]) == int(node):
                    if target_atomic_symbol == "N":
                        N_ring_sizes.append(cycle)
                        # self.N_ring_size = len(cycle)
                    elif target_atomic_symbol == "C":
                        C_ring_sizes.append(cycle)
                        # self.C_ring_size = len(cycle)
                    else:
                        warnings.warn("!!!Warning!!! --> the source or the target are not the donor atoms --> Exiting Program")
                        raise ValueError

        if len(C_ring_sizes) == 0:
            self.C_ring_size = 0
            self.C_aromatic = False
        else:
            self.C_ring_size = len(min(C_ring_sizes, key=len))
            not_aromatic = False
            for node in min(C_ring_sizes, key=len):
                node_degree = ligand.graph.degree[node]
                # If the node is not one of the donor atoms
                if (node != ligand.graph_index_to_atomic_index[target_atomic_index]) and (node != ligand.graph_index_to_atomic_index[source_atomic_index]):
                    #if a carbon atom in the ring has exactly 3 connections then it is aromatic
                    if (node_degree == 3) and (ligand.graph.nodes[node]['node_label'] == 'C'):
                        pass
                    #If there is a hetro atom it is a little more ambiguous but we still just assume its aromatic
                    elif ligand.graph.nodes[node]['node_label'] != 'C':
                        self.C_ring_hetero_atoms.append(ligand.graph.nodes[node]['node_label'])
                    #If none of these conditions are met then a carbon atom has more than 3 bonds and is no longer aromatic
                    else:
                        not_aromatic = True
                else:
                    pass
            if not_aromatic:
                self.C_aromatic = False
            else:
                self.C_aromatic = True




        if len(N_ring_sizes) == 0:
            self.N_ring_size = 0
            self.N_aromatic = False
        else:
            self.N_ring_size = len(min(N_ring_sizes, key=len))
            not_aromatic = False
            for node in min(N_ring_sizes, key=len):
                node_degree = ligand.graph.degree[node]
                # If the node is not one of the donor atoms
                if (node != ligand.graph_index_to_atomic_index[target_atomic_index]) and (node != ligand.graph_index_to_atomic_index[source_atomic_index]):
                    if (node_degree == 3) and (ligand.graph.nodes[node]['node_label'] == 'C'):
                        pass
                    elif ligand.graph.nodes[node]['node_label'] != 'C':
                        self.N_ring_hetero_atoms.append(ligand.graph.nodes[node]['node_label'])
                    else:
                        not_aromatic = True
                else:
                    pass
                if not_aromatic:
                    self.N_aromatic = False
                else:
                    self.N_aromatic = True


        for metal in ligand.count_metals.keys():
            if metal == "Au":
                self.gold_binding_precedence = True
            else:
                pass
        if self.gold_binding_precedence is None:
            self.gold_binding_precedence = False

    def get_cpu_hours(self):
        with open(self.log, "r") as f:
            lines = f.readlines()

        cpu_time_found = False
        cpu_time_total_in_hours = None
        for idx, line in enumerate(lines):
            if "Job cpu time:" in str(line):
                if not cpu_time_found:
                    words = str(line).split()
                    num_days = float(words[3])
                    num_hours = float(words[5])
                    num_minutes = float(words[7])
                    num_seconds = float(words[9])
                    cpu_time_total_in_hours = (num_days * 24.0) + num_hours + (num_minutes / 60.0) + (num_seconds / (60 * 60))
                    pass
                cpu_time_found = True
        assert cpu_time_total_in_hours is not None
        return cpu_time_total_in_hours

    @staticmethod
    def angle_between_three_points(point1, point2, point3):
        # Calculate vectors between the three points
        ba = point1 - point2
        bc = point3 - point2
        # Calculate angle between the vectors in radians
        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        return np.degrees(angle)

    def get_Molecular_weight(self):
        MW = 0
        for atom in self.log_atoms:
            MW += pt.Element(atom).atomic_mass
        return MW

    def get_longest_distance_between_atoms(self):

        max_distance = 0
        for x1, y1, z1 in zip(self.log_x, self.log_y, self.log_z):
            for x2, y2, z2 in zip(self.log_x, self.log_y, self.log_z):
                p1 = np.array([x1, y1, z1])
                p2 = np.array([x2, y2, z2])
                distance = np.linalg.norm(p2 - p1)
                if distance > max_distance:
                    max_distance = distance
                else:
                    pass
        self.longest_distance = max_distance

    @staticmethod
    def is_finished_calc_dir(calc_path):
        """
        Checks if the directory contains a calculation.
        """
        files = os.listdir(calc_path)
        if any([file.endswith(".com") for file in files]) and \
                any([file.endswith(".log") for file in files]) and \
                any([file.endswith(".chk") for file in files]):
            return True

    def get_homo_lumo(self) -> Tuple[float, float, float]:
        """
        Returns the HOMO, LUMO and HOMO-LUMO gap.
        """
        try:
            # Get HOMO and LUMO energies
            homo_index = self.fchk_parser.homos[0]
            lumo_index = homo_index + 1  # The orbital just above the HOMO
            homo_energy = self.fchk_parser.moenergies[0][homo_index]
            lumo_energy = self.fchk_parser.moenergies[0][lumo_index]

            # Calculate HOMO-LUMO gap
            gap = lumo_energy - homo_energy
        except AttributeError as e:
            print(f'Encountered error {e} in calculation {self.calc_path}')
            homo_energy, lumo_energy, gap = np.nan, np.nan, np.nan

        return homo_energy, lumo_energy, gap

    def get_metal_mulliken_charge(self):
        # Make sure charges have been calculated; you may use 'mulliken' or other methods depending on your calculation
        try:
            if 'mulliken' in self.fchk_parser.atomcharges:
                mulliken_charges = self.fchk_parser.atomcharges['mulliken']

                # Get the index of the metal atom (Pd or Ni)
                atom_numbers = self.fchk_parser.atomnos
                metal_index = [idx for idx, Z in enumerate(atom_numbers) if DART_Element(Z).is_transition_metal][0]

                # Get the charge on the metal atom
                metal_charge = mulliken_charges[metal_index]
            else:
                raise ValueError('Mulliken charges have not been found.')
        except AttributeError as e:
            print(f'Encountered error {e} in calculation {self.calc_path}')
            metal_charge = np.nan

        return metal_charge


if __name__ == "__main__":
    root = "../data_from_Cian"
    expected_donors = ['P', 'N', 'Br', 'C']
    data = []

    calcdirs = [calcdir for calcdir, _, _ in os.walk(root) if GaussianCalculation.is_finished_calc_dir(calcdir)]
    # calcdirs = calcdirs[:50]  # TODO: DEBUGGING

    # for calcdir in tqdm(calcdirs):
    #
    #     calc = Calculation(calc_path=calcdir, expected_donors=expected_donors)
    #     calc_data = {'metal': calc.metal_symbol}
    #     calc_data.update({f'dist_{donor}': distance for donor, distance in zip(calc.donor_elements, calc.donor_distances)})
    #     calc_data.update({
    #                         'P-N bite angle': calc.get_bite_angle(['P', 'N']),
    #                         'HOMO': calc.homo,
    #                         'LUMO': calc.lumo,
    #                         'HL_gap': calc.hl_gap,
    #                         'Metal charge': calc.metal_mulliken_charge,
    #                         'dir': calcdir
    #                         })
    #     calc_data.update({})
    #     data.append(calc_data)
    #     # NBO.update({calc.calc_name: [calc.Au_N, calc.Au_C]})
    #     # NBO.update({calc.calc_name: {'HTS_Data': {'X-value': calc.Au_N,
    #     #                                           'Y-value': calc.Au_C},
    #     #
    #     #                              'Structure': {'symbol': calc.log_atoms,
    #     #                                            'x_value': calc.log_x,
    #     #                                            'y_value': calc.log_y,
    #     #                                            'z_value': calc.log_z},
    #     #
    #     #                              'Other_Data': {'bond_distance [Au-N]:': calc.Au_N_bond_distance,
    #     #                                             'bond_distance [Au-C]:': calc.Au_C_bond_distance,
    #     #                                             'CSD_codes': calc.CSD_list,
    #     #                                             'Total_Lewis': 1,
    #     #                                             'Min_edges_between_donor_atoms': calc.shortest_distance_between_donor_atoms,
    #     #                                             'CPU_usage_in_cpu_hours': calc.cpu_time_in_hours,
    #     #                                             'Num_Atoms': calc.num_atoms,
    #     #                                             'Bite_Angle': calc.bidentate_bite_angle,
    #     #                                             'C_ring_size': calc.C_ring_size,
    #     #                                             'Is_C_ring_aromatic': calc.C_aromatic,
    #     #                                             'Is_N_ring_aromatic': calc.N_aromatic,
    #     #                                             'C_hetro_atom_list': calc.C_ring_hetero_atoms,
    #     #                                             'N_hetro_atom_list': calc.N_ring_hetero_atoms,
    #     #                                             'N_ring_size': calc.N_ring_size,
    #     #                                             'Num_atoms': calc.num_atoms,
    #     #                                             'Molecular_Weight': calc.molecular_weight,
    #     #                                             'longest distance': calc.longest_distance,
    #     #                                             'Gold_binding_precedence': calc.gold_binding_precedence}
    #     #                              }
    #     #             })
    # df = pd.DataFrame(data)
    # df.to_csv('../stats/data.csv', index=False)

    # %% Plot histograms
    sns.set_style("ticks")
    df = pd.read_csv('../stats/data.csv')

    # Remove complexes with NaN values
    df = df.dropna(axis='rows')

    for donor in expected_donors:
        df[f'dist_norm_{donor}'] = df[f'dist_{donor}'] - df['metal'].apply(lambda el: DART_Element(el).covalent_radius_angstrom) - DART_Element(donor).covalent_radius_angstrom

    # Population plots
    norm_distances = [f'norm_{donor}' for donor in expected_donors]
    others = ['P-N bite angle', 'HOMO', 'LUMO', 'HL_gap', 'Metal charge']
    cols = expected_donors + norm_distances + others
    for donor in cols:
        name = f'dist_{donor}' if donor not in others else donor
        # data = df[[name, 'metal']]
        # data['bins'] = pd.cut(data[name], 20)
        # data = data.groupby(by=['metal', 'bins']).size().reset_index(name='count')
        # data = data.rename(columns={'bins': name})
        # data = data.pivot(index=name, columns='metal', values='count').reset_index()
        # data['Pd'] = - data['Pd']
        # data[name] = data[name].apply(lambda x: round(float(x.mid), 3)).astype(float)

        fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)

        bins = np.linspace(start=df[name].min(), stop=df[name].max(), num=20)
        for ax, metal, invert, color in zip(axes.ravel(), ['Pd', 'Ni'], [False, True], ['C0', 'C1']):
            metal_df = df[df['metal'] == metal]
            sns.histplot(data=metal_df, x=name, bins=bins, color=color, ax=ax, alpha=0.7, legend=True, label=metal)

        max_ytick = max([tick for ax in axes.ravel() for tick in ax.get_yticks()])
        min_ytick = min([tick for ax in axes.ravel() for tick in ax.get_yticks()])
        for ax in axes.ravel():
            ax.set_ylim(min_ytick, max_ytick)

        axes[1].invert_yaxis()
        plt.subplots_adjust(hspace=0)
        plt.legend(loc='lower left')

        plt.savefig(f'../stats/pop_{donor}.svg')
        plt.close()

    # Histograms
    for donor in cols:
        plt.figure()
        name = f'dist_{donor}' if donor not in others else donor
        sns.histplot(data=df, x=name, hue='metal')
        plt.xlabel(donor)
        plt.savefig(f'../stats/hist_{donor}.svg')
        plt.close()

    # 2D scatter plots
    for donor1, donor2 in [('P', 'N'), ('P', 'Br'), ('P', 'C'), ('N', 'Br'), ('N', 'C'), ('Br', 'C')]:
        for dist in ['dist', 'dist_norm']:
            plt.figure()
            sns.scatterplot(data=df, x=f'{dist}_{donor1}', y=f'{dist}_{donor2}', hue='metal', alpha=0.4)
            label = 'Distance' if dist == 'dist' else 'Normalized distance'
            plt.xlabel(f'{label} to {donor1} in Angstroms')
            plt.ylabel(f'{label} to {donor2} in Angstroms')
            plt.savefig(f'../stats/scatter_{dist}_{donor1}_{donor2}.svg')
            plt.close()

    # 2D scatter plots with Mulliken charge as hue
    for donor1, donor2 in [('P', 'N'), ('P', 'Br'), ('P', 'C'), ('N', 'Br'), ('N', 'C'), ('Br', 'C')]:
        for dist in ['dist', 'dist_norm']:
            plt.figure()
            palette = sns.color_palette("coolwarm", as_cmap=True)#'coolwarm'
            markers = ['o', '^']
            ax = sns.scatterplot(data=df, x=f'{dist}_{donor1}', y=f'{dist}_{donor2}', hue='Metal charge', style='metal', palette=palette, alpha=1, markers=markers)

            norm = plt.Normalize(df['Metal charge'].min(), df['Metal charge'].max())#-0.5, -0.2)#
            sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
            sm.set_array([])

            # Remove the legend and add a colorbar
            # ax.get_legend().remove()
            ax.figure.colorbar(sm)

            h, l = ax.get_legend_handles_labels()
            plt.legend(h[6:], l[6:])

            label = 'Distance' if dist == 'dist' else 'Normalized distance'
            plt.xlabel(f'{label} to {donor1} in Angstroms')
            plt.ylabel(f'{label} to {donor2} in Angstroms')
            plt.savefig(f'../stats/scatter_metal_charge_{dist}_{donor1}_{donor2}.svg')
            plt.close()

    # Data analytics
    df_stats = df.groupby('metal').agg(['max', 'min', 'mean', 'std']).T
    print(df_stats)
    df_corr = df.corr(method='spearman')
    print('Done!')



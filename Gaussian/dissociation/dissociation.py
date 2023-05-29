import os
import pymatgen.core as pt
import re
import pickle
import numpy as np
class Calculation:
    def __init__(self, calc_path, calc_name):
        self.calc_name = calc_name
        self.calc_path = calc_path
        self.file_path_list = self.get_files()
        ###Keeps track of how many times a calculation has been restarted
        self.restart_status = self.calculation_restarted()  # Dictionary That contains info on how many restarts
        self.restart_tag = self.get_continuation_stamp()    # i.e. _cont or _cont_cont or _cont_cont_cont
        self.com, self.log, self.chk, self.pkl, = self.get_latest_files()
        self.log_atoms, self.log_x, self.log_y, self.log_z = self.log_to_xyz()
        self.Au_C, self.Au_N = self.load_pkl()


    def get_files(self):
        file_list = os.listdir(self.calc_path)
        file_path_list = []
        for file in file_list:
            file_path_list.append(self.calc_path+"/"+file)
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
            elif (str(file_path).endswith(".pkl")) and (str(file_path).count("_cont") == self.restart_status[".pkl"]):
                tmp_pkl = file_path
            elif (str(file_path).endswith(".chk")) and (str(file_path).count("_cont") == self.restart_status[".chk"]):
                tmp_chk = file_path
            else:
                pass
        assert (tmp_com is not None) and (tmp_log is not None) and (tmp_pkl is not None) and (tmp_chk is not None)
        return tmp_com, tmp_log, tmp_chk, tmp_pkl


    def calculation_restarted(self):
        #This also serves as a check
        #to make sure all the correct
        #files are present
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
            elif str(file_path).endswith(".pkl"):
                pkl_detected.append(str(file_path).count("_cont"))
            elif str(file_path).endswith(".chk"):
                chk_detected.append(str(file_path).count("_cont"))
        tmp.update({".com": max(com_detected),
                    ".log": max(log_detected),
                    ".pkl": max(pkl_detected),
                    ".chk": max(chk_detected)})
        assert (len(com_detected) != 0) and (len(chk_detected) != 0) and (len(log_detected) != 0) and (len(pkl_detected) != 0)
        assert tmp[".com"] == tmp[".chk"]
        return tmp

    def get_continuation_stamp(self):
        i = 0
        stamp = "_cont"
        final_stamp = ""
        for i in range(self.restart_status[".com"] + 1):
            final_stamp = final_stamp + stamp
        return final_stamp



    def log_to_xyz(self):
        log_file = open(self.log, 'r')
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
                    ("                        Standard orientation:                        " in rline[i - 4]) and not already_got_coords and ((Standard_orientation_count_total - Standard_orientation_count == 1) or Standard_orientation_count_total == 1):
                already_got_coords = True
                detected_start_coord = True
            else:
                pass
            i = i + 1
        return symbol_list, x_coord_list, y_coord_list, z_coord_list

    def load_pkl(self):
        file = open(self.pkl, 'rb')
        ligand = pickle.load(file)
        elem = ligand.local_elements
        indexes = ligand.ligand_to_metal
        metal_coords = np.array([self.log_x[0], self.log_y[0], self.log_z[0]])
        carbon_atom_coords = np.array([self.log_x[indexes[elem.index('C')]+1], self.log_y[indexes[elem.index('C')]+1], self.log_z[indexes[elem.index('C')]+1]])
        nitrogen_atom_coords = np.array([self.log_x[indexes[elem.index('N')]+1], self.log_y[indexes[elem.index('N')]+1], self.log_z[indexes[elem.index('N')]+1]])
        Au_N_bond_lenght = np.linalg.norm(metal_coords-nitrogen_atom_coords)
        Au_C_bond_lenght = np.linalg.norm(metal_coords-carbon_atom_coords)

        print(f"[{str(self.log).split('/')[-1]}] Au-C --> [{Au_C_bond_lenght}]")
        print(f"[{str(self.log).split('/')[-1]}] Au-N --> [{Au_N_bond_lenght}]")
        return Au_N_bond_lenght, Au_C_bond_lenght


    def check_distances(self):
        if self.Au_C or self.Au_N > 3:
            print(f"calculation [{self.calc_name}] has failed Au_C has ")
        pass




"/home/users/clarkc18/CASINI_GOLD_CALCS/CASINI_230523_calcs"
if __name__ == "__main__":
    root = "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/CASINI_test_220523"
    list_of_directories = os.listdir(root)
    bond_lenght_dist = {}
    for calculation in list_of_directories:
        calculation_path = root+"/"+calculation
        calc = Calculation(calc_path=calculation_path, calc_name=calculation)
        bond_lenght_dist.update({calc.calc_name: [calc.Au_N, calc.Au_C]})
    print(bond_lenght_dist)
    pass

import os
import warnings
from pathlib import Path
import pymatgen.core as pt
import re
class Post_Process:
    def __init__(self,master_path, calc_prefix, file_name):
        self.master_path = master_path
        self.calc_prefix = calc_prefix
        self.list_of_all_directories = os.listdir(self.master_path)
        self.list_of_all_cals = self.get_cals()
        self.file_name = file_name
        self.current_working_directory = self.master_path + "/" + self.file_name
        self.com_file_path = self.current_working_directory + "/" + self.file_name + ".com"
        self.pkl_file_path = self.current_working_directory + "/" + self.file_name + ".pkl"
        self.log_file_path = self.current_working_directory + "/" + self.file_name + ".log"
        self.chk_file_path = self.current_working_directory + "/" + self.file_name + ".chk"

    def get_cals(self):
        tmp_calc_list = []
        for directory in self.list_of_all_directories:
            if not str(directory).startswith("AuCl2"):
                warnings.warn(f"!!!Warning!!! -> directory [{directory}] is not a calculation -> continuing program")
            else:
                tmp_calc_list.append(directory)
        return tmp_calc_list

    def check_all_files_present(self):
        ###We check to make sure all these files exist###
        path1_check = Path(self.com_file_path).is_file()
        path2_check = Path(self.pkl_file_path).is_file()
        path3_check = Path(self.log_file_path).is_file()
        path4_check = Path(self.chk_file_path).is_file()
        if path1_check and path2_check and path3_check and path4_check:
            return True
        else:
            print(f"!!!Fatal Error!!!-> There are files missing in the calculation [{self.current_working_directory}] -> Aborting Program")
            print(f".com = [{path1_check}]\n.pkl = [{path2_check}]\n.log = [{path3_check}]\n.chk = [{path4_check}]\n")
            raise ValueError

    def Optimisation_Completed(self):
        # Returns True if NBO completed successfully
        # Returns False if not
        log_file = open(self.log_file_path, 'r')
        for line in log_file.readlines():
            if "-- Stationary point found" in line:
                return True
            else:
                pass
        return False

    def NBO_Completed(self):
        # Returns True if NBO completed successfully
        # Returns False if not
        log_file = open(self.log_file_path, 'r')
        for line in log_file.readlines():
            if "NBO analysis completed" in line:
                return True
            else:
                pass
        return False



    def logtoxyz(self):
        log_file = open(self.log_file_path, 'r')
        start = None
        i = 0
        rline = log_file.readlines()
        detected_start_coord = False
        already_got_coords = False
        symbol_list = []
        x_coord_list = []
        y_coord_list = []
        z_coord_list = []
        for line in rline:

            if "---------------------------------------------------------------------" in line:
                detected_start_coord = False

            if detected_start_coord == True:
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

            if ("---------------------------------------------------------------------" in line) and \
                    ("Number     Number       Type             X           Y           Z" in rline[i - 1]) and \
                    ("Center     Atomic      Atomic             Coordinates (Angstroms)" in rline[i - 2]) and \
                    ("---------------------------------------------------------------------" in rline[i - 3]) and \
                    ("                        Standard orientation:                        " in rline[i - 4]) and not already_got_coords:
                already_got_coords = True
                detected_start_coord = True
            else:
                pass

            i = i + 1

        return symbol_list, x_coord_list, y_coord_list, z_coord_list


if __name__ == "__main__":
    root = "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/test"
    list_of_directories = os.listdir(root)
    for calculation in list_of_directories:
        # we assume that we only need to create cont once
        PP = Post_Process(master_path=root,
                          calc_prefix="AuCl2",
                          file_name=calculation)

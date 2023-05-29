import os
import pymatgen.core as pt
import re
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
        self.Au_C, self.Au_N = self.get_NBO_data()


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

            if start_found and not end_found:
                line_list = line.split()
                # NITROGEN
                # Here we make sure that LP comes before N which comes before Au which should come before Cl
                if (")Au" in line) and ("-Cl" in line) and ("LP" in line) and (len(line_list) > 3) and ("N" in line) and (line.find("LP") < line.find("N")) and (line.find("N") < line.find("Au")) and (
                        line.find("N") < line.find("Cl")):
                    interaction_energy_N = line_list[-3]
                    nitrogen_interaction_sum = interaction_energy_N
                    #print(line.strip() + f"         Interaction_N: [{interaction_energy_N}]")
                    result_N_found = True
                else:
                    pass

                # CARBON
                if (")Au" in line) and ("-Au" in line) and ("BD" in line) and (len(line_list) > 3) and ("C" in line) and ("Cl" not in line) and (line.find("BD") < line.find("C")) and (
                        line.find("C") < line.find("-Au")) and (line.find("-Au") < line.find(")Au")):
                    interaction_energy_C = float(line_list[-3])

                    if interaction_energy_C >= 5.0:
                        carbon_interaction_sum += interaction_energy_C
                    #print(line.strip() + f"         Interaction_C: [{interaction_energy_C}]")
                    result_C_found = True
                else:
                    pass

        if not result_N_found:
            self.Au_N = 1000
            print("!!!Fatal_Error!!!-> No value found for Au-N")
            #raise ValueError

        if not result_C_found:
            self.C = -1000
            print("!!!Fatal_Error!!!-> No value found for Au-C")
            #raise ValueError

        self.Au_N = nitrogen_interaction_sum
        self.Au_C = carbon_interaction_sum





"/home/users/clarkc18/CASINI_GOLD_CALCS/CASINI_230523_calcs"
if __name__ == "__main__":
    root = "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/CASINI_test_220523"
    list_of_directories = os.listdir(root)
    NBO = {}
    for calculation in list_of_directories:
        calculation_path = root+"/"+calculation
        calc = Calculation(calc_path=calculation_path, calc_name=calculation)
        NBO.update({calc.calc_name: [calc.Au_N, calc.Au_C]})
    print(NBO)
    pass





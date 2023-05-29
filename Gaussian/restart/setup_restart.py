# This function is to be called in the parent directory of the calculations
import os
import pymatgen.core.periodic_table as pt
import re
from copy import deepcopy
from pathlib import Path
import warnings


def Optimisation_Completed(log_path: str = None):
    # Returns True if NBO completed successfully
    # Returns False if not
    log_file = open(log_path, 'r')
    for line in log_file.readlines():
        if "-- Stationary point found" in line:
            return True
        else:
            pass
    return False


def NBO_Completed(log_path: str = None):
    # Returns True if NBO completed successfully
    # Returns False if not
    log_file = open(log_path, 'r')
    for line in log_file.readlines():
        if "NBO analysis completed" in line:
            return True
        else:
            pass
    return False


def logtoxyz(log_path):
    log_file = open(log_path, 'r')
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


def cont_com_file(old_com_path, cont_dir_path, sym_lst, x_lst, y_lst, z_lst):
    com_file = open(old_com_path, 'r')
    rline = com_file.readlines()
    com_list = []
    i = 0
    for line in rline:
        if ".chk" in line:
            line_list = line.split(".")
            line_list[-1] = "_cont." + deepcopy(line_list[-1])
            string = "".join(map(str, line_list))
            com_list.append(string)

        elif ("#p" in line) and ("guess=read" not in line):
            com_list.append(str(line).strip() + " guess=read" + "\n")

        elif len(line.split("  ")) == 4:
            atoms_stats = line.split("  ")
            atoms_stats[1] = "  " + str(deepcopy(x_lst[i]))
            atoms_stats[2] = "  " + str(deepcopy(y_lst[i]))
            atoms_stats[3] = "  " + str(deepcopy(z_lst[i])) + "\n"
            cont_com_xyz = "".join(map(str, atoms_stats))
            com_list.append(cont_com_xyz)
            i += 1


        else:
            com_list.append(line)

    old_com_file_name = old_com_path.split("/")[-1]
    name_list = old_com_file_name.split(".")
    name_list[-1] = "_cont." + deepcopy(name_list[-1])
    cont_com_name = "".join(map(str, name_list))
    os.system(f"touch {cont_dir_path + '/' + cont_com_name}")

    with open(cont_dir_path + '/' + cont_com_name, 'w') as f:
        for line in com_list:
            f.write(f"{line}")
    print("done")


def cont_chk_file(old_chk_path, cont_dir_path):
    old_chk_name = str(old_chk_path.split("/")[-1])
    old_chk_name_lst = old_chk_name.split(".")
    old_chk_name_lst[-1] = "_cont." + deepcopy(old_chk_name_lst[-1])
    new_chk_name = "".join(map(str, old_chk_name_lst))
    os.system(f"cp {old_chk_path} {cont_dir_path}/{new_chk_name}")
    pass


path = "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/test"

if __name__ == "__main__":
    list_of_directories = os.listdir(path)
    ###We make sure there are no extra directories that we are not aware of
    for calculation in list_of_directories:
        if not str(calculation).startswith("AuCl2"):
            print(f"!!!Fatal Error!!! -> directory [{calculation}] that should not be there, please fix this issue -> Aborting Program")
            raise ValueError
    else:
        pass

    for calculation in list_of_directories:
        input_file_path = path + "/" + str(calculation) + "/" + str(calculation) + ".com"
        info_file_path = path + "/" + str(calculation) + "/" + str(calculation) + ".pkl"
        output_file_path = path + "/" + str(calculation) + "/" + str(calculation) + ".log"
        chk_file_path = path + "/" + str(calculation) + "/" + str(calculation) + ".chk"

        ###We check to make sure all these files exist###
        path1_check = Path(input_file_path).is_file()
        path2_check = Path(info_file_path).is_file()
        path3_check = Path(output_file_path).is_file()
        path4_check = Path(chk_file_path).is_file()

        if path1_check and path2_check and path3_check and path4_check:
            Opt_status = Optimisation_Completed(output_file_path)
            if Opt_status:  # If Optimisation was successful
                NBO_status = NBO_Completed(output_file_path)
                if NBO_status:
                    print(f"NBO Completed Successfully for [{calculation}]")
                    pass
                else:
                    print(f"!!!Fatal Error!!! -> {calculation} NBO calculation has completely failed. This should not have happened -> Aborting Program")
                    raise NotImplementedError


            elif not Opt_status:  # If optimisation was not successful we prepare to restart the calculation
                warnings.warn(f"!!!Warning!!! -> Optimisation failed for calculation [{calculation}] -> Setting up restart")
                ###Path of continuation folder###
                continuation_directory_path = path + '/' + str(calculation) + '/' + 'cont'
                ###Make continuation folder###
                os.system(f"mkdir -p {continuation_directory_path}")
                ###Get xyz from log file###
                atom_symbol_list, x, y, z = logtoxyz(log_path=output_file_path)
                ###create restart com file with updated coordinates###
                cont_com_file(old_com_path=input_file_path,
                              cont_dir_path=continuation_directory_path,
                              sym_lst=atom_symbol_list,
                              x_lst=x, y_lst=y, z_lst=z)
                ###copy chk file to cont and change its name###
                cont_chk_file(old_chk_path=chk_file_path,  # todo this should work but cannot test it yet
                              cont_dir_path=continuation_directory_path)
        else:
            print(f"!!!Fatal Error!!!-> There are files missing in the calculation [{calculation}] -> Aborting Program")
            print(f".com = [{path1_check}]\n.pkl = [{path2_check}]\n.log = [{path3_check}]\n.chk = [{path4_check}]\n")
            raise ValueError

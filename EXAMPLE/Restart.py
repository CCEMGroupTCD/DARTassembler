###############################
# The Goal of this
# program is to provide a
# method to restart calculations
# that need to be restarted.
# This will involve the
# generation of a restart script
# and a new .com file
###############################
import pathlib
import pymatgen.core as pt
from copy import deepcopy
import re
import shutil
from src05_Assembly_Refactor.Submission import SumbitCalc


class RestartCalc:
    def __init__(self, global_batch_directory, previous_restarts: int = 0):

        self.cont_com_file_name = None
        assert global_batch_directory is not None
        self.global_path = pathlib.Path(global_batch_directory)
        assert self.global_path.is_dir()
        self.batch_dir = []
        self._get_batch_dir()
        assert self.batch_dir != []
        self.current_tag = str("" + ("_cont" * previous_restarts))
        self.restart_tag = str("" + ("_cont" * (previous_restarts + 1)))
        self.log_symbol = None
        self.log_x = None
        self.log_y = None
        self.log_z = None

    def _get_batch_dir(self):
        for child in self.global_path.iterdir():
            if child.is_file():
                pass
            elif child.is_dir():
                self.batch_dir.append(child)

    def iter_calcs(self):
        # Iterate Batches
        counter = 0
        for batch in self.batch_dir:
            # Iterate entries
            for child in batch.iterdir():
                if child.is_file():
                    pass
                if child.is_dir():
                    for entries in child.iterdir():
                        if entries.is_dir() and str(entries.name) == "complexes":
                            for calculation in entries.iterdir():
                                for file in calculation.iterdir():
                                    if file.is_dir():
                                        raise ValueError
                                    elif str(file.name).endswith(self.current_tag + ".log"):
                                        if self._check_progress(file):
                                            pass
                                        else:
                                            print(f"{str(file.name)}")
                                            # Extract updated structure from the .log and then this updates the relevant self.log_x ... and so on attributes
                                            self._get_final_structure(file)
                                            # We then generate the restart file
                                            self._gen_restart_com_file(file.parent)
                                            # Then we need to generate a new chk
                                            self._gen_new_chk(file)
                                            # Generate the new submission script
                                            sh_file_string = SumbitCalc(calc_name=self.cont_com_file_name.split("/")[-1].split(".com")[0], queue="amd", num_cores="32", requested_time="72:00:00", node_preference=[6, 7]).gen_submission_string()

                                            with open(str(file.parent)+"/"+"run_cont.sh", 'w') as f:
                                                for line in sh_file_string:
                                                    f.write(f"{line}")

                                            a=3
                                            pass

                        else:
                            pass


                else:
                    print(f"!!!Fatal_Error!!! --> Unexpected directory [{child.name}] found in [{child.parent}]")
                    raise ValueError

    @staticmethod
    def _check_progress(file: pathlib):
        # Here we check to see if the calculation successfully converged
        calc_completed = False
        with open(file, "r") as log_file:
            for line in log_file:
                if "-- Stationary point found" in str(line):
                    calc_completed = True
                    break
                else:
                    # Stationary point not found
                    pass
        return calc_completed

    def _get_final_structure(self, file):
        log_file = open(file, 'r')
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
            self.log_symbol = symbol_list
            self.log_x = x_coord_list
            self.log_y = y_coord_list
            self.log_z = z_coord_list

    def _gen_restart_com_file(self, dir: pathlib):
        # file must be a .com file
        assert dir.is_dir()
        com_file_path = None
        for file in dir.iterdir():
            if str(file.name).endswith(self.current_tag + ".com"):
                com_file_path = file
            else:
                pass

        assert com_file_path is not None
        com_file = open(com_file_path, 'r')
        rline = com_file.readlines()
        com_list = []  # This list contains all the lines for our new folder
        i = 0
        for line in rline:
            if ".chk" in line:
                line_list = line.split(".")
                line_list[-1] = self.restart_tag + "." + deepcopy(line_list[-1])
                string = "".join(map(str, line_list))
                com_list.append(string)

            elif ("#p" in line) and ("guess=read" not in line):
                com_list.append(str(line).strip() + " guess=read" + "\n")

            # Here we add the updated coordinates from the log file
            elif len(line.split("  ")) == 4:
                atoms_stats = line.split("  ")
                atoms_stats[1] = "  " + str(deepcopy(self.log_x[i]))
                atoms_stats[2] = "  " + str(deepcopy(self.log_y[i]))
                atoms_stats[3] = "  " + str(deepcopy(self.log_z[i])) + "\n"
                cont_com_xyz = "".join(map(str, atoms_stats))
                com_list.append(cont_com_xyz)
                i += 1
            else:
                com_list.append(line)
        restart_com_name = str(com_file_path.parent) + "/" + str(com_file_path.stem) + self.restart_tag + ".com"
        # Here we write the new com file
        self.cont_com_file_name = restart_com_name
        with open(restart_com_name, 'w') as f:
            for line in com_list:
                f.write(f"{line}")
        return restart_com_name

    def _gen_new_chk(self, log_file: pathlib):
        source_path = str(log_file).replace("gaussian.log", "gaussain.chk")
        destination_path = str(log_file).replace("gaussian.log", "gaussain_cont.chk")
        print(source_path)
        print(destination_path)
        shutil.copy(source_path, destination_path)


if __name__ == "__main__":
    restart_calc = RestartCalc(global_batch_directory="/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/EXAMPLE/DATA/data_before_restarts/DART_Example_Pd_Ni_Complexes", previous_restarts=0)
    restart_calc.iter_calcs()

# Calcs that need to be restarted
"""
UTAWATIM_PN_Pd_gaussian.log
XOLAVOKO_PN_Pd_gaussian.log
ADIXUDUR_PN_Pd_gaussian.log
LORUSEXO_PN_Pd_gaussian.log
NUREWEVE_PN_Pd_gaussian.log
ERIRECUF_PN_Pd_gaussian.log
UMOGALUJ_PN_Pd_gaussian.log
BOSASEBE_PN_Pd_gaussian.log
AKIMIWER_PN_Ni_gaussian.log
HUQAKOHU_PN_Ni_gaussian.log
ETEXIDOY_PN_Ni_gaussian.log
EMUBUBIC_PN_Ni_gaussian.log
AHUCIWOF_PN_Ni_gaussian.log
EQUXUXUD_PN_Ni_gaussian.log
FUGAZIHE_PN_Ni_gaussian.log
UVOMIVUN_PN_Ni_gaussian.log
WAQANOTU_PN_Ni_gaussian.log
OGAYOMUM_PN_Ni_gaussian.log
DINOPEWA_PN_Ni_gaussian.log
OPUSIMUK_PN_Ni_gaussian.log
URICAPIW_PN_Ni_gaussian.log
UHEQILUS_PN_Ni_gaussian.log
NEGACASI_PN_Ni_gaussian.log
ADITUMIG_PN_Ni_gaussian.log
AJEGOGAJ_PN_Ni_gaussian.log
OYAPALUQ_PN_Ni_gaussian.log
SOSOWIFE_PN_Ni_gaussian.log
NUPAKOTU_PN_Ni_gaussian.log
"""

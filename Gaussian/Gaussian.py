from pathlib import Path
import yaml
import datetime
import json
import re

now = datetime.datetime.now()  # Get the current date and time
date = now.date()  # Extract the date from the datetime object
time = now.time()  # Extract the time from the datetime object


class FormatGaussianInput:
    def __init__(self, yaml_path, output_path: str = None, instruction: str = None):

        #
        #
        # Global Attributes
        self.yaml_path = None
        self.directory_path = None
        self.instruction = instruction
        assert (self.instruction == "purge_input_files") or (self.instruction == "generate_input_files")
        if instruction != "purge_input_files":
            assert yaml_path is not None
        else:
            pass
        assert output_path is not None
        assert instruction is not None
        self.yaml_path = yaml_path
        self.directory_path = Path(output_path)

        #
        #
        # Current Calculation Attributes
        self.current_calculation_parent_directory = None
        self.current_json_path = None
        self.current_json_data = None
        self.filename = None
        self.xyz = None
        self.metal_type = None
        self.ligands = None
        self.num_processors = None
        self.memory_GB = None
        self.charge = None
        self.multiplicity = None
        self.calc_instruction = None
        self.functional = None
        self.basis_set_instruction = None
        self.pseudo_potential_instruction = None
        self.spacer_message = None
        self.basis_sets_dict = None
        self.basis_set_seperator = None

        #
        #
        # Open Yaml file
        if instruction != "purge_input_files":
            self._read_yaml()
        else:
            pass

        #
        #
        # Error Catching
        self._assert_input_directory()


    def _assert_input_directory(self):
        directories_found = 0
        for item in self.directory_path.iterdir():
            if item.is_dir():
                directories_found += 1
            else:
                file_extension = item.suffix
                assert (file_extension == ".csv") or (file_extension == ".xyz")
        assert directories_found == 1

    def _assert_calculation_directory(self):
        directories_found = 0
        for item in self.current_calculation_parent_directory.iterdir():
            if item.is_dir():
                directories_found += 1
            else:
                file_extension = item.suffix
                assert (file_extension == ".json") or (file_extension == ".xyz") or (file_extension == ".yml") or (file_extension == ".com")
        assert directories_found == 0

    def _read_json(self):
        with open(self.current_json_path, 'r') as f:
            # Load the JSON data from the file
            self.current_json_data = json.load(f)
        self.ligands = self.current_json_data["complex"]["ligand_props"]
        self.metal_type = self.current_json_data["complex"]["metal"]
        self.charge = self.current_json_data["complex"]["total_charge"]
        self.multiplicity = self.current_json_data["complex"]["ligand_props"]["0"]["Metal Spin:"]  # TODO: !!!WARNNG!!!WARNNG!!!WARNNG!!!WARNNG!!!WARNNG!!!WARNNG!!!WARNNG!!! This is unstable and needs to be urgently replaced
        self.filename = self.current_calculation_parent_directory.name
        self.xyz = self.current_json_data["xyz_structure"]

    def _read_yaml(self):
        with open(self.yaml_path, 'r') as file:
            self.config_data = yaml.safe_load(file)
        self.num_processors = self.config_data["num_processors"]
        self.memory_GB = self.config_data["memory_GB"]
        self.calc_instruction = self.config_data["calc_instruction"]
        self.functional = self.config_data["functional"]  # rwb97xd
        self.basis_set_instruction = self.config_data["basis_set_instruction"]  # gen This could also be, for example 6-31G(d)
        self.pseudo_potential_instruction = self.config_data["pseudo_potential_instruction"]  # pseudo/read
        self.spacer_message = str(self.config_data["spacer_message"]) + ". This file was generated on the " + str(date) + " at " + str(time).split(".")[
            0]  # f"This Gaussian Input files was generated using the exceptionally Brilliant DART program at {now.date()} /// {now.time()}"
        self.basis_sets_dict = self.config_data["basis_sets"]
        self.basis_set_seperator = "\n****\n"


    def purge_input_files(self):
        pass

    def _traverse_directory(self):
        for global_item in self.directory_path.iterdir():
            # we find the batches folder
            if global_item.is_dir():
                for batches in Path(global_item).iterdir():
                    if batches.is_dir():
                        for calculation_folder in Path(batches).iterdir():
                            if calculation_folder.is_dir():
                                for calculation in Path(calculation_folder).iterdir():
                                    self.current_calculation_parent_directory = calculation
                                    self._assert_calculation_directory()
                                    for calculation_data in Path(calculation).iterdir():
                                        if calculation_data.is_file():
                                            if calculation_data.suffix == ".json":
                                                # Now we need to decide what do depending on the instruction
                                                if self.instruction == "purge_input_files":
                                                    self.purge_input_files()

                                                elif self.instruction == "generate_input_files":
                                                    print(calculation_data)
                                                    self.current_json_path = calculation_data
                                                    self._read_json()
                                                    file_path = calculation / (str(self.filename) + ".com")
                                                    file_path.write_text(self.Generate_Gaussian_com())
                                                    pass
                                            else:
                                                pass
                                        else:
                                            pass
                            else:
                                pass
                    else:
                        pass
            else:
                pass

    def _gen_chk_name_(self):
        # This will ultimately be the first line in the input file
        # Its purpose is to specify the name of the checkpoint file
        line1 = f"%chk={self.filename}.chk\n"
        return str(line1)

    def _gen_num_proc(self):
        # Here we specify the number of processors for the calculation
        line2 = f"%nprocshared={str(self.num_processors)}\n"
        return str(line2)

    def _gen_mem(self):
        # This specifies the number of processors that are to be used in the calculation
        line3 = f"%mem={self.memory_GB}GB\n"
        return str(line3)

    def _gen_calc_type_and_theory(self):
        # Here the #p command specifies the type of calculation and the level of theory being used
        line4 = f"#p {self.calc_instruction} {self.functional}/{self.basis_set_instruction} {self.pseudo_potential_instruction}\n"
        return str(line4)

    def _gen_spacer_message(self):
        # Here specify a spacer message that is necessary in a Gaussian calculation
        line5 = f"\n{self.spacer_message}\n"
        return str(line5)

    def _gen_multiplicity_charge(self):
        # Here we specify the total charge and multiplity of the complex
        line6 = f"\n{self.charge} {self.multiplicity}"
        return str(line6)

    def _gen_atomic_coords(self):
        coordinates = self.xyz.split("\n\n")[1]
        return "\n" + str(coordinates) + "\n"

    def _gen_basi_sets(self):
        full_atom_str = ""
        for ligand in self.ligands.values():
            for character in ligand["stoichiometry"]:
                if character.isnumeric():
                    # C2H4O3 --> we want to skip the numbers
                    pass
                else:
                    # we append all the atoms to a string
                    full_atom_str = full_atom_str + character

            pass
        # Now we use regex to split the string up everytime it encounters a capital letter to yield a list of elements
        full_atom_list = re.split('(?<=.)(?=[A-Z])', full_atom_str)
        # Now we get the unique list of all elements
        # However this list will not contain our metal
        reduced_atom_list = list(set(full_atom_list))
        reduced_atom_list.insert(0, str(self.metal_type))
        basis_set_string = ""
        for atom in reduced_atom_list:

            try:
                basis_set_string_tmp = str(self.basis_sets_dict[str(atom)])
            except KeyError:
                basis_set_string_tmp = self.basis_sets_dict["other"]
                str(basis_set_string_tmp).replace("x", atom)
            basis_set_string = basis_set_string + basis_set_string_tmp + self.basis_set_seperator
        return basis_set_string

    def _gen_ecp(self):
        line8 = ""
        new_line = "\n"
        for atom_str, basis_set in self.basis_sets_dict.items():
            if basis_set.count("\n") == 3 and atom_str == self.metal_type:
                line8 = line8 + "\n" + f'{basis_set.split(new_line)[0]}{new_line}{basis_set.split(new_line)[1]}{new_line}'
            else:
                pass
        return line8

    @staticmethod
    def _gen_link1():
        line9 = "\n__Link1__\n"
        return line9

    def _gen_link1_header(self):
        line10 = f"#p Geom=AllCheck pseudo=read guess=read {self.functional}/{self.basis_set_instruction} pop=nbo7read\n\n"
        return line10

    @staticmethod
    def _gen_footer():
        line11 = "\n$nbo aonbo=c $end\n"
        return line11

    def Generate_Gaussian_com(self):
        file_string = self._gen_chk_name_() + self._gen_num_proc() + self._gen_mem() + self._gen_calc_type_and_theory() + self._gen_spacer_message() \
                      + self._gen_multiplicity_charge() + self._gen_atomic_coords() + self._gen_basi_sets() + self._gen_ecp() + self._gen_link1() \
                      + self._gen_link1_header() + self._gen_basi_sets() + self._gen_ecp() + self._gen_footer()

        print("\n\n\n")
        print(file_string)
        print("\n\n\n")
        return file_string


if __name__ == "__main__":
    test = FormatGaussianInput(yaml_path="/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/src05_Assembly_Refactor/Gaussian_config.yml",
                               output_path="/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/src14_Assembly_Unit_Test/output",
                               instruction="generate_input_files")

    test._traverse_directory()

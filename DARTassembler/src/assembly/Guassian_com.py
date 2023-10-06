from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand
import re
import datetime
import yaml
import pandas as pd

now = datetime.datetime.now()  # Get the current date and time
date = now.date()  # Extract the date from the datetime object
time = now.time()  # Extract the time from the datetime object


class Generate_Gaussian_input_file:
    def __init__(self, xyz: str = None,
                 complex_charge: int = None,
                 spin: int = None,
                 path_to_Gaussian_input_file: str = 'Gaussian_config.yml',
                 filename: str = None,
                 ligands: dict[RCA_Ligand] = None,
                 metal_type: str = None):
        #
        #
        # Open the yml file for the inputs
        with open(path_to_Gaussian_input_file, 'r') as file:
            config_data = yaml.safe_load(file)

        #
        #
        # Input parameters
        self.filename = filename
        self.xyz = xyz
        self.metal = metal_type
        self.ligands = ligands
        self.num_processors = config_data["num_processors"]
        self.memory_GB = config_data["memory_GB"]
        self.charge = int(complex_charge)
        self.multiplicity = spin
        self.calc_instruction = config_data["calc_instruction"]
        self.functional = config_data["functional"]  # rwb97xd
        self.basis_set_instruction = config_data["basis_set_instruction"]  # gen This could also be, for example 6-31G(d)
        self.pseudo_potential_instruction = config_data["pseudo_potential_instruction"]  # pseudo/read
        self.spacer_message = str(config_data["spacer_message"])  # f"This Gaussian Input files was generated using the exceptionally Brilliant DART program"
        self.basis_sets_dict = config_data["basis_sets"]
        self.basis_set_seperator = "\n****\n"
        self.all_elements = [self.metal] + list(pd.unique([el for lig in self.ligands.values() for el in lig.atomic_props['atoms']]))  # list of all elements in the complex

        #
        #
        # Generate Gaussian Input file
        self.Generate_Gaussian_com()

        #
        #
        # These assert statements check the validity of the input
        assert not str(self.filename).endswith(".com")
        assert isinstance(self.num_processors, int)

    def _gen_chk_name_(self):
        # This will ultimately be the first line in the input file
        # Its purpose is to specify the name of the checkpoint file
        line1 = f"%chk={self.filename}_gaussain.chk\n"
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
        basis_set_string = ""
        for atom in self.all_elements:
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
            if basis_set.count("\n") == 3 and atom_str == self.metal:
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
        return file_string

    def Generate_Gaussian_com_without_NBO(self):
        file_string = self._gen_chk_name_() + self._gen_num_proc() + self._gen_mem() + self._gen_calc_type_and_theory() + self._gen_spacer_message() \
                      + self._gen_multiplicity_charge() + self._gen_atomic_coords() + self._gen_basi_sets() + self._gen_ecp() + "\n" + "\n"
        return file_string

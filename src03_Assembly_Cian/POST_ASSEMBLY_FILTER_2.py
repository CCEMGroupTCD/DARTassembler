import stk
from openbabel import openbabel as ob
from openbabel import pybel
import ast
import os


class OPTIMISE:
    @staticmethod
    def movie(input_file, working_directory):
        # This function takes in an input of a
        movie_input_file = working_directory.format("movie_input.xyz")
        movie_ouput_file = working_directory.format("movie.xyz")
        input_file_ = working_directory.format(input_file)
        #os.system('obabel .mol {} .xyz -O  {}'.format(input_file_, movie_input_file))
        os.system(f"obabel .mol {input_file_} .xyz -O  {movie_input_file}")
        os.system('touch  {}'.format(movie_ouput_file))
        os.system('cat {} >> {}'.format(movie_input_file, movie_ouput_file))

    @staticmethod
    def reset_coords(file_read, file_write, indexes, working_directory):
        # file_read is the original_unoptimised_file and file_write is the file that is being optimised these are .mol files
        with open(file_read, "r") as f1, open(file_write, "r+") as f2:
            #print(indexes)
            i = 0
            data = []
            for line_read1, line_read2 in zip(f1.readlines(), f2.readlines()):
                if (indexes.count(i - 2) == 1) or (i == 4):                         # i == 4 corresponds to fixing the metal in place
                    data.append(line_read1)
                    #print(i)
                    #print(line_read1)
                else:
                    data.append(line_read2)
                i = i + 1
        temp_file = "tmp_uff_tmp.mol".format(working_directory)
        os.system('touch  {}'.format(temp_file))
        with open(temp_file, "w") as f3:
            for item in data:
                f3.write(item)
        f1.close()
        f2.close()
        f3.close()
        output_file = working_directory.format("tmp_uff.mol")
        os.system('mv {} {}'.format(temp_file, output_file))
        return None

    @staticmethod
    def get_energy(input_file):
        mol = next(pybel.readfile("mol", input_file))
        obmol = mol.OBMol
        ff = ob.OBForceField_FindType("uff")
        assert (ff.Setup(obmol))
        kj_to_kcal = 1.0 / 4.184
        ff.SetCoordinates(mol.OBMol)
        uffE = ff.Energy(False) * kj_to_kcal
        return uffE

    @staticmethod
    def optimiser(input_file, output_file, n, working_directory):
        #print("The input directory: " + str(input_file))
        #print("The output directory: " + str(output_file))
        mol = next(pybel.readfile("mol", str(input_file)))
        mol.localopt(forcefield='uff', steps=n)
        mol.write("mol", output_file, overwrite="True")

    @staticmethod
    def optimise_complex_step1(directory, func_group_lines, setting):
        if setting == "lock":
            i = 0
            for i in range(100):
                #print(i)
                if i == 0:
                    OPTIMISE.optimiser(directory.format("tmp_unoptimised.mol"), directory.format("tmp_uff.mol"), n=1, working_directory=directory)
                    OPTIMISE.reset_coords(directory.format("tmp_unoptimised.mol"), directory.format("tmp_uff.mol"), func_group_lines, directory)

                else:
                    OPTIMISE.optimiser(directory.format("tmp_uff.mol"), directory.format("tmp_uff.mol"), n=1, working_directory=directory)
                    OPTIMISE.reset_coords(directory.format("tmp_unoptimised.mol"), directory.format("tmp_uff.mol"), func_group_lines, directory)

                #OPTIMISE.movie(input_file="tmp_uff.mol", working_directory=directory)
            #OPTIMISE.optimiser(directory.format("tmp_uff.mol"), directory.format("tmp_uff_2.mol"), n=300, working_directory=directory)
            #OPTIMISE.movie(input_file="tmp_uff_2.mol", working_directory=directory)
            # Now we want to move tmp_uff_2.mol back to tmp_uff.mol
            directory_and_input_file = directory.format("tmp_uff.mol")
            directory_and_output_file = directory.format("tmp_optimisation_complete.mol")
            os.system("mv " + directory_and_input_file + " " + directory_and_output_file)
            del_file = directory.format("tmp_uff.mol")
            #os.system("rm {}".format(del_file))
        else:
            print("quickly all you need to do is fix up this bit of code with the correct directories")
            exit()
            OPTIMISE.optimiser("tmp.mol", "tmp_uff.mol", n=1000, working_directory=directory)

    @staticmethod
    def Optimise_STK_Constructed_Molecule(constructed_molecule, func_groups_str_input, func_groups_index_input, func_groups_type_input):
        working_directory = "../tmp/{}"
        constructed_molecule.write(working_directory.format("tmp_unoptimised_.mol"))
        os.system('obabel .mol {} .xyz -O  {}'.format(working_directory.format("tmp_unoptimised_.mol"), working_directory.format("tmp_unoptimised.mol")))
        os.system("rm {}".format(working_directory.format("tmp_unoptimised_.mol")))
        lines_to_be_saved = OPTIMISE.process_functional_groups(func_groups_index_input, func_groups_str_input, working_directory.format("tmp_unoptimised_.mol"), working_directory.format("tmp_unoptimised_.xyz"))
        OPTIMISE.optimise_complex_step1(working_directory, lines_to_be_saved, "lock")
        return stk.BuildingBlock.init_from_file(working_directory.format("tmp_optimisation_complete.mol"))

    @staticmethod
    def process_functional_groups(func_indexes, func_str, mol_file, xyz_file):
        os.system('obabel .mol {} .xyz -O  {}'.format(mol_file, xyz_file))
        # This function should take some input and return the lines in the line numbers file that need to be kept
        total_num_atoms_line = 1
        space = 1
        metal = 1  # assuming mono-metallic
        total_atoms = []  # each entry corresponds to a ligand
        for key in func_str.keys():
            #print(func_str[key])
            atom_counter = 0
            for line in str(func_str[key]).split("\n"):
                if str(line) != "":
                    if str(line)[0].isalpha():
                        atom_counter += 1
                    else:
                        pass
                else:
                    pass
            total_atoms.append(atom_counter)
        ligand_num = 0
        all_line_numbers = []
        for key in func_indexes.keys():
            index_list = ast.literal_eval(str(func_indexes[key]))
            #print(index_list)
            for selected_index in index_list:
                line_number = total_num_atoms_line + space + metal + int(selected_index)
                if ligand_num != 0:
                    i = 1
                    while (ligand_num - i) >= 0:
                        line_number = line_number + total_atoms[ligand_num - i]
                        i = i + 1
                else:
                    pass
                all_line_numbers.append(line_number)
            ligand_num += 1

        #print("This is all the lines to be saved")
        #print(all_line_numbers)
        return all_line_numbers  # THIS IS SUCCESSFULL!!!!

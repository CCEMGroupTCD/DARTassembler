import os


import stk
from src01.Molecule import RCA_Molecule
from src01.DataBase import LigandDB
from src01.utilities_graph import view_graph

from src03_Assembly_Cian.RandomComplexAssembler import RandomComplexAssembler
from src03_Assembly_Cian.TransitionMetalComplex import TransitionMetalComplex as TMC
import ast
import time
from openbabel import openbabel as ob
from timeit import default_timer as timer
from random import choice
from string import ascii_uppercase
import random

import matplotlib
from rdkit.Chem import rdmolfiles, Draw
from rdkit import Chem
from rdkit.Chem import AllChem

from openbabel import pybel
from rdkit.Chem import PyMol

matplotlib.use('TkAgg')


test_batch_list = [
    {'Name': 'first test batch', "Input_Path": "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/data/Filtered_Jsons/filteredLigDB_Cian_Ct_ff_graph.json",
     'Output_Path': '/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/', 'MAX_num_complexes': '123', 'Topology_1': '[2, 1, 1, ["1", "2", "2"]]',
     'Topology_2': '[ 2, 2, [ "1", "1" ] ]',
     'Metal_1': ['Cr', '+4', 'Low'], 'Metal_2': ['Cu', '+7', 'High'], 'Metal_3': ['V', '+2', 'Low'], "Isomers": "Generate All"},
    {'Name': 'second test batch', 'Output_Path': '/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/', 'MAX_num_complexes': '4', 'Topology_1': '[2, 1, 0]',
     'Topology_2': '[2, 1, 0]', 'Topology_3': '[2, 1, 0]',
     'Topology_4': '[2, 1, 1, ["1", "2", "2"]]', 'Metal_1': ['V', '+2', 'High'], 'Metal_2': ['Mn', '+2', 'High'], 'Metal_3': ['Mn', '+2', 'High'], "Isomers": "Generate Lowest Energy"},
    {'Name': 'third test batch', 'Output_Path': '/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/', 'MAX_num_complexes': '26395625',
     'Topology_1': '[4, 1, 1, ["1", "2", "2"]]', 'Metal_1': ['Cr', '+1', 'Low'],
     'Metal_2': ['Fe', '+1', 'Low'], "Isomers": "Generate Lowest Energy"}]


test_batch_list_2 = [{'Name': 'first test batch', "Input_Path": "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/data/Filtered_Jsons/filteredLigDB_Cian_Ct_ff_graph.json",
                      'Output_Path': '/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/', 'MAX_num_complexes': '100', 'Topology_1': '[ 2, 2, [ "1", "1" ] ]',
                      'Metal_1': ['Fe', '+2', 'High'], "Isomers": "Generate All", "Optimisation_Choice": "False", "Random_Seed": "10"}]


def movie(input_file, working_directory):
    # This function takes in an input of a
    movie_input_file = str(working_directory).format("movie_input.xyz")
    movie_ouput_file = str(working_directory).format("movie.xyz")
    input_file_ = working_directory.format(input_file)
    # os.system('obabel .mol {} .xyz -O  {}'.format(input_file_, movie_input_file))
    os.system(f"obabel .mol {input_file_} .xyz -O {movie_input_file}")
    os.system('touch  {}'.format(movie_ouput_file))
    os.system('cat {} >> {}'.format(movie_input_file, movie_ouput_file))


class Optimise:
    def __init__(self, complex_, ligands_, metal_, bb, metal_ox):
        self.complex = complex_
        self.ligands = ligands_  # todo need to switch from using ligands to using the correctly rotated stk_building blocks
        self.metal = metal_
        self.stk_building_blocks = bb
        self.metal_ox =metal_ox

    def optimiser_experimental(self):
        ligand_total_num_atoms = []  # Contains the total number of atoms for each ligand
        list_of_ligand_mol_objects = []  # List of rdkit mol representations of each ligand
        for i in range(len(self.ligands.keys())):  # Loop through all ligands
            RCA_Molecule_instance = RCA_Molecule(atomic_props=self.ligands[i].atomic_props)  # rdkit instance of the ligand
            ligand_total_num_atoms.append(len(self.ligands[i].atomic_props["atoms"]))
            #####HERE I MAKE A VERY DANGEROUS ASSUMPTION THAT THE ORDER OF THE LIGANDS IN THE LIGAND FILE IS THE SAME AS THE ORDER OF THE LIGANDS IN THE BB LOST######
            tmp_lig_mol = self.stk_building_blocks[0][i].to_rdkit_mol()
            tmp_xyz_string = rdmolfiles.MolToXYZBlock(tmp_lig_mol)
            conv = ob.OBConversion()  # Openbabel functionality to convert xyz to mol
            conv.SetInAndOutFormats("xyz", "mol")  # Set conversion settings
            mol_ = ob.OBMol()  # Initiate a rdkit mol object
            conv.ReadString(mol_, tmp_xyz_string)  # Read in xyz string
            string_ouput = conv.WriteString(mol_)  # Output mol string
            stk_mol = rdmolfiles.MolFromMolBlock(string_ouput, removeHs=False, sanitize=False, strictParsing=False)  # Generate mol object from mol string with all sanitization turned completely off
            stk_mol_without_Hg = Chem.RWMol()  # initialise mol object
            stk_mol_without_Hg.InsertMol(stk_mol)
            idx = []
            for atom in stk_mol_without_Hg.GetAtoms():
                if atom.GetAtomicNum() == 80:
                    print("here")
                    idx.append(int(atom.GetIdx()))
                else:
                    print("passing")
                    pass
            for item in idx:
                stk_mol_without_Hg.RemoveAtom(item)
            assert len(self.ligands[i].atomic_props["atoms"]) == stk_mol_without_Hg.GetNumAtoms()
            Draw.MolToFile(stk_mol_without_Hg, "image.png", size=(600, 600))
            list_of_ligand_mol_objects.append(stk_mol_without_Hg)  # append the list with the mol object

        combo = Chem.RWMol()  # Initiate editable rdkit mol object
        for ligand_mol_object in list_of_ligand_mol_objects:  # Combine all ligand mol objects
            combo = Chem.CombineMols(combo, ligand_mol_object)
        metal_xyz = f'1\n \n{self.metal} 0 0 0 '  # xyz string of just metal
        mol_b = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("xyz", "mol")
        conv.ReadString(mol_b, metal_xyz)
        metal_string_output = conv.WriteString(mol_b)
        mol_metal = rdmolfiles.MolFromMolBlock(metal_string_output, removeHs=False, sanitize=False, strictParsing=False)  # Created rdkit object of just metal atom
        combo = Chem.CombineMols(mol_metal, combo)  # Add metal with the rest of the ligands
        assembled_complex_mol = Chem.RWMol()  # create a new editable mol object
        assembled_complex_mol.InsertMol(combo)  # add comboned ligands and metal mols
        # Now we need to iterate through each coordinating atom and add a dative bond to the metal
        previous_num_atoms = 0  # This variable is updated with the sum of all atoms from all previous ligands
        for i in range(len(self.ligands.keys())):  # we iterate through all the ligands in the complex
            indexes = self.ligands[i].ligand_to_metal  # These are the indexes at which each coordinating atom occurs
            for index in indexes:
                assembled_complex_mol.AddBond(0, int(1 + int(index) + int(previous_num_atoms)), Chem.BondType.DATIVE)  # we add a bond between the metal and the coordinating groups
            previous_num_atoms += ligand_total_num_atoms[i]
        #
        #
        # optimisation
        loops = 100
        assembled_complex_mol.GetAtomWithIdx(0).SetFormalCharge(int(self.metal_ox))
        Chem.SanitizeMol(assembled_complex_mol, Chem.SANITIZE_ALL)
        for i in range(loops):
            AllChem.UFFOptimizeMolecule(assembled_complex_mol, 1)
            rdmolfiles.MolToXYZFile(assembled_complex_mol, "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/movie_frame.xyz")
            os.system(
                'cat {} >> {}'.format("/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/movie_frame.xyz", "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/movie.xyz"))

            pass
        Draw.MolToFile(assembled_complex_mol, "image.png", size=(600, 600))
        print("done")


class Assembly:

    def __init__(self, list_of_batches):
        self.progress_bar = 0  # This may be used in the future to return a progress bar to the gui
        self.list_of_batches = list_of_batches  # This is the main input from the gui to the assembly process
        self.total_assembled_complexes = 0  # This counter keeps track of the number of assembled complexes (including isomers)
        self.current_assembled_complex_name = None  # What's the name of the complex that is currently being assembled
        self.current_working_Directory = None  # This specifically handles the outputs

    @staticmethod
    def input_controller(batch_input):
        metal_list = []  # This will eventually contain all our selected metals
        topology_list = []  # This will eventually contain all our Topologies
        for key in batch_input.keys():
            if key.startswith("Metal"):
                metal_list.append(batch_input[key])
            elif key.startswith("Topology"):
                topology_list.append(ast.literal_eval(batch_input[key]))
            else:
                pass
        return metal_list, topology_list

    def output_controller(self, list_of_complexes_wih_isomers, ligands, metal, metal_ox_state, output_path, metal_multiplicity, name_of_batch, complex_time_start, topology, rot_bb, random_seed, input_batch):

        write_to_file = True  # Do we want to make an output
        isomer_counter = 0  # Counts the number of isomers that have been constructed
        for complex_ in list_of_complexes_wih_isomers:  # We loop through all the created isomers
            if (complex_ is not None) and write_to_file:
                Assembled_complex = TMC(compl=complex_, ligands=ligands, metal=metal, metal_charge=int(metal_ox_state))

                #Optimiser_instance.optimiser_experimental()
                # From here we need to optimise if the user asks, which will involve using the above graph to generate a .mo file and from the .mol file use stks built in optimiser uff
                # view_graph(assembled_complex_graph)
                #print("we have passed the graph merger")
                A = Assembled_complex.functional_groups
                #print("These_are_the_functional_group types" + str(A))
                if isomer_counter == 0:  # If this is the first isomer of a complex then ...

                    # todo: Outcommented, maybe we need another way to catch that error
                    # if str(output_path)[-1] != "/":
                    #     print("!!!Fatal_Error!!! -> your path should end in a '/' character -> Exiting Program")
                    #     exit()
                    # else:

                    self.current_working_Directory = "RCA_" + str(''.join(choice(ascii_uppercase) for j in range(12)))  # Generate random name of complex
                    os.system("mkdir " + str(output_path) + "{}".format(self.current_working_Directory))  # Set-up directory that will store all isomers/info on a particular complex
                else:
                    pass
                timer_end = timer()

                #Assembled_complex.mol.view_3d()  # This is purely for debugging

                self.total_assembled_complexes += 1  # We increment our complex counter (its placement here is arbitrary)
                isomer_counter += 1  # counts just the total number of isomers
                if write_to_file:
                    Assembled_complex.mol.print_to_xyz(str(output_path) + "tmp_output_test_cian.xyz")  # Print to a temporary file
                    Assembled_complex.mol.print_to_xyz(
                        str(str(output_path)) + "{}/{}{}".format(self.current_working_Directory, str(self.current_working_Directory) + "_isomer_" + str(isomer_counter + 1), ".xyz"))
                    os.system('touch  ' + str(output_path) + 'all.xyz')
                    os.system('cat ' + str(output_path) + 'tmp_output_test_cian.xyz >> {}/all.xyz'.format(str(output_path)))  # all.xyz contains a concatenated list of XYZs for ASE
                    # print("These are the functional groups" + str(Assembled_complex.functional_groups()))

                    # This is where you can write anything you want in the file
                    info = {"Batch_Name": name_of_batch,  # This dictionary is written into the info file line by line
                            "Time_of_Creation": str(time.strftime("%H:%M:%S", time.localtime())),
                            "Time_to_Create_all_isomers": timer_end - complex_time_start,
                            "Random_Seed": random_seed,
                            "Metal_Type": str(metal),
                            "Metal_Oxidation_State": str(metal_ox_state),
                            "Multiplicity": metal_multiplicity,
                            "Number_of_isomers": len(list_of_complexes_wih_isomers),
                            "Topology": str(topology),
                            }

                    if isomer_counter == len(list_of_complexes_wih_isomers):  # If this is the last isomer of this complex the ...
                        for keys in Assembled_complex.ligand_props.keys():
                            info.update({"Ligand_{}_Name".format(keys + 1): ligands[keys].name})
                            # todo: No original metal anymore
                            # info.update({"Ligand_{}_Original_Metal".format(keys + 1): ligands[keys].original_metal})
                            info.update({"Ligand_{}_Graph_Hash".format(keys + 1): ligands[keys].graph_hash})
                            #info.update({"Ligand_{}_Original_Metal_Symbol".format(keys + 1): ligands[keys].original_metal_symbol})
                            info.update({"Ligand_{}_Global_Props".format(keys + 1): ligands[keys].global_props})
                            info.update({"Ligand_{}_Atomic_Props".format(keys + 1): ligands[keys].coordinates})
                            info.update({"Ligand_{}_Denticity".format(keys + 1): ligands[keys].denticity})
                            info.update({"Ligand_{}_Coordinating_Atoms".format(keys + 1): ligands[keys].local_elements})
                            info.update({"Ligand_{}_Stoichiometry".format(keys + 1): ligands[keys].stoichiometry})

                        os.system('touch ' + str(output_path) + '{}/info.txt'.format(self.current_working_Directory))  # Create info.txt file
                        with open(str(output_path) + '{}/info.txt'.format(self.current_working_Directory), 'w') as f:
                            for key, value in info.items():
                                f.write('%s:%s\n' % (key, value))
                else:
                    pass
            else:
                pass

    def assembly_main(self):
        # This is the entry point to the assembly process

        # only for Felix, can be deleted in the future
        total_assembled_complexes = []

        for batch in self.list_of_batches:  # We loop through all the batches that have been received
            F = LigandDB.from_json(json_=batch["Input_Path"], type_="Ligand")  # We initiate the database in the RandomComplexAssembler
            RCA = RandomComplexAssembler(database=F)
            self.total_assembled_complexes = 0  # Initialise as 0 for each batch
            Metals, Topologies = self.input_controller(batch)  # Use the input controller to extract a list of selected metals and topologies
            MAX = int(batch["MAX_num_complexes"])  # MAX complexes we want to make (if "Generate all isomers" is selected we make slightly more)
            isomer_instruction = batch["Isomers"]  # To select whether to generate all isomers or not
            creation_path = batch["Output_Path"]  # Path to the output
            batch_name = batch["Name"]
            opt_choice = batch["Optimisation_Choice"]
            i = 0
            while self.total_assembled_complexes < MAX:  # We loop until the target number of complexes has been surpassed
                random.seed(int(batch["Random_Seed"]) + i)

                print("Assembling Complex No." + str(int(batch["Random_Seed"]) + i))
                self.progress_bar = round((self.total_assembled_complexes / MAX) * 100)  # May or may not use in the future
                complex_timer_start = timer()
                list_of_complexes_wih_isomers, ligands, metal, metal_ox_state, multiplicity, chosen_topology, rotated_building_blocks = RCA.create_random_TMC(Metals, Topologies, isomer_instruction,
                                                                                                                                                              opt_choice)  # !!!ASSEMBLY!!!
                # The following If statement is for error handling
                if list_of_complexes_wih_isomers is not None:
                    self.output_controller(list_of_complexes_wih_isomers, ligands, metal, metal_ox_state, creation_path, multiplicity, batch_name, complex_timer_start, chosen_topology,
                                           rotated_building_blocks, (str(int(batch["Random_Seed"]) + i)), batch)

                    for complex_ in list_of_complexes_wih_isomers:
                        Assembled_complex = TMC(compl=complex_,
                                                ligands=ligands,
                                                metal=metal,
                                                metal_charge=int(metal_ox_state)
                                                )
                        total_assembled_complexes.append(Assembled_complex)

                elif list_of_complexes_wih_isomers is None:
                    print("!!!Warning!!! -> None type complex list encountered -> Skipping to next complex")
                    pass
                i += 1
        print("Assembly Completed Successfully")
        return total_assembled_complexes
        #exit()  # Once the assembly is complete we exit


# Uncomment the following lines of code if you wish to do Debugging without the GUI
#instance = Assembly(test_batch_list_2)
#instance.assembly_main()

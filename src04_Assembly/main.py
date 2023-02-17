import os
from src01.DataBase import LigandDB
from src04_Assembly.RandomComplexAssembler import RandomComplexAssembler
from src04_Assembly.TransitionMetalComplex import TransitionMetalComplex as TMC
import ast
import time
from timeit import default_timer as timer
from random import choice
from string import ascii_uppercase
import random

test_batch_list = [
    {'Name': 'first test batch', "Input_Path": "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/data/Filtered_Jsons/filteredLigDB_OER_170223.json",
     'Output_Path': '/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/', 'MAX_num_complexes': '50',
     'Topology_1': '[4, 1, 0]', 'Topology_2': '[3, 2, 0]',
     'Metal_1': ['Fe', '+3', 'Low'], 'Metal_2': ['Co', '+3', 'High'], 'Metal_3': ['Mn', '+6', 'Low'], "Isomers": "Generate Lowest Energy", "Optimisation_Choice": "True", "Random_Seed": "0"}]

test_batch_list_2 = [{'Name': 'first test batch', "Input_Path": "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/data/Filtered_Jsons/filteredLigDB_Cian_Ct_ff_graph.json",
                      'Output_Path': '/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/', 'MAX_num_complexes': '100', 'Topology_1': '[ 2, 2, [ "1", "1" ] ]',
                      'Metal_1': ['Fe', '+2', 'High'], "Isomers": "Generate Lowest Energy", "Optimisation_Choice": "True", "Random_Seed": "10"}]


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

    def output_controller(self, list_of_complexes_wih_isomers, ligands, metal, metal_ox_state, output_path, metal_multiplicity, name_of_batch, complex_time_start, topology, random_seed, input_batch):

        write_to_file = True  # Do we want to make an output
        isomer_counter = 0  # Counts the number of isomers that have been constructed
        for complex_ in list_of_complexes_wih_isomers:  # We loop through all the created isomers
            if (complex_ is not None) and write_to_file:
                Assembled_complex = TMC(compl=complex_, ligands=ligands, metal=metal, metal_charge=int(metal_ox_state), spin=int(metal_multiplicity))
                # Assembled_complex.mol.view_3d()  # This is purely for debugging
                print("charge = "+str(Assembled_complex.total_charge))
                if Assembled_complex.total_charge == 0:
                    Assembled_complex.mol.print_to_xyz(str(output_path) + "tmp_output_test_cian.xyz")  # Print to a temporary file
                    os.system('touch  ' + str(output_path) + 'all.xyz')
                    os.system('cat ' + str(output_path) + 'tmp_output_test_cian.xyz >> {}/all.xyz'.format(str(output_path)))  # all.xyz contains a concatenated list of XYZs for ASE
                    print("done")
                    self.total_assembled_complexes += 1
                else:
                    print("Wrong Charge")
                    pass
                """
                # Optimiser_instance.optimiser_experimental()
                # From here we need to optimise if the user asks, which will involve using the above graph to generate a .mol file and from the .mol file use stks built in optimiser uff
                # view_graph(assembled_complex_graph)
                # print("we have passed the graph merger")
                A = Assembled_complex.functional_groups
                #print("These_are_the_functional_group types" + str(A))
                if isomer_counter == 0:  # If this is the first isomer of a complex then ...
                    if str(output_path)[-1] != "/":
                        print("!!!Fatal_Error!!! -> your path should end in a '/' character -> Exiting Program")
                        exit()
                    else:

                        self.current_working_Directory = "RCA_" + str(''.join(choice(ascii_uppercase) for j in range(12)))  # Generate random name of complex
                        os.system("mkdir " + str(output_path) + "{}".format(self.current_working_Directory))  # Set-up directory that will store all isomers/info on a particular complex
                else:
                    pass
                timer_end = timer()

                #Assembled_complex.mol.view_3d()  # This is purely for debugging

                #next = input("Press return")

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
                            info.update({"Ligand_{}_Original_Metal".format(keys + 1): ligands[keys].original_metal})
                            info.update({"Ligand_{}_Graph_Hash".format(keys + 1): ligands[keys].graph_hash})
                            info.update({"Ligand_{}_Original_Metal_Symbol".format(keys + 1): ligands[keys].original_metal_symbol})
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
                """
            else:
                print("!!!Warning!!! -> None type complex detected in list -> skipping to next complex")
                pass
        print("done")

    def assembly_main(self):
        # This is the entry point to the assembly process

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
                print("########################################################################################################################")
                random.seed(int(batch["Random_Seed"]) + i)
                print("Assembling Complex No." + str(int(batch["Random_Seed"]) + i))
                self.progress_bar = round((self.total_assembled_complexes / MAX) * 100)  # May or may not use in the future
                complex_timer_start = timer()
                list_of_complexes_wih_isomers, ligands, metal, metal_ox_state, multiplicity, chosen_topology = RCA.create_random_TMC(Metals, Topologies, isomer_instruction,
                                                                                                                                     opt_choice)  # !!!ASSEMBLY!!!
                # The following If statement is for error handling
                if list_of_complexes_wih_isomers is not None:
                    self.output_controller(list_of_complexes_wih_isomers, ligands, metal, metal_ox_state, creation_path, multiplicity, batch_name, complex_timer_start, chosen_topology,
                                           (str(int(batch["Random_Seed"]) + i)), batch)
                elif list_of_complexes_wih_isomers is None:
                    print("!!!Warning!!! -> None type complex list encountered -> Skipping to next complex")
                    pass
                i += 1
        print("Assembly Completed Successfully")
        exit()  # Once the assembly is complete we exit


# Uncomment the following lines of code if you wish to do Debugging without the GUI
# instance = Assembly(test_batch_list)
# instance.assembly_main()

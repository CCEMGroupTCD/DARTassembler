import os
from src01.DataBase import LigandDB
from src04_Assembly.RandomComplexAssembler import RandomComplexAssembler
from src04_Assembly.TransitionMetalComplex import TransitionMetalComplex as TMC
import ast
from timeit import default_timer as timer
import random
import matplotlib
from src04_Assembly.test_batches import test_batch_list, test_batch_list_2

matplotlib.use('TkAgg')


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

    def output_controller(self,
                          list_of_complexes_wih_isomers,
                          ligands,
                          metal,
                          metal_ox_state,
                          output_path,
                          metal_multiplicity
                          ):

        write_to_file = True  # Do we want to make an output

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

            else:
                print("!!!Warning!!! -> None type complex detected in list -> skipping to next complex")
                pass
        print("done")

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
                print("########################################################################################################################")
                random.seed(int(batch["Random_Seed"]) + i)
                print("Assembling Complex No." + str(int(batch["Random_Seed"]) + i))
                self.progress_bar = round((self.total_assembled_complexes / MAX) * 100)  # May or may not use in the future
                complex_timer_start = timer()
                list_of_complexes_wih_isomers, ligands, metal, metal_ox_state, multiplicity, chosen_topology = RCA.create_random_TMC(Metals, Topologies, isomer_instruction,
                                                                                                                                     opt_choice)  # !!!ASSEMBLY!!!
                # The following If statement is for error handling
                if list_of_complexes_wih_isomers is not None:
                    self.output_controller(list_of_complexes_wih_isomers,
                                           ligands,
                                           metal,
                                           metal_ox_state,
                                           creation_path,
                                           multiplicity)

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


# Uncomment the following lines of code if you wish to do Debugging without the GUI
# instance = Assembly(test_batch_list)
# instance.assembly_main()


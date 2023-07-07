from src01.DataBase import LigandDB
from src05_Assembly_Refactor.Assemble import PlacementRotation
from src05_Assembly_Refactor.ligands import ChooseRandomLigands, ChooseIterativeLigands
from src05_Assembly_Refactor.Isomer import BuildIsomers
from src05_Assembly_Refactor.Optimise import OPTIMISE
from src05_Assembly_Refactor.Post_Filter import PostFilter
import random
import stk, os
import warnings
from constants.Paths import project_path

warnings.simplefilter("always")


def visualize(input_complex):
    # This method allows mw to visualize in a blocking way during debug but is not essential at all
    print("initializing visualization")
    stk.MolWriter().write(input_complex, 'input_complex.mol')
    os.system('obabel .mol input_complex.mol .xyz -O  output_complex.xyz')
    os.system("ase gui output_complex.xyz")
    os.system("rm -f input_complex.mol")
    os.system("rm -f output_complex.xyz")
    print("visualization complete")


# These are our input dictionaries. they contain all the input information that the program requires to run
Integration_test = [
    {'Name': 'Integration_test_1',
     "Input_Path": project_path().extend("data", "Filtered_Jsons", "INTEGRATION_TEST_LIGAND_DATABASE_270623.json"),
     'Output_Path': project_path().extend("src14_Assembly_Unit_Test"),
     'MAX_num_complexes': '1',
     'Topology_1': '[2, 2]--[1, 2]',
     'Topology_2': '[2, 2]--[1, 1]',
     'Metal_1': ['Fe', '+2', '0'],
     "Isomers": "Generate Lowest Energy",
     "Optimisation_Choice": "False",
     "Total_Charge": 0,
     "Random_Seed": "0"},

    {'Name': 'Integration_test_2',
     "Input_Path": project_path().extend("data", "Filtered_Jsons", "INTEGRATION_TEST_LIGAND_DATABASE_270623.json"),
     'Output_Path': project_path().extend("src14_Assembly_Unit_Test"),
     'MAX_num_complexes': '1',
     'Topology_1': '[2, 1, 1]--[1, 2, 3]',
     'Metal_2': ['Fe', '+2', '0'],
     "Isomers": "Generate Lowest Energy",
     "Optimisation_Choice": "False",
     "Total_Charge": 0,
     "Random_Seed": "0"},

    {'Name': 'Integration_test_3',
     "Input_Path": project_path().extend("data", "Filtered_Jsons", "INTEGRATION_TEST_LIGAND_DATABASE_270623.json"),
     'Output_Path': project_path().extend("src14_Assembly_Unit_Test"),
     'MAX_num_complexes': '10',
     'Topology_1': '[3, 2, 1]--[1, 2, 3]',
     'Metal_2': ['Fe', '+2', '0'],
     "Isomers": "Generate All",
     "Optimisation_Choice": "False",
     "Total_Charge": 0,
     "Random_Seed": "0"},

    {'Name': 'Integration_test_4',
     "Input_Path": project_path().extend("data", "Filtered_Jsons", "INTEGRATION_TEST_LIGAND_DATABASE_270623.json"),
     'Output_Path': project_path().extend("src14_Assembly_Unit_Test"),
     'MAX_num_complexes': '5',
     'Topology_1': '[4, 1, 1]--[1, 2, 3]',
     'Metal_2': ['Fe', '+4', '0'],
     "Isomers": "Generate All",
     "Optimisation_Choice": "True",
     "Total_Charge": 0,
     "Random_Seed": "0"},

    {'Name': 'Integration_test_5',
     "Input_Path": project_path().extend("data", "Filtered_Jsons", "INTEGRATION_TEST_LIGAND_DATABASE_270623.json"),
     'Output_Path': project_path().extend("src14_Assembly_Unit_Test"),
     'MAX_num_complexes': '2',
     'Topology_1': '[5, 1]--[1, 2]',
     'Metal_2': ['Fe', '+2', '0'],
     "Isomers": "Generate Lowest Energy",
     "Optimisation_Choice": "True",
     "Total_Charge": 0,
     "Random_Seed": "23"},
]



if __name__ == "__main__":

    # Here the user specifies which input they want
    USER_INPUT = Integration_test
    RCA = PlacementRotation

    for i in range(len(USER_INPUT)):
        # We then take our input dictionary and create all the input variables from it
        Batch_name, Ligand_json, Assembled_Complex_json, Max_Num_Assembled_Complexes, Generate_Isomer_Instruction, Optimisation_Instruction, Random_Seed, Total_Charge, metal_list, topology_list = RCA.input_controller(
            USER_INPUT[i])
        print(f"\n\n\n\n**********Batch: {Batch_name}**********\n\n\n\n")
        # Here we load the ligand database
        F = LigandDB.from_json(json_=Ligand_json, type_="Ligand")  # We initiate the database in the RandomComplexAssembler
        RCA = PlacementRotation(database=F)
        ########################################################################
        # This section of code needs to be commented out if you do not want to
        # generate every combination of ligands but would rather build in a
        # random manner

        """IterativeLigands = ChooseIterativeLigands(database=RCA.ligand_dict,
                                                  top_list=topology_list,
                                                  charge=Total_Charge,
                                                  metal=metal_list,
                                                  random_seed=Random_Seed)
        IterativeLigands.obtain_all_combinations()"""
        ########################################################################

        j = 0  # Assembly iteration we are on
        Sum_Assembled_Complexes = 0  # Number of assembled complexes produced
        # Note: j is not always equal to Sum_Assembled_Complexes
        while Sum_Assembled_Complexes < Max_Num_Assembled_Complexes:
            print(f"###############################__Attempting_Assembly_of_Complex_#_{str(Random_Seed + j)}__###############################")
            #
            #
            # 1. Choose Random seed
            random.seed(int(Random_Seed) + j)

            #
            #
            # 3. Random Choice of Inputs
            assembly_instruction = random.choice(topology_list)
            Topology, Similarity = RCA.format_topologies(assembly_instruction)
            metal_centre = random.choice(metal_list)
            metal_type = metal_centre[0]
            metal_ox_state = str(metal_centre[1])
            metal_spin = int(metal_centre[2])

            #
            #
            # 3a. Choose Ligands Randomly
            # This section of code needs to be uncommented if you want complexes assembled randomly
            ligands = ChooseRandomLigands(database=RCA.ligand_dict,
                                          topology=Topology,
                                          instruction=Similarity,
                                          metal_oxidation_state=int(metal_ox_state),
                                          total_complex_charge=Total_Charge,
                                          max_attempts=1000000).choose_ligands()

            #
            #
            # 3b. Choose Ligands Iteratively
            # This section of code needs to be commented if you want complexes assembled randomly
            """ligands = IterativeLigands.choose_ligands(iteration=Random_Seed + j)"""

            #
            #
            # 4. Detect Tridentate Non-Planar
            if not RCA.planar_check_(ligands):
                warnings.warn("!!!Warning!!! -> The program has encountered a non_planar tridentate ligand, The program as of yet can not assemble with this ligand, skipping to next assembly")
                j += 1
                continue
            else:
                pass

            hydride_found = False
            for ligand in ligands.values():
                if ((ligand.atomic_props['atoms'][0] == "H") or (ligand.atomic_props['atoms'][0] == "Se")) and (len(ligand.atomic_props['atoms']) == 1):
                    hydride_found = True
                else:
                    pass

            if hydride_found:
                warnings.warn("!!!Warning!!! -> The program has encountered a hydride, The program as of yet can not assemble with this ligand, skipping to next assembly")
                j += 1
                continue

            #
            #
            # 5. Obtain rotated building blocks
            # Here we pass in our ligands and get out our stk building blocks
            stk_ligand_building_blocks_list, denticities = RCA.convert_ligand_to_building_block_for_complex(ligands=ligands,
                                                                                                            topology=Topology,
                                                                                                            metal=metal_type)

            #
            #
            # 6. Generate Isomers
            Isomers = BuildIsomers(topology=assembly_instruction,
                                   building_blocks_list=stk_ligand_building_blocks_list,
                                   metal_input=metal_type,
                                   charge_input=metal_ox_state,
                                   denticity_list=denticities,
                                   return_all_isomers=Generate_Isomer_Instruction,
                                   opt_choice=Optimisation_Instruction,
                                   ligand_list=ligands)

            Assembled_Complex_list, Building_Block_list = Isomers.Generate()

            #
            #
            # 7. Post-Process
            # Post process includes error detection and optimization
            Post_Process_Complex_List = []
            print("entering post process")
            for complex, building_blocks in zip(Assembled_Complex_list, Building_Block_list):
                complex_built = PostFilter(isomer=complex,
                                           metal_centre=metal_type,
                                           metal_offset=-0.9,  # The more negative the number the more crap gets through
                                           ligand_atom_offset=0.5,
                                           building_blocks=building_blocks,
                                           instruction=Optimisation_Instruction).closest_distance()

                output = OPTIMISE(isomer=complex_built,
                                  ligand_list=ligands,
                                  building_blocks=building_blocks,
                                  instruction=Optimisation_Instruction).Optimise_STK_Constructed_Molecule()

                complex_built = output[0]
                building_blocks_list_opt = output[1]

                complex_normal_built = PostFilter(isomer=complex_built,
                                                  metal_centre=metal_type,
                                                  building_blocks=building_blocks_list_opt,
                                                  metal_offset=-0.9,
                                                  ligand_atom_offset=0.5,
                                                  instruction=Optimisation_Instruction).post_optimisation_filter()
                Post_Process_Complex_List.append(complex_normal_built)

                if complex_normal_built is not None:
                    Sum_Assembled_Complexes += 1
            print("Leaving post process")

            #
            #
            # 8. Format Outputs
            RCA.output_controller_(list_of_complexes_wih_isomers=Post_Process_Complex_List,
                                   ligands=ligands,
                                   metal=metal_type,
                                   metal_ox_state=int(metal_ox_state),
                                   output_path=Assembled_Complex_json,
                                   metal_multiplicity=metal_spin,
                                   view_complex=False,
                                   concatonate_xyz=True,
                                   concatonate_xyz_name="INTEGRATION_TEST.xyz",
                                   write_gaussian_input_files=False,
                                   output_directory=USER_INPUT[i]['Output_Path'],
                                   frames=1)

            if j == 138746543956439563475683496736:
                exit()
            j += 1

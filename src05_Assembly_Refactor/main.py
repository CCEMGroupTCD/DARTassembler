from src01.DataBase import LigandDB
from src05_Assembly_Refactor.Assemble import PlacementRotation
from src05_Assembly_Refactor.ligands import ChooseRandomLigands, ChooseIterativeLigands
from src05_Assembly_Refactor.Isomer import BuildIsomers
from src05_Assembly_Refactor.Optimise import OPTIMISE
from src05_Assembly_Refactor.Post_Filter import PostFilter
import random
import warnings

warnings.simplefilter("always")

test_batch_list = [
    {'Name': 'first test batch',
     "Input_Path": "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/data/Filtered_Jsons/filteredLigDB_200423_FRANK_NN_donor_CHN_composition_Cu_precedence_with_pentafluorophenyl.json",
     'Output_Path': '/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/',
     'MAX_num_complexes': '100',
     'Topology_1': '[2, 0]--[1, 2]',
     'Metal_1': ['Cu', '+1', '0'],
     "Isomers": "Generate Lowest Energy",
     "Optimisation_Choice": "False",
     "Random_Seed": "0",
     "Total_Charge": 0}]

CASINI_Batch = [
    {'Name': 'CASINI',
     "Input_Path": "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/data/Filtered_Jsons/filteredLigDB_170423_CASINI_FINAL_ligand_db_V_1_6.json",
     'Output_Path': '/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/',
     'MAX_num_complexes': '100',
     'Topology_1': '[2, 1, 1]--[1, 2, 2]',
     'Metal_1': ['Au', '+3', '0'],
     "Isomers": "Generate Lowest Energy",
     "Optimisation_Choice": "False",
     "Total_Charge": 0,
     "Random_Seed": "0"}]
#todo: add a very small box to the inside of the horse shoe


if __name__ == "__main__":

    USER_INPUT = CASINI_Batch
    RCA = PlacementRotation

    for i in range(len(USER_INPUT)):
        Batch_name, Ligand_json, Assembled_Complex_json, Max_Num_Assembled_Complexes, Generate_Isomer_Instruction, Optimisation_Instruction, Random_Seed, Total_Charge, metal_list, topology_list = RCA.input_controller(USER_INPUT[i])

        F = LigandDB.from_json(json_=Ligand_json, type_="Ligand")  # We initiate the database in the RandomComplexAssembler
        RCA = PlacementRotation(database=F)

        IterativeLigands = ChooseIterativeLigands(database=RCA.ligand_dict,
                                                  top_list=topology_list,
                                                  charge=Total_Charge,
                                                  metal=metal_list,
                                                  random_seed=Random_Seed)
        IterativeLigands.obtain_all_combinations()

        j = 0
        Sum_Assembled_Complexes = 0
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
            """ligands = ChooseRandomLigands(database=RCA.ligand_dict,
                                          topology=Topology,
                                          instruction=Similarity,
                                          metal_oxidation_state=int(metal_ox_state),
                                          total_complex_charge=Total_Charge,
                                          max_attempts=1000000).choose_ligands()"""

            #
            #
            # 3b. Choose Ligands Iteratively
            ligands = IterativeLigands.choose_ligands(iteration=Random_Seed + j)


            #
            #
            # 4. Detect Tridentate Non-Planar
            if not RCA.planar_check_(ligands):
                warnings.warn("!!!Warning!!! -> The program has encountered a non_planar tridentate ligand, The program as of yet can not assemble with this ligand, skipping to next assembly")
                j += 1
                continue
            else:
                pass

            #
            #
            # 5. Obtain rotated building blocks
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
            Post_Process_Complex_List = []
            print("entering post process")
            for complex, building_blocks in zip(Assembled_Complex_list, Building_Block_list):

                complex_built = PostFilter(isomer=complex,
                                           metal_centre=metal_type,
                                           metal_offset=-0.2,
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
                                                  metal_offset=-0.2,
                                                  instruction=Optimisation_Instruction).post_optimisation_filter()
                
                Post_Process_Complex_List.append(complex_built)
                #todo delte line below
                #Post_Process_Complex_List.append(complex)
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
                                   concatonate_xyz_name="CASINI_COMPLEX_all.xyz",
                                   write_gaussian_input_files=True,
                                   output_directory="/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/CASINI_030423_calcs",
                                   frames=1)

            if j == 138746543956439563475683496736:
                exit()
            j += 1

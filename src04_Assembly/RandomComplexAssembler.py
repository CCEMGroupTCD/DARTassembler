import random
import mendeleev as atom
from src03_Assembly_Cian.building_block_utility import mercury_remover
import logging
import stk
import numpy as np
import os
import re
from pymatgen.core.periodic_table import Element
from src03_Assembly_Cian.building_block_utility import rotate_tridentate_bb, rotate_tetradentate_bb, penta_as_tetra, \
    get_optimal_rotation_angle_tridentate, Bidentate_Rotator, nonplanar_tetra_solver, get_energy_stk
from src03_Assembly_Cian.stk_utils import create_placeholder_Hg_bb
import src03_Assembly_Cian.stk_extension as stk_e
from copy import deepcopy
from src01.DataBase import LigandDB
from src01.Molecule import RCA_Ligand
from rdkit.Chem import rdmolfiles
from openbabel import openbabel as ob


class RandomComplexAssembler:
    """
    is kind of the hub to handle the random assembly.
    Stores the ligand database and stores the configuration for the random assembly
    """

    # def __init__(self, database: LigandDB, store_path: str = "../data/Assembled_Molecules"): If somethig breaks in the assembly it may be because I have replaced this line with the one below
    def __init__(self, database: LigandDB, store_path: str = "../data/Assembled_Molecules"):
        self.ligand_dict = database.get_lig_db_in_old_format()
        self.store_path = store_path

    @staticmethod
    def create_metal_building_block(metal, charge) -> stk.BuildingBlock:  # Creates the metal building block

        # build the metal block with the new metal atom
        # This may not work for neutral metals
        if str(charge) != "0":
            smiles_str = f"[{metal}{charge}]"
        elif str(charge) == "0":
            smiles_str = f"[{metal}]"
        stk_metal_func = getattr(stk, metal)
        functional_groups = (stk.SingleAtom(stk_metal_func(0, charge=charge)) for i in range(6))
        final_metal_bb = stk.BuildingBlock(smiles=smiles_str,
                                           functional_groups=functional_groups,
                                           position_matrix=np.ndarray([0, 0, 0])
                                           )

        return final_metal_bb

    @staticmethod
    def planar_check(ligands):  # Check if ligands are planar or not
        """
        If a tri or tetradentate ligand is contained in the topology, we need to evaluate whether this is
        a planar one or not.
        returns True if yes and false otherwise
        """
        for key, lig in ligands.items():
            if lig.denticity == 4 or lig.denticity == 3:
                return lig.planar_check()

        return True

    def convert_ligand_to_building_block_for_complex(self, ligands: dict[RCA_Ligand], topology) -> dict:
        # Here we pick and choose are ligands and rotate and place them based on our topology
        topology_determining_ligand_planar = self.planar_check(ligands)  # Check are either the tetra or tri ligands planar
        # print("topology_determining_ligand_planar " + str(topology_determining_ligand_planar))
        topology_list = topology  # I changed this because in the previous version it seemed that the topologies were being
        # print("topology_list " + str(topology_list))
        if type(topology_list[-1]) == list:
            topology_list.pop()
        else:
            pass

        # This ensures we don't enter the same if statement twice if we have to place a ligand of the same denticity twice
        first_lig0_placed = False
        first_lig1_placed = False
        first_lig2_placed_a = False
        first_lig2_placed_b = False

        ligand_buildingblocks = {}  # This will be updated with our ligand building blocks
        ligand_denticities = {}
        for i, ligand in enumerate(ligands.values()):

            # print("The number of atoms in this ligand are: " + str(len(ligand.atomic_props["atoms"])))
            if ligand.denticity == 0:
                # print("in denticity 0")
                # This is the intermediate of interest, a user predefined monodentate ligand that is a reaction intermediate
                # Below we dictate the rotations that must be carried out based on the selected topology
                monodentate_topology = stk_e.Monodentate(metals=create_placeholder_Hg_bb(),
                                                         ligands=ligand.to_stk_bb())
                if (topology_list == [4, 1, 0] and (topology_determining_ligand_planar is True)) or (topology_list == [3, 2, 0]):
                    monodentate_topology = stk_e.Monodentate_Top(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                    single_atom_position = [0, 0, 1.9]
                    # monodentate_topology.set_ligand_coordinates(coordinates=np.array([0, 0, 1.9]))

                elif (topology_list == [4, 1, 0]) and (topology_determining_ligand_planar is False):
                    monodentate_topology = stk_e.Monodentate_Back_Left(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                    single_atom_position = [-1.2, 1.2, 0]
                    # monodentate_topology.set_ligand_coordinates(coordinates=np.array([-1.2, 1.2, 0]))
                    first_lig0_placed = True

                elif ((topology_list == [4, 1, 0]) and (topology_determining_ligand_planar is False) and (first_lig0_placed == True)) or (topology_list == [2, 1, 0]):
                    monodentate_topology = stk_e.Monodentate_Front_Left(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                    single_atom_position = [-1.2, -1.2, 0]
                    # monodentate_topology.set_ligand_coordinates(coordinates=np.array([-1.2, -1.2, 0]))

                elif topology_list == [5, 0]:
                    monodentate_topology = stk_e.Monodentate_Bottom(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                    single_atom_position = [0, 0, 1.9]

                else:
                    single_atom_position = None
                    print("!!!Fatal Error!!! -> Your newly created topology {}, has not been accounted for in the assembly process (denticity = 0) -> Exiting Program ...".format(topology_list))
                    exit()

                # print("before")

                bb_for_complex = stk.BuildingBlock.init_from_molecule(stk.ConstructedMolecule(
                    topology_graph=monodentate_topology),
                    functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())]
                )
                # This function below is to handle monodentate ligands with only one atom

                # print("after")
                if len(ligand.atomic_props["atoms"]) == 1:
                    # print("we have decided we have 1 atom")

                    old_position_matrix = bb_for_complex.get_position_matrix()
                    old_position_matrix[1] = single_atom_position
                    new_position_matrix = old_position_matrix
                    bb_for_complex = bb_for_complex.with_position_matrix(new_position_matrix)
                else:

                    # print("we have decided we have multiple atoms")

                    pass



            elif ligand.denticity == 1:
                # If our ligand has denticity of 1 we enter this if statement

                # print("in denticity 1")
                # print(type(topology_list))

                if ((((topology_list == [4, 1, 1]) and (topology_determining_ligand_planar is False)) or (topology_list == [2, 1, 1])) and (first_lig1_placed == False)) or (
                        topology_list == [2, 1, 0]):
                    single_atom_position = [-1.2, 1.2, 0]
                    # The class after stk_e is called Monodentate_Back_Left which describes where the monodentate ligand will be placed
                    monodentate_topology = stk_e.Monodentate_Back_Left(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                    first_lig1_placed = True  # The first monodentate has been placed, so the next time we place a monodentate we don't want to place it here

                elif ((((topology_list == [4, 1, 1]) or (topology_list == [4, 1, 0])) and (topology_determining_ligand_planar is False)) or (topology_list == [2, 1, 1])) and (
                        first_lig1_placed == True):
                    single_atom_position = [-1.2, -1.2, 0]
                    monodentate_topology = stk_e.Monodentate_Front_Left(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                    # monodentate_topology.set_ligand_coordinates(coordinates=np.array([-1.2, 1.2, 0])
                elif ((topology_list == [4, 1, 1]) and (topology_determining_ligand_planar is True)) and (first_lig1_placed == False) or (topology_list == [3, 2, 1]):
                    single_atom_position = [0, 0, 1.9]
                    monodentate_topology = stk_e.Monodentate_Top(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                    # monodentate_topology.set_ligand_coordinates(coordinates=np.array([0, 0, 1.9]))
                    first_lig1_placed = True

                elif ((topology_list == [4, 1, 1] and (first_lig1_placed == True)) or (topology_list == [4, 1, 0])) and (topology_determining_ligand_planar is True) or (topology_list == [5, 1]):
                    monodentate_topology = stk_e.Monodentate_Bottom(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                    single_atom_position = [0, 0, -1.9]
                    # monodentate_topology.set_ligand_coordinates(coordinates=np.array([0, 0, -1.9]))

                else:
                    single_atom_position = None
                    monodentate_topology = None
                    print("!!!Fatal Error!!! -> Your newly created topology {}, has not been accounted for in the assembly process (denticity = 1) -> Exiting Program ...".format(topology_list))
                    exit()

                bb_for_complex = stk.BuildingBlock.init_from_molecule(stk.ConstructedMolecule(
                    topology_graph=monodentate_topology),
                    functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())]
                )
                if len(ligand.atomic_props["atoms"]) == 1:
                    # todo: get rid of this for single atoms

                    # Please Note that this is purely a temporary fix for ligands with only one atom, this is not intended to be in the final version
                    # print("we have decide we have 1 atom")

                    old_position_matrix = bb_for_complex.get_position_matrix()
                    old_position_matrix[1] = single_atom_position
                    new_position_matrix = old_position_matrix
                    bb_for_complex = bb_for_complex.with_position_matrix(new_position_matrix)
                else:

                    # print("we have decide we have multiple atoms")

                    pass




            elif ligand.denticity == 4:

                # print("in denticity 4")
                building_block = ligand.to_stk_bb()

                if topology_determining_ligand_planar is True:
                    # print("planar tetradentate")

                    # Then some rotation needs to be done
                    building_block = rotate_tetradentate_bb(building_block, ligand_=ligand)

                    tetra_topology_graph = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(),
                                                                       ligands=building_block
                                                                       )

                    bb_for_complex = stk.BuildingBlock.init_from_molecule(
                        stk.ConstructedMolecule(topology_graph=tetra_topology_graph),
                        functional_groups=[
                            stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), )]
                    )

                    # print(list(bb_for_complex.get_atoms()))
                    # print(bb_for_complex.get_position_matrix())
                    # print("finished qefiuhbvlar")


                elif topology_determining_ligand_planar is False:
                    # todo: This throws still errors
                    new_bb_path = nonplanar_tetra_solver(bb=building_block,
                                                         lig=ligand)

                    bb_for_complex = stk.BuildingBlock.init_from_file(new_bb_path, functional_groups=[
                        stk.SmartsFunctionalGroupFactory(
                            smarts='[Hg]', bonders=(0,),
                            deleters=(), ), ], )
                else:
                    bb_for_complex = None
                    print("!!!Fatal Error!!! -> the denticity of the tetradentate ligand is unresolved")
                    exit()




            elif ligand.denticity == 3:

                # print("in denticity 3")

                building_block = ligand.to_stk_bb()

                if topology_determining_ligand_planar is False:
                    print("Not implemented yet, we should not end up here")
                    raise NotImplementedError

                elif topology_determining_ligand_planar is True:
                    building_block = rotate_tridentate_bb(tridentate_bb_=building_block, ligand_=ligand)

                    tridentate_toplogy = stk_e.Tridentate(metals=create_placeholder_Hg_bb(),
                                                          ligands=building_block
                                                          )

                    compl_constructed_mol = stk.ConstructedMolecule(topology_graph=tridentate_toplogy)

                    compl_constructed_mol = compl_constructed_mol.with_rotation_about_axis(
                        axis=np.array((0, 0, 1)),
                        angle=float(np.radians(
                            get_optimal_rotation_angle_tridentate(compl_constructed_mol, 10.0, 0.0, 0.0, ligand))),
                        origin=np.array((0, 0, 0))
                    )

                    position_matrix = compl_constructed_mol.get_position_matrix()

                    position_matrix[0] = [-0.8, 0, 0]  # was previously -0.6

                    compl_constructed_mol = compl_constructed_mol.with_position_matrix(position_matrix=position_matrix)

                    if topology_list == [3, 2, 0] or [3, 2, 1]:
                        bb_for_complex = stk.BuildingBlock.init_from_molecule(compl_constructed_mol, functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])
                    elif topology_list == [3, 3]:  # This arrangement yet is not implemented so this elif statement should never be entrered
                        # todo: Might implement this in the future but also might not
                        bb_for_complex_pre = stk.BuildingBlock.init_from_molecule(compl_constructed_mol,
                                                                                  functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])
                        bb_for_complex = bb_for_complex_pre.with_rotation_about_axis(angle=90.0 * (np.pi / 180.0), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)), )
                    else:
                        bb_for_complex = None
                        print("!!!Fatal Error!!! -> the geometry of the tridentate ligand is unresolved")
                        exit()
                else:
                    bb_for_complex = None
                    print("!!!Fatal Error!!! -> the geometry of the tridentate ligand is unresolved")
                    exit()



            elif ligand.denticity == 2:

                # print("in denticity 2")

                if topology_list == [3, 2, 0]:
                    bidentate_topology = stk_e.Bidentate(metals=create_placeholder_Hg_bb(),
                                                         ligands=ligand.to_stk_bb()
                                                         )
                elif (topology_list == [2, 2] or [2, 1, 1] or [2, 1, 0]) and (first_lig2_placed_a == False):
                    bidentate_topology = stk_e.Bidentate_Planar_Right(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                    first_lig2_placed_a = True

                elif (topology_list == [2, 2] or [2, 1, 1] or [2, 1, 0]) and (first_lig2_placed_a == True):
                    bidentate_topology = stk_e.Bidentate_Planar_Left(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
                else:
                    bidentate_topology = None
                    print("!!!Fatal_Error!!! something wrong with the bidentate stuff (a)")
                    exit()

                complex_bidentate = stk.ConstructedMolecule(topology_graph=bidentate_topology)

                if not first_lig2_placed_b:
                    new_mol_path = Bidentate_Rotator(ligand_bb=complex_bidentate, ligand=ligand, top_list=topology_list, bool_placed=first_lig2_placed_b)
                    first_lig2_placed_b = True
                elif first_lig2_placed_b:
                    new_mol_path = Bidentate_Rotator(ligand_bb=complex_bidentate, ligand=ligand, top_list=topology_list, bool_placed=first_lig2_placed_b)

                else:
                    print("!!!Fatal_Error!!! -> Something wrong with the bidentate stuff (b)")
                    exit()
                bb_for_complex = stk.BuildingBlock.init_from_file(new_mol_path, functional_groups=[
                    stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ], )




            elif ligand.denticity == 5:
                # This elif statement shouldn't need much modifying.

                # print("in denticity 5")

                tetra_bb_for_penta, position_index = penta_as_tetra(ligand=ligand)

                tetra_bb_for_penta = rotate_tetradentate_bb(tetra_bb_for_penta, ligand)

                tip_position = list(tetra_bb_for_penta.get_atomic_positions(atom_ids=[int(position_index), ]))

                if float(tip_position[0][2]) > 0:
                    # Additional rotation is required
                    tetra_bb_for_penta = tetra_bb_for_penta.with_rotation_about_axis(angle=np.radians(180), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)))
                elif float(tip_position[0][2]) < 0:
                    # No rotation is required
                    pass
                else:
                    print("Tip position for penta = 0, unexpected")
                    raise ValueError

                penta_topology = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(),
                                                             ligands=tetra_bb_for_penta
                                                             )

                bb_for_complex = stk.BuildingBlock.init_from_molecule(
                    stk.ConstructedMolecule(topology_graph=penta_topology),
                    functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ]
                )

            else:

                print("!!!Fatal Error!!! -> Unknown Ligand Denticity -> Aborting Program")
                raise ValueError

            ligand_buildingblocks[i] = bb_for_complex
            ligand_denticities[i] = ligand.denticity

        return ligand_buildingblocks, ligand_denticities  # We essentially build a list of building blocks that we need to assemble our complex

    def create_random_TMC(self, Metal_input_list, topology_input_list, isomer_build_status, opt_decision):
        # From The Assembly Class this is now the entry point in this file
        chosen_metal_list = random.choice(Metal_input_list)  # We pick a metal
        metal = chosen_metal_list[0]  # We pick a metal
        metal_ox = chosen_metal_list[1]  # The metals OX
        # print("The Metal OX state in create random tmc is :" + str(metal_ox))
        spin_state = str(chosen_metal_list[2])  # extract spin state
        multiplicity = int(self.calculate_spin(metal, metal_ox, spin_state))
        # print("multiplicity: " + str(multiplicity))
        comp = deepcopy((random.choice(topology_input_list)))  # comp is short for composition and is just another word for topology
        complete_topology_nomenclature = deepcopy(str(comp))  # This may be redundant now

        # print("comp: " + str(comp))
        # print("complete_topology_nomenclature: " + str(complete_topology_nomenclature))

        # The following code  is used to process the topology string, its main role is to ensure that the same random
        # ligand is used, for example [2,2,[1,1]]
        if type(comp[-1]) == list:
            similarity_list = comp[-1]
            comp.pop()
            ligands = {i: random.choice(self.ligand_dict[index]) for i, index in enumerate(comp)}  # ligands chosen based on the toplogy
            passed_similarity_codes = []
            i = 0
            # This for loop looks at the final list to determine if any ligands should be identical
            for similarity_code_1 in similarity_list:
                j = 0
                for similarity_code_2 in similarity_list:
                    if (similarity_code_1 == similarity_code_2) and (i != j) and not (
                            bool(any(item in passed_similarity_codes for item in passed_similarity_codes))):
                        ligands[similarity_list.index(similarity_code_2)] = ligands[int(similarity_code_1)]
                        passed_similarity_codes.append(similarity_code_1)
                    else:
                        pass
                    j = j + 1
                i = i + 1
        elif type(int(comp[-1])) == int:
            if (0 not in self.ligand_dict.keys()) and (0 in comp):
                print("!!!Warning!!! -> Reaction intermediate not specified but is necessitated by the chosen topology -> replacing with monodentate ligand")

                i = 0
                for item in comp:
                    if item == 0:
                        comp[i] = 1
                    else:
                        pass

                    i = i + 1
            ligands = {i: random.choice(self.ligand_dict[index]) for i, index in enumerate(comp)}  # ligands chosen based on the toplogy

        else:
            ligands = None
            print("!!!Fatal_Error!!! -> I would advise you check the .yml file and make sure everything is okay there")
            exit()

        planar = self.planar_check(ligands)  # is True for 5-0!

        if (comp == [3, 2, 0] and planar is False) or (comp == [4, 1, 0] and planar is False) or (comp == [4, 1, 1] and planar is False) or (comp == [3, 2, 1] and planar is False):
            # account for tetrahedral planar strucure ligands here
            if ((comp == [3, 2, 0]) or (comp == [3, 2, 1])) and planar is False:
                print("!!!Warning!!! -> The program has encountered a non_planar tridentate ligand, The program as of yet can not assemble with this ligand, skipping to next assembly")
            elif (comp == [4, 1, 0] and planar is False) or (comp == [4, 1, 1] and planar is False):
                # Currently the issue with tetradentate non-planar ligands is with the assembly.
                # The molecule class prevents the self.coordinates attribute being updated with the relevant
                # information from the add_atom method. The molecule class will need to be updated in the future to
                # allow for the number of atoms to be updated
                print("!!!Warning!!! -> The program has encountered a non_planar tetradentate ligand, The program as of yet can not assemble with this ligand, skipping to next assembly")
            complex_ = None
            logging.info("Not implemented yet")

            return complex_, ligands, metal, int(metal_ox), multiplicity, complete_topology_nomenclature

        else:
            # and thus we got everything to assemble the complex
            # Here we extract information regarding the functional groups (currently used exclusively for the optimization process)
            functional_groups_str = {key_: lig.get_assembly_dict()["str"] for key_, lig in ligands.items()}
            functional_groups_index = {key_: lig.get_assembly_dict()["index"] for key_, lig in ligands.items()}
            functional_groups_type = {key_: lig.get_assembly_dict()["type"] for key_, lig in ligands.items()}
            # Next we start picking and choosing our ligands
            building_blocks, denticities = self.convert_ligand_to_building_block_for_complex(ligands, comp, )

            # print("comp: " + str(comp))
            # print("complete_topology_nomenclature: " + str(complete_topology_nomenclature))
            complex_ = (
                self.isomer_handler(topology=str(complete_topology_nomenclature), building_blocks_list=building_blocks, metal_input=metal, charge_input=metal_ox, denticity_list=denticities,
                                    return_all_isomers=isomer_build_status, opt_choice=opt_decision, ligand_list=ligands))
            return complex_, ligands, metal, int(metal_ox), multiplicity, complete_topology_nomenclature

    def isomer_handler(self, topology, building_blocks_list, metal_input, charge_input, denticity_list, return_all_isomers, opt_choice, ligand_list):
        # ("This is the topology " + str(topology))
        # The goal of this function is to take in a complex and return all possible isomers of that complex

        #
        #
        # 1
        if (str(topology) == '[2, 1, 0]') or (str(topology) == "[2, 1, 1, ['1', '2', '2']]") or (str(topology) == "[2, 1, 1, ['1', '2', '3']]") or (str(topology) == "[2, 2, ['1', '1']]") or (
                str(topology) == "[2, 2, ['1', '2']]"):
            # print("we are going to rotate one of the bidentates")
            # I expect only one isomer from this if statement (so a list of two constructed molecules). Only a single bidentate needs to be flipped 180
            bidentate_already_rotated = False  # This is to ensure we don't enter the same if statement twice
            building_blocks_rotated = {}  # This will contain all the ligands for the isomer of our complex
            for j in range(len(building_blocks_list)):
                # We iterate through all the complexes in our input complex
                if (denticity_list[j] == 2) and (bidentate_already_rotated == False):

                    # If we come across a bidentate ligand and it's the first time we have come across one then ...

                    bidentate_already_rotated = True
                    # Then we flip 180 degrees
                    building_blocks_rotated[j] = building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)), )
                else:
                    # Otherwise we just copy over the ligand from the previous complex
                    building_blocks_rotated[j] = building_blocks_list[j]

            if len(building_blocks_list) == 2:
                # If there are only two ligands that comprise our complex then ...
                # This is our input complex
                complex_normal = stk_e.complex_topology_two(metals=self.create_metal_building_block(metal_input, charge_input),
                                                            ligands={building_block: (i,) for i, building_block in
                                                                     building_blocks_list.items()}
                                                            )

                # This is an isomer
                complex_flipped = stk_e.complex_topology_two(metals=self.create_metal_building_block(metal_input, charge_input),
                                                             ligands={building_block: (i,) for i, building_block in
                                                                      building_blocks_rotated.items()}
                                                             )

            elif len(building_blocks_list) == 3:
                # If there are 3 ligands in our complex then ...
                complex_normal = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal_input, charge_input),
                                                              ligands={building_block: (i,) for i, building_block in
                                                                       building_blocks_list.items()}
                                                              )

                complex_flipped = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal_input, charge_input),
                                                               ligands={building_block: (i,) for i, building_block in
                                                                        building_blocks_rotated.items()}
                                                               )

            else:
                complex_normal = None
                complex_flipped = None
                print("!!!Fatal Error!!! -> you need to account for this topology {} in the isomer handler".format(str(topology)))
                exit()

            # Now we build our isomers, remove any dummy atoms and optimise if we so desire
            complex_normal_built = stk.ConstructedMolecule(topology_graph=complex_normal)
            complex_normal_built = mercury_remover(complex_normal_built)
            complex_normal_built, building_blocks_list = post_filter(isomer=complex_normal_built, building_blocks=building_blocks_list, metal=metal_input).closest_distance()
            complex_normal_built, building_blocks_list_opt = OPTIMISE(isomer=complex_normal_built, ligand_list=ligand_list, building_blocks=building_blocks_list).Optimise_STK_Constructed_Molecule() if (opt_choice == "True") else complex_normal_built, building_blocks_list
            complex_normal_built = complex_normal_built[0]
            complex_normal_built = post_filter(isomer=complex_normal_built, building_blocks=building_blocks_list_opt, metal=metal_input).post_optimisation_filter()

            complex_flipped_built = stk.ConstructedMolecule(topology_graph=complex_flipped)
            complex_flipped_built = mercury_remover(complex_flipped_built)
            complex_flipped_built, building_blocks_rotated = post_filter(isomer=complex_flipped_built, building_blocks=building_blocks_rotated, metal=metal_input).closest_distance()
            complex_flipped_built, building_blocks_rotated_opt = OPTIMISE(isomer=complex_flipped_built, ligand_list=ligand_list, building_blocks=building_blocks_rotated).Optimise_STK_Constructed_Molecule() if (opt_choice == "True") else complex_flipped_built, building_blocks_rotated
            complex_flipped_built = complex_flipped_built[0]
            complex_flipped_built = post_filter(isomer=complex_flipped_built, building_blocks=building_blocks_rotated_opt, metal=metal_input).post_optimisation_filter()

            if return_all_isomers == "Generate Lowest Energy":

                if (complex_flipped_built is not None) and (complex_normal_built is not None):
                    if get_energy_stk(complex_flipped_built) <= get_energy_stk(complex_normal_built):
                        return [complex_flipped_built]
                    elif get_energy_stk(complex_flipped_built) > get_energy_stk(complex_normal_built):
                        return [complex_normal_built]

                elif (complex_flipped_built is not None) and (complex_normal_built is None):
                    return [complex_flipped_built]

                elif (complex_flipped_built is None) and (complex_normal_built is not None):
                    return [complex_normal_built]

                elif (complex_flipped_built is None) and (complex_normal_built is None):
                    return None, None


            elif return_all_isomers == "Generate All":
                if (complex_flipped_built is not None) and (complex_normal_built is not None):
                    #print("Returning All Isomers")
                    #print("Cis Energy: " + str(get_energy_stk(complex_normal_built)))
                    #print("Trans Energy: " + str(get_energy_stk(complex_flipped_built)))
                    isomer_list = [complex_normal_built, complex_flipped_built]
                    return isomer_list
                elif (complex_flipped_built is None) and (complex_normal_built is not None):
                    return [complex_normal_built]
                elif (complex_flipped_built is not None) and (complex_normal_built is None):
                    return [complex_flipped_built]


            else:
                print("!!!Fatal error!!! -> Isomer Build status not as expected -> Exiting program ...")
                exit()


        #
        #
        # 2

        elif str(topology) == '[3, 2, 0]':
            # building_blocks_list -> good
            # These are the empty dictionaries that will contain all the building blocks for each isomer
            # This is by far the most 'involved' isomer generation processes
            building_blocks_rotated_bi = {}
            building_blocks_rotated_tri = {}
            building_blocks_rotated_bi_and_tri = {}
            for j in range(len(building_blocks_list)):
                if denticity_list[j] == 2:
                    building_blocks_rotated_bi[j] = building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((-1, 0, -1)), origin=np.array((0, 0, 0)), )
                    building_blocks_rotated_bi_and_tri[j] = building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((-1, 0, -1)), origin=np.array((0, 0, 0)), )
                    building_blocks_rotated_tri[j] = building_blocks_list[j]
                    pass
                elif denticity_list[j] == 3:
                    building_blocks_rotated_tri[j] = building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)), )
                    building_blocks_rotated_bi_and_tri[j] = building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)), )
                    building_blocks_rotated_bi[j] = building_blocks_list[j]
                else:
                    building_blocks_rotated_bi[j] = building_blocks_list[j]
                    building_blocks_rotated_tri[j] = building_blocks_list[j]
                    building_blocks_rotated_bi_and_tri[j] = building_blocks_list[j]
                    pass

            # Here are the four potential isomers
            complex_normal = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal_input, charge_input),
                                                          ligands={building_block: (i,) for i, building_block in
                                                                   building_blocks_list.items()}
                                                          )

            complex_rotated_bi = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal_input, charge_input),
                                                              ligands={building_block: (i,) for i, building_block in
                                                                       building_blocks_rotated_bi.items()}
                                                              )

            complex_rotated_tri = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal_input, charge_input),
                                                               ligands={building_block: (i,) for i, building_block in
                                                                        building_blocks_rotated_tri.items()}
                                                               )

            complex_rotated_bi_and_tri = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal_input, charge_input),
                                                                      ligands={building_block: (i,) for i, building_block in
                                                                               building_blocks_rotated_bi_and_tri.items()}
                                                                      )

            complex_normal_built = stk.ConstructedMolecule(topology_graph=complex_normal)
            complex_normal_built = mercury_remover(complex_normal_built)
            complex_normal_built, building_blocks_list = post_filter(isomer=complex_normal_built, building_blocks=building_blocks_list, metal=metal_input).closest_distance()
            complex_normal_built, building_blocks_list_opt = OPTIMISE(isomer=complex_normal_built, ligand_list=ligand_list, building_blocks=building_blocks_list).Optimise_STK_Constructed_Molecule() if (opt_choice == "True") else complex_normal_built, building_blocks_list
            complex_normal_built = complex_normal_built[0]
            complex_normal_built = post_filter(isomer=complex_normal_built, building_blocks=building_blocks_list_opt, metal=metal_input).post_optimisation_filter()

            complex_rotated_bi_built = stk.ConstructedMolecule(topology_graph=complex_rotated_bi)
            complex_rotated_bi_built = mercury_remover(complex_rotated_bi_built)
            complex_rotated_bi_built, building_blocks_rotated_bi = post_filter(isomer=complex_rotated_bi_built, building_blocks=building_blocks_rotated_bi, metal=metal_input).closest_distance()
            complex_rotated_bi_built, building_blocks_rotated_bi_opt = OPTIMISE(isomer=complex_rotated_bi_built, ligand_list=ligand_list, building_blocks=building_blocks_rotated_bi).Optimise_STK_Constructed_Molecule() if (opt_choice == "True") else complex_rotated_bi_built, building_blocks_rotated_bi
            complex_rotated_bi_built = complex_rotated_bi_built[0]
            complex_rotated_bi_built = post_filter(isomer=complex_rotated_bi_built, building_blocks=building_blocks_rotated_bi_opt, metal=metal_input).post_optimisation_filter()

            complex_rotated_tri_built = stk.ConstructedMolecule(topology_graph=complex_rotated_tri)
            complex_rotated_tri_built = mercury_remover(complex_rotated_tri_built)
            complex_rotated_tri_built, building_blocks_rotated_tri = post_filter(isomer=complex_rotated_tri_built, building_blocks=building_blocks_rotated_tri, metal=metal_input).closest_distance()
            complex_rotated_tri_built, building_blocks_rotated_tri_opt = OPTIMISE(isomer=complex_rotated_tri_built, ligand_list=ligand_list, building_blocks=building_blocks_rotated_tri).Optimise_STK_Constructed_Molecule() if (opt_choice == "True") else complex_rotated_tri_built, building_blocks_rotated_tri
            complex_rotated_tri_built = complex_rotated_tri_built[0]
            complex_rotated_tri_built = post_filter(isomer=complex_rotated_tri_built, building_blocks=building_blocks_rotated_tri_opt, metal=metal_input).post_optimisation_filter()

            complex_rotated_bi_and_tri_built = stk.ConstructedMolecule(topology_graph=complex_rotated_bi_and_tri)
            complex_rotated_bi_and_tri_built = mercury_remover(complex_rotated_bi_and_tri_built)
            complex_rotated_bi_and_tri_built, building_blocks_rotated_bi_and_tri = post_filter(isomer=complex_rotated_bi_and_tri_built, building_blocks=building_blocks_rotated_bi_and_tri, metal=metal_input).closest_distance()
            complex_rotated_bi_and_tri_built, building_blocks_rotated_bi_and_tri_opt = OPTIMISE(isomer=complex_rotated_bi_and_tri_built, ligand_list=ligand_list, building_blocks=building_blocks_rotated_bi_and_tri).Optimise_STK_Constructed_Molecule() if (opt_choice == "True") else complex_rotated_bi_and_tri_built, building_blocks_rotated_bi_and_tri
            complex_rotated_bi_and_tri_built = complex_rotated_bi_and_tri_built[0]
            complex_rotated_bi_and_tri_built = post_filter(isomer=complex_rotated_bi_and_tri_built, building_blocks=building_blocks_rotated_bi_and_tri_opt, metal=metal_input).post_optimisation_filter()

            isomer_list = [complex_normal_built, complex_rotated_bi_built, complex_rotated_tri_built, complex_rotated_bi_and_tri_built]
            building_blocks = [building_blocks_list, building_blocks_rotated_bi, building_blocks_rotated_tri, building_blocks_rotated_bi_and_tri]

            isomer_energy_list = [get_energy_stk(complex_normal_built), get_energy_stk(complex_rotated_bi_built), get_energy_stk(complex_rotated_tri_built), get_energy_stk(complex_rotated_bi_and_tri_built)]

            #for isomer, building_block, isomer_energy in zip(isomer_list, building_blocks, isomer_energy_list):
            del_list = []
            for j in range(len(isomer_list)):
                # We delete any Nones in our list
                if (isomer_list[j] is None) or (building_blocks[j] is None) or (isomer_energy_list[j] is None):
                    del_list.append(j)

            for del_ in reversed(del_list):
                del isomer_energy_list[del_]
                del isomer_list[del_]
                del building_blocks[del_]
            if not isomer_list:
                return [None]

            if return_all_isomers == "Generate Lowest Energy":
                smallest_energy = min(isomer_energy_list)

                for complex_, building_block in zip(isomer_list, building_blocks):
                    if abs(get_energy_stk(complex_) - smallest_energy) < 0.0001:  # are they approximately = (just in case there are floating point errors)
                        return [complex_]
                    else:
                        pass

            if return_all_isomers == "Generate All":

                return isomer_list

            else:
                print("!!!Fatal error!!! -> Isomer Build status not as expected")
                exit()


        #
        #
        # 3

        elif (str(topology) == "[4, 1, 1, ['1', '2', '3']]") or (str(topology) == "[4, 1, 0]"):
            # Here we are just flipping the tetradentate ligand
            building_blocks_rotated_tetra = {}
            for j in range(len(building_blocks_list)):
                if denticity_list[j] == 4:

                    # print(building_blocks_list[j].get_position_matrix())
                    building_blocks_rotated_tetra[j] = building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)), )
                else:
                    # print(building_blocks_list[j].get_position_matrix())
                    # print(list(building_blocks_list[j].get_atoms()))

                    building_blocks_rotated_tetra[j] = building_blocks_list[j]

            complex_normal = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal_input, charge_input),
                                                          ligands={building_block: (i,) for i, building_block in
                                                                   building_blocks_list.items()}
                                                          )

            complex_rotated_tetra = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal_input, charge_input),
                                                                 ligands={building_block: (i,) for i, building_block in
                                                                          building_blocks_rotated_tetra.items()}
                                                                 )

            complex_normal_built = stk.ConstructedMolecule(topology_graph=complex_normal)
            complex_normal_built = mercury_remover(complex_normal_built)
            complex_normal_built, building_blocks_list = post_filter(isomer=complex_normal_built, building_blocks=building_blocks_list, metal=metal_input).closest_distance()
            complex_normal_built, building_blocks_list_opt = OPTIMISE(isomer=complex_normal_built, ligand_list=ligand_list, building_blocks=building_blocks_list).Optimise_STK_Constructed_Molecule() if (opt_choice == "True")else complex_normal_built, building_blocks_list
            complex_normal_built = complex_normal_built[0]
            complex_normal_built = post_filter(isomer=complex_normal_built, building_blocks=building_blocks_list_opt, metal=metal_input).post_optimisation_filter()

            complex_rotated_tetra_built = stk.ConstructedMolecule(topology_graph=complex_rotated_tetra)
            complex_rotated_tetra_built = mercury_remover(complex_rotated_tetra_built)
            complex_rotated_tetra_built, building_blocks_rotated_tetra = post_filter(isomer=complex_rotated_tetra_built, building_blocks=building_blocks_rotated_tetra, metal=metal_input).closest_distance()
            complex_rotated_tetra_built, building_blocks_rotated_tetra_built_opt = OPTIMISE(isomer=complex_rotated_tetra_built, ligand_list=ligand_list, building_blocks=building_blocks_rotated_tetra).Optimise_STK_Constructed_Molecule() if (opt_choice == "True") else complex_rotated_tetra_built, building_blocks_rotated_tetra
            complex_rotated_tetra_built = complex_rotated_tetra_built[0]
            complex_rotated_tetra_built = post_filter(isomer=complex_rotated_tetra_built, building_blocks=building_blocks_rotated_tetra, metal=metal_input).post_optimisation_filter()

            if return_all_isomers == "Generate Lowest Energy":

                if (complex_rotated_tetra_built is not None) and (complex_normal_built is not None):
                    if get_energy_stk(complex_rotated_tetra_built) <= get_energy_stk(complex_normal_built):
                        return [complex_rotated_tetra_built]
                    elif get_energy_stk(complex_rotated_tetra_built) > get_energy_stk(complex_normal_built):
                        return [complex_normal_built]

                elif (complex_rotated_tetra_built is not None) and (complex_normal_built is None):
                    return [complex_rotated_tetra_built]

                elif (complex_rotated_tetra_built is None) and (complex_normal_built is not None):
                    return [complex_normal_built]

                elif (complex_rotated_tetra_built is None) and (complex_normal_built is None):
                    return None, None


            elif return_all_isomers == "Generate All":
                if (complex_rotated_tetra_built is not None) and (complex_normal_built is not None):
                    isomer_list = [complex_normal_built, complex_rotated_tetra_built]
                    return isomer_list
                elif (complex_rotated_tetra_built is None) and (complex_normal_built is not None):
                    return [complex_normal_built]
                elif (complex_rotated_tetra_built is not None) and (complex_normal_built is None):
                    return [complex_rotated_tetra_built]

        #
        #
        # 4

        else:  # This is for complexes that have no isomers
            # print("This particular complex has no isomers")

            if len(building_blocks_list) == 3:
                complex_top = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal_input, charge_input),
                                                           ligands={building_block: (i,) for i, building_block in
                                                                    building_blocks_list.items()}
                                                           )

            elif len(building_blocks_list) == 2:
                complex_top = stk_e.complex_topology_two(metals=self.create_metal_building_block(metal_input, charge_input),
                                                         ligands={building_block: (i,) for i, building_block in
                                                                  building_blocks_list.items()}
                                                         )

            else:
                complex_top = None

                # print("!!!Fatal Error!!! -> You need to update your Octahedral Constructor")

                exit()

            complex_built = stk.ConstructedMolecule(topology_graph=complex_top)
            complex_built = mercury_remover(complex_built)
            complex_built, building_blocks_list = post_filter(isomer=complex_built, building_blocks=building_blocks_list, metal=metal_input).closest_distance()
            complex_built, building_blocks_list = OPTIMISE(isomer=complex_built, ligand_list=ligand_list, building_blocks=building_blocks_list).Optimise_STK_Constructed_Molecule() if (opt_choice == "True") else complex_built, building_blocks_list
            complex_built = complex_built[0]

            if complex_built is not None:
                return [complex_built]
            else:
                return [None]

    @staticmethod
    def calculate_spin(metal_input, oxidation_state_input, spin_input):
        # This function is purely for calculating the multiplicity of a complex
        # print("Calculating Spin ...")
        # The key in the spins dictionary corresponds to d-electron count, the value corresponds to the spin
        spins = {0: [1],
                 1: [2],
                 2: [3],
                 3: [4],
                 4: [3, 5],
                 5: [2, 6],
                 6: [1, 5],
                 7: [2, 4],
                 8: [3],
                 9: [2],
                 10: [1]}
        # print("The Metal is: " + str(metal_input))
        # print("The oxidation state is: " + str(oxidation_state_input))
        # print("The group is: " + str(int(atom.element(str(metal_input)).group_id)))
        d_electron = int(atom.element(str(metal_input)).group_id) - int(oxidation_state_input)
        # print("The d-electron count is: " + str(d_electron))
        potential_multiplicities = spins[d_electron]
        # print("The potential_Multiplicities are:" + str(potential_multiplicities))
        if len(potential_multiplicities) == 1:
            # print("only one option for choosing multiplicity")
            spin_output = int(potential_multiplicities[0])
            pass
        else:
            if spin_input == "High":
                # print("high spin was chosen")
                spin_output = int(potential_multiplicities[1])
            elif spin_input == "Low":
                # print("Low spin was chosen")
                spin_output = int(potential_multiplicities[0])
            else:
                spin_output = None
                print("!!!Fatal Error!!! -> Spin State unaccounted for --> Exiting program")
                exit()
        # print("spin_output :" + str(spin_output))
        return spin_output


class post_filter:
    def __init__(self, isomer, building_blocks, metal):
        self.building_blocks = building_blocks
        self.isomer = isomer
        self.failed_isomers = []
        self.threshold = 0.3    #Angstrom todo: needs to be tuned

    def closest_distance(self):
        # This function will detect and filter out complexes that have clashing between atoms in different ligands but it WILL NOT detect clashing between atoms of the same
        # ligand or the metal centre
        # todo: detect clashing between atoms and the metal centre without taking into account the functional groups maybe an idea for the future
        for keys_1, values_1 in self.building_blocks.items():  # loop through the building blocks within each isomer
            values_1 = mercury_remover(values_1)  # Eliminate the temporary Mercury
            for keys_2, values_2 in self.building_blocks.items():  # loop through the building blocks for within each isomer
                values_2 = mercury_remover(values_2)  # Eliminate the temporary Mercury
                if keys_1 == keys_2:  # Don't compare anything if they are the same ligand
                    pass
                elif keys_1 != keys_2:  # Compare distance if the ligands are not the same ligand
                    for point_1, atom_1 in zip(values_1.get_position_matrix(), values_1.get_atoms()):  # we loop through all the positions of the atoms
                        atom_1_type = [str(atom_1).split("(")][0][0]
                        cov_1 = Element(atom_1_type).atomic_radius
                        for point_2, atom_2 in zip(values_2.get_position_matrix(), values_2.get_atoms()):  # we loop through all the positions of the atoms
                            atom_2_type = [str(atom_2).split("(")][0][0]
                            cov_2 = Element(atom_2_type).atomic_radius
                            distance = np.linalg.norm(point_1 - point_2)  # Calculate distance
                            if distance < (cov_1 + cov_2):  # This function shouldn't take into account ligand metal distances
                                print("PRE OPTIMISATION FILTER FAILED, RETURNING NONE")  # if clashing is present
                                return None, None
                            else:
                                pass
        return self.isomer, self.building_blocks

    def post_optimisation_filter(self):
        if self.building_blocks is None and self.isomer is None:
            print("NONE DETECTED AS INPUT IN POST OPTIMISATION FILTER")
            return None
        for keys_1, values_1 in self.building_blocks.items():
            #print(keys_1)
            for bond in list(values_1.get_bonds()):
                atom_1_id = bond.get_atom1().get_id()
                atom_2_id = bond.get_atom2().get_id()
                atom_1_pos = list(values_1.get_atomic_positions(atom_1_id))
                atom_2_pos = list(values_1.get_atomic_positions(atom_2_id))
                atom_1_AN = bond.get_atom1().get_atomic_number()
                atom_2_AN = bond.get_atom2().get_atomic_number()
                distance = np.linalg.norm(atom_1_pos[0] - atom_2_pos[0])
                #print(distance)
                if distance > ((Element.from_Z(atom_1_AN).atomic_radius + Element.from_Z(atom_2_AN).atomic_radius) + self.threshold):
                    print("POST FILTER FAILED, RETURNING NONE")
                    return None
                else:
                    pass
        return self.isomer



class OPTIMISE:
    def __init__(self, isomer, ligand_list, building_blocks):
        self.isomer = isomer
        self.ligands = ligand_list
        self.building_blocks = building_blocks

    @staticmethod
    def movie(input_file, working_directory):
        os.system('touch {}{}'.format(working_directory, "movie_loop.xyz"))
        os.system('cat {}{} >> {}{}'.format(working_directory, input_file, working_directory, "movie_loop.xyz"))

    def Optimise_STK_Constructed_Molecule(self):
        #print("1 " + str(type(self.isomer)))
        if (self.isomer is None) and (self.building_blocks is None):
            print("NONE DETECTED IN OPTIMISER")
            return None, None
        else:
            print("Beginning Optimisation")
           # print("2 " + str(type(self.isomer)))
            debug = False
            # todo: Return the rotated stk building blocks as well, they may still be of use to someone
            # REFERENCE: https://github.com/hjkgrp/molSimplify/blob/c3776d0309b5757d5d593e42db411a251b558e59/molSimplify/Scripts/structgen.py#L658
            # REFERENCE: https://gist.github.com/andersx/7784817

            # stk to xyz string
            complex_mol = self.isomer.to_rdkit_mol()
            xyz_string = rdmolfiles.MolToXYZBlock(complex_mol)

            # setup conversion
            conv = ob.OBConversion()
            conv.SetInAndOutFormats('xyz', 'xyz')
            mol = ob.OBMol()
            conv.ReadString(mol, xyz_string)

            # Define constraints
            # OPEN BABEL INDEXING SARTS AT ONE !!!!!!!!!!!!!!!!!!!!!!!!!!!
            constraints = ob.OBFFConstraints()
            constraints.AddAtomConstraint(1)  # Here we lock the metal


            # we constrain the all the coordinating atoms to the metal
            sum_of_atoms = []
            for index in self.ligands:
                ligand = self.ligands[index]
                coord_indexes = ligand.ligand_to_metal
                for atom_index in coord_indexes:
                    constraints.AddAtomConstraint(1 + 1 + atom_index + sum(sum_of_atoms))  # The one is to account for open babel indexing starting at 1 and to account for the metal
                sum_of_atoms.append(len(ligand.atomic_props["atoms"]))

            # Set up the force field with the constraints
            forcefield = ob.OBForceField.FindForceField("Uff")
            forcefield.Setup(mol, constraints)
            forcefield.SetConstraints(constraints)

            # Do a 500 steps conjugate gradient minimiazation
            # and save the coordinates to mol.
            # The below function is very slow -> better off just setting i = 200 if you are not debugging
            if debug:
                for i in range(50):
                    forcefield.ConjugateGradients(i)
                    forcefield.GetCoordinates(mol)
                    conv.WriteFile(mol, "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/loop_xyz")
                    self.movie(input_file="loop_xyz", working_directory="/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/")

            elif not debug:
                forcefield.ConjugateGradients(200)
                forcefield.GetCoordinates(mol)
            #
            #
            # UPDATE THE COORDINATES OF THE STK BUILDING BLOCK ISOMER WITH THE NEW COORDINATES
            xyz_string_output = conv.WriteString(mol)
            list_of_nums = re.findall(r"[-+]?(?:\d*\.*\d+)", f"Current Level: {xyz_string_output}")
            num_of_atoms = int(list_of_nums[0])
            del list_of_nums[0]  # we remove the number that corresponds to the number of atoms
            i = 0
            new_position_matrix = []
            for coord in range(num_of_atoms):
                new_position_matrix.append([float(list_of_nums[0 + i]), float(list_of_nums[1 + i]), float(list_of_nums[2 + i])])
                i += 3
            new_position_matrix = np.array(new_position_matrix)
            self.isomer = self.isomer.with_position_matrix(new_position_matrix)
            #print("3 "+str(type(self.isomer)))
            #
            #
            # UPDATE THE COORDINATES OF THE STK BUILDING BLOCK WITH THE NEW COORDINATES
            i = 0
            num_of_lig_atoms = []
            for bb in self.building_blocks.values():
                bb = mercury_remover(bb)
                pos_matrix = bb.get_position_matrix()
                #print(pos_matrix)
                for atom in range(bb.get_num_atoms()):
                    pos_matrix[atom] = new_position_matrix[atom+1+sum(num_of_lig_atoms)]
                num_of_lig_atoms.append(bb.get_num_atoms())
                bb = bb.with_position_matrix(pos_matrix)
                self.building_blocks[i] = bb
                i += 1
        #print("4 "+str(type(self.isomer)))
        return self.isomer, self.building_blocks

import shutil
from typing import Any, Union
from stk import BuildingBlock
import pickle
import numpy as np

from DARTassembler.src.assembly.building_block_utility import rotate_tridentate_bb, rotate_tetradentate_bb, penta_as_tetra, \
    get_optimal_rotation_angle_tridentate, Bidentate_Rotator
from DARTassembler.src.assembly.stk_utils import create_placeholder_Hg_bb
import DARTassembler.src.assembly.stk_extension as stk_e
from DARTassembler.src.assembly.stk_extension import monodentate_coordinating_distance, Bidentate_coordinating_distance
from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand
import ast
from DARTassembler.src.assembly.TransitionMetalComplex import TransitionMetalComplex as TMC
from rdkit import Chem

from DARTassembler.src.ligand_extraction.DataBase import LigandDB

import stk, os
import warnings
from pathlib import Path
import logging

warnings.simplefilter("always")



class PlacementRotation:
    def __init__(self):
        pass

    @staticmethod
    def visualize(input_complex):
        """
        This method allows to visualize in a blocking way during debug but is not essential at all.
        """
        logging.debug("initializing visualization")
        stk.MolWriter().write(input_complex, 'input_complex.mol')
        os.system('obabel .mol input_complex.mol .xyz -O  output_complex.xyz')
        os.system("ase gui output_complex.xyz")
        os.system("rm -f input_complex.mol")
        os.system("rm -f output_complex.xyz")
        logging.debug("visualization complete")

    def touch_file(self, file_path: str):
        # Convert the file path to a Path object
        file_path = Path(file_path)
        # Touch the file by setting its modification time to the current time
        file_path.touch()

    def concatenate_files(self,file1_path, file2_path, output_path):
        """
        Concatenates the contents of two input files into a new output file.

        Args:
            file1_path (str or pathlib.Path): The path to the first input file.
            file2_path (str or pathlib.Path): The path to the second input file.
            output_path (str or pathlib.Path): The path to the output file.
        """
        # Open the output file for writing
        with open(output_path, 'w') as f_out:
            # Open the first input file for reading
            with open(file1_path, 'r') as f_in:
                # Copy the contents of the first input file to the output file
                shutil.copyfileobj(f_in, f_out)
            # Open the second input file for reading
            with open(file2_path, 'r') as f_in:
                # Copy the contents of the second input file to the output file
                shutil.copyfileobj(f_in, f_out)

    def output_controller_(
            self,
            list_of_complexes_wih_isomers: list = None,
            ligands: dict = None,
            metal: str = None,
            metal_ox_state: int = None,
            view_complex: bool = True,
            write_gaussian_input_files: bool = False,
            output_directory: Union[str,Path,None] = None,
            ):

        # in order to check for duplicates we may need to append a list here
        for complex_ in list_of_complexes_wih_isomers:  # We loop through all the created isomers
            if (complex_ is not None) and (complex_ != (None, None)):

                Assembled_complex = TMC.from_stkBB(compl=complex_, ligands=ligands, metal=metal,metal_idx=0, metal_charge=metal_ox_state)
                #
                #
                # 1.
                if view_complex:
                    Assembled_complex.mol.view_3d()
                    input("Press Enter to Continue")

                #
                #
                # 2.
                if write_gaussian_input_files:
                    ###---NAME---###
                    bidentate_ligand = None
                    bidentate_name = None
                    for ligand in ligands.values():
                        if ligand.denticity == 2:
                            bidentate_ligand = ligand
                            bidentate_name = str(ligand.unique_name)
                    name_list = bidentate_name.split("_")
                    name_list[0] = "AuCl2_"
                    name = "".join(map(str, name_list))

                    ###---GAUSSIAN_STRING---###

                    gaussian_string = Assembled_complex.to_gaussian_string(filename=name, num_processors=20, memory=40, charge=0,
                                                                           multiplicity=1, metal_basis_set="lanl2dz",
                                                                           output_directory=output_directory)

                    """gaussian_string = Assembled_complex.to_gaussian_string_Frank(filename=name, num_processors=8, memory=12, charge=0,
                                                                           multiplicity=1, metal_basis_set="SDD",
                                                                           output_directory=output_directory)"""

                    ###---GAUSSIAN_FILE---###
                    os.system(f"mkdir -p {output_directory}")
                    os.system(f"mkdir -p {output_directory}/{name}")
                    com_file = open(f"{output_directory}/{name}/{name}.com", "wt")
                    com_file.write(gaussian_string)
                    com_file.close()

                    ###---BIDENTATE_JSON--###
                    with open(f"{output_directory}/{name}/{name}.pkl", "wb") as outfile:
                        pickle.dump(bidentate_ligand, outfile)
                else:
                    pass
            else:
                logging.debug("!!!Warning!!! -> None type complex detected in list -> skipping to next complex")

    @staticmethod
    def process_single_atom(bb_for_complex):
        # This is here to remove the hydrogen atoms for mono-atomic ligands
        stk_mol = bb_for_complex.to_rdkit_mol()
        edit_mol = Chem.RWMol(stk_mol)
        num_removed_atoms = []
        for atom in stk_mol.GetAtoms():
            pass
            if atom.GetSymbol() == "H":
                #logging.debug(atom.GetIdx())
                edit_mol.RemoveAtom(atom.GetIdx() - len(num_removed_atoms))
                num_removed_atoms.append(1)
            else:
                pass
        output_mol = edit_mol.GetMol()
        output_bb = stk.BuildingBlock.init_from_rdkit_mol(output_mol)
        return output_bb

    @staticmethod
    def format_topologies(top_string):
        #This simply splits the topology and similarity(instruction) lists
        output_list = str(top_string).split("--")
        topology = output_list[0]
        instruction = output_list[1]
        topology = ast.literal_eval(topology)
        instruction = ast.literal_eval(instruction)
        return topology, instruction

    @staticmethod
    def planar_check_(ligands):  # Check if ligands are planar or not
        for ligand in ligands.values():
            if ligand.denticity in [3, 4]:
                return ligand.if_donors_planar(with_metal=True)
            else:
                pass
        return True

    def Process_Monodentate(self, ligand: RCA_Ligand = None, coordinates: list = None):
        # coords specifies the position of the donor atom relative to the metal
        stk_e.Monodentate._ligand_vertex_prototypes[0]._position = np.array(coordinates)
        if len(ligand.atomic_props["atoms"]) == 1:
            if ligand.atomic_props["atoms"] == ["H"]:
                # RDKIT really does not like H- as ligand. This will throw an error in this case.
                raise NotImplementedError('H as a ligand is not supported yet.')

            # This bit of code is a little unstable. Errors are possible here.
            single_atom = ligand.atomic_props['atoms'][0]
            ligand_bb = stk.BuildingBlock(smiles=f"{single_atom}", functional_groups=[stk.SmartsFunctionalGroupFactory(smarts=f"{single_atom}", bonders=(0,), deleters=(), )])
            complex = stk.ConstructedMolecule(
                topology_graph=stk_e.Monodentate(
                    metals=create_placeholder_Hg_bb(),
                    ligands=ligand_bb,
                ),
            )
            bb_for_complex = self.process_single_atom(complex) #stk seems to always generate a building block of a single atom with a hydrogen attached so we need to remove it
            bb_for_complex = stk.BuildingBlock.init_from_molecule(bb_for_complex, functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])
            return bb_for_complex
        else:
            ligand_bb = ligand.to_stk_bb()


            monodentate_topology = stk_e.Monodentate(metals=create_placeholder_Hg_bb(), ligands=ligand_bb)
            bb_for_complex = stk.BuildingBlock.init_from_molecule(stk.ConstructedMolecule(
                topology_graph=monodentate_topology),
                functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())]
            )

            return bb_for_complex

    @staticmethod
    def Process_Bidentate(ligand: RCA_Ligand = None, coordinates: list = None, bidentate_placed: bool = None, top_list: list = None, direction: str = None, build_options: str = None):
        if direction == "Right":
            stk_e.Bidentate_Planar_Right._ligand_vertex_prototypes[0]._position = np.array(coordinates)
            bidentate_topology = stk_e.Bidentate_Planar_Right(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
        elif direction == "Left":
            stk_e.Bidentate_Planar_Left._ligand_vertex_prototypes[0]._position = np.array(coordinates)
            bidentate_topology = stk_e.Bidentate_Planar_Left(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())
        else:
            bidentate_topology = stk_e.Bidentate(metals=create_placeholder_Hg_bb(), ligands=ligand.to_stk_bb())

        complex_bidentate = stk.ConstructedMolecule(topology_graph=bidentate_topology)
        final_bb = stk.BuildingBlock.init_from_molecule(Bidentate_Rotator(ligand_bb=complex_bidentate,
                                                                          ligand=ligand,
                                                                          top_list=top_list,
                                                                          bool_placed=bidentate_placed,
                                                                          build_options=build_options),
                                                        functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]',
                                                                                                            bonders=(0,),
                                                                                                            deleters=(), ), ], )
        return final_bb

    def convert_ligand_to_building_block_for_complex(self, ligands: dict[RCA_Ligand], topology, metal: str = None, build_options: dict = None) -> tuple[dict[int, BuildingBlock], dict[int, Any]]:
        """
        Here we pick and choose our ligands and rotate and place them based on our topology
        """
        if build_options is None:
            build_options = {}

        topology_determining_ligand_planar = self.planar_check_(ligands)  # Check are either the tetra or tri ligands planar

        topology_list = topology
        # This ensures we don't enter the same if statement twice if we have to place a ligand of the same denticity twice
        first_lig0_placed = False
        first_lig1_placed = False
        first_lig2_placed = False

        ligand_buildingblocks = {}
        ligand_denticities = {}
        for i, ligand in enumerate(ligands.values()):
            if ligand.denticity == 0:
                if (topology_list == [4, 1, 0] and (topology_determining_ligand_planar)) or (topology_list == [3, 2, 0]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Top()

                elif (topology_list == [4, 1, 0]) and (not topology_determining_ligand_planar):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Back_Left()
                    first_lig0_placed = True

                elif ((topology_list == [4, 1, 0]) and (not topology_determining_ligand_planar) and (first_lig0_placed)) or (topology_list == [2, 1, 0]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Front_Left()


                elif topology_list == [5, 0]:
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Bottom()

                elif topology_list == [2, 0]:
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Middle_Left()

                else:
                    raise ValueError(f"Your newly created geometry {topology_list}, has not been accounted for in the assembly process (denticity = 0).")
                bb_for_complex = self.Process_Monodentate(ligand=ligand, coordinates=coords)

            elif ligand.denticity == 1:

                if ((((topology_list == [4, 1, 1]) and (not topology_determining_ligand_planar)) or (topology_list == [2, 1, 1])) and (not first_lig1_placed)) or (
                        topology_list == [2, 1, 0]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Back_Left()
                    first_lig1_placed = True

                elif ((((topology_list == [4, 1, 1]) or (topology_list == [4, 1, 0])) and (not topology_determining_ligand_planar)) or (topology_list == [2, 1, 1])) and (
                        first_lig1_placed):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Front_Left()

                elif ((topology_list == [4, 1, 1]) and topology_determining_ligand_planar) and (not first_lig1_placed) or (topology_list == [3, 2, 1]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Top()
                    first_lig1_placed = True

                elif ((topology_list == [4, 1, 1] and first_lig1_placed) or (topology_list == [4, 1, 0])) and topology_determining_ligand_planar or (topology_list == [5, 1]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Bottom()

                else:
                    raise ValueError(f"Your newly created topology {topology_list}, has not been accounted for in the assembly process (denticity = 1).")

                bb_for_complex = self.Process_Monodentate(ligand=ligand, coordinates=coords)


            elif ligand.denticity == 4:
                #
                # If our ligand has denticity of 4 we enter this if statement
                #
                building_block = ligand.to_stk_bb()
                if topology_determining_ligand_planar:
                    # Then some rotation needs to be done
                    building_block = rotate_tetradentate_bb(building_block, ligand_=ligand)
                    tetra_topology_graph = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(), ligands=building_block)
                    bb_for_complex = stk.BuildingBlock.init_from_molecule(stk.ConstructedMolecule(topology_graph=tetra_topology_graph),
                                                                          functional_groups=[
                                                                              stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), )]
                                                                          )
                elif not topology_determining_ligand_planar:
                    raise NotImplementedError("Non-planar tetradentate complexes are not yet supported")
                    # bb_for_complex = nonplanar_tetra_solver(stk_bb=building_block, ligand=ligand)
                    # bb_for_complex = stk.BuildingBlock.init_from_molecule(bb_for_complex, functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg]', bonders=(0,), deleters=(), ), ], )
                else:
                    raise ValueError("Program unable to determine if the tetradentate ligand is planar or not.")


            elif ligand.denticity == 3:
                #
                # If our ligand has denticity of 3 we enter this if statement
                #
                building_block = ligand.to_stk_bb()

                if topology_determining_ligand_planar:
                    building_block = rotate_tridentate_bb(tridentate_bb_=building_block, ligand_=ligand)
                    tridentate_toplogy = stk_e.Tridentate(metals=create_placeholder_Hg_bb(), ligands=building_block)
                    compl_constructed_mol = stk.ConstructedMolecule(topology_graph=tridentate_toplogy)
                    compl_constructed_mol = compl_constructed_mol.with_rotation_about_axis(
                        axis=np.array((0, 0, 1)),
                        angle=float(np.radians(
                            get_optimal_rotation_angle_tridentate(compl_constructed_mol, 10.0, 0.0, 0.0, ligand))),
                        origin=np.array((0, 0, 0))
                    )

                    # Here we essentially shift the tridentate ligand back in the negative x direction by 0.8 A to give a better placement
                    position_matrix = compl_constructed_mol.get_position_matrix()
                    position_matrix[0] = [-0.8, 0, 0]
                    compl_constructed_mol = compl_constructed_mol.with_position_matrix(position_matrix=position_matrix)

                    if topology_list == [3, 2, 0] or [3, 2, 1]:
                        bb_for_complex = stk.BuildingBlock.init_from_molecule(compl_constructed_mol, functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])
                    else:
                        raise ValueError(f"Your newly created topology {topology_list}, has not been accounted for in the assembly process (denticity = 1).")
                else:
                    raise ValueError("The planarity of the tridentate ligand is not clear.")



            elif ligand.denticity == 2:
                #
                # If our ligand has denticity of 2 we enter this if statement
                #
                bidentate_box_choice_instruction = build_options['bidentate_rotator']
                # If this option is `auto`, choose the box shape based on the planarity of the ligand
                logging.debug("rot123:" + str(build_options['bidentate_rotator']))
                if build_options['bidentate_rotator'] == 'auto':
                    if ligand.check_bidentate_planarity():
                        logging.debug("h chosen")
                        bidentate_box_choice_instruction = 'horseshoe'
                    else:
                        bidentate_box_choice_instruction = 'slab'
                        logging.debug("s chosen")

                if topology_list == [3, 2, 0] or topology_list == [3, 2, 1]:
                    coord = Bidentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Bottom()
                    bb_for_complex = self.Process_Bidentate(ligand=ligand, coordinates=coord, direction="Bottom", bidentate_placed=first_lig2_placed, top_list=topology_list, build_options=bidentate_box_choice_instruction)

                elif (topology_list == [2, 2] or [2, 1, 1] or [2, 1, 0] or [2, 0]) and (not first_lig2_placed):
                    coord = Bidentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Right()
                    bb_for_complex = self.Process_Bidentate(ligand=ligand, coordinates=coord, direction="Right", bidentate_placed=first_lig2_placed, top_list=topology_list, build_options=bidentate_box_choice_instruction)
                    first_lig2_placed = True

                elif (topology_list == [2, 2] or [2, 1, 1] or [2, 1, 0]) and first_lig2_placed:
                    coord = Bidentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Left()
                    bb_for_complex = self.Process_Bidentate(ligand=ligand, coordinates=coord, direction="Left", bidentate_placed=first_lig2_placed, top_list=topology_list, build_options=bidentate_box_choice_instruction)

                else:
                    raise ValueError("Geometry not accounted for in the context of bidentate ligands.")





            elif ligand.denticity == 5:
                #
                # If our ligand has denticity of 5 we enter this if statement
                #
                tetra_bb_for_penta, position_index = penta_as_tetra(ligand=ligand)

                tetra_bb_for_penta = rotate_tetradentate_bb(tetra_bb_for_penta, ligand)

                tip_position = list(tetra_bb_for_penta.get_atomic_positions(atom_ids=[int(position_index), ]))

                if float(tip_position[0][2]) > 0:
                    # Additional rotation is required so that the out of plane coordinating atom is facing down (-Z)
                    tetra_bb_for_penta = tetra_bb_for_penta.with_rotation_about_axis(angle=np.radians(180), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)))
                elif float(tip_position[0][2]) < 0:
                    # No rotation is required
                    pass
                else:
                    raise ValueError("!!!Fatal_Error!!! -> Error involving the orientation of the pentadenate ligand-> Exiting Program")

                penta_topology = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(), ligands=tetra_bb_for_penta)

                bb_for_complex = stk.BuildingBlock.init_from_molecule(
                    stk.ConstructedMolecule(topology_graph=penta_topology),
                    functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ])

            else:
                raise ValueError("!!!Fatal Error!!! -> Unknown Ligand Denticity -> Exiting Program")

            #
            # Here we store all the ligand building blocks and their denticities
            #
            ligand_buildingblocks[i] = bb_for_complex
            ligand_denticities[i] = ligand.denticity

        return ligand_buildingblocks, ligand_denticities





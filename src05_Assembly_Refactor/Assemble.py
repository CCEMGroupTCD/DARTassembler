import shutil
from typing import Any
from stk import BuildingBlock
import pickle
import numpy as np
from src05_Assembly_Refactor.building_block_utility import rotate_tridentate_bb, rotate_tetradentate_bb, penta_as_tetra, \
    get_optimal_rotation_angle_tridentate, Bidentate_Rotator, nonplanar_tetra_solver
from src05_Assembly_Refactor.stk_utils import create_placeholder_Hg_bb
import src05_Assembly_Refactor.stk_extension as stk_e
from src05_Assembly_Refactor.stk_extension import monodentate_coordinating_distance, Bidentate_coordinating_distance
from src01.Molecule import RCA_Ligand
import ast
from TransitionMetalComplex import TransitionMetalComplex as TMC
from rdkit import Chem
from src01.DataBase import LigandDB
import stk, os
import warnings
from pathlib import Path

warnings.simplefilter("always")



class PlacementRotation:
    def __init__(self, database: LigandDB = None, store_path: str = Path('..', 'data', 'Assembled_Molecules')):
        self.ligand_dict = database.get_lig_db_in_old_format()
        self.store_path = store_path

    @staticmethod
    def visualize(input_complex):
        """
        This method allows to visualize in a blocking way during debug but is not essential at all.
        """
        print("initializing visualization")
        stk.MolWriter().write(input_complex, 'input_complex.mol')
        os.system('obabel .mol input_complex.mol .xyz -O  output_complex.xyz')
        os.system("ase gui output_complex.xyz")
        os.system("rm -f input_complex.mol")
        os.system("rm -f output_complex.xyz")
        print("visualization complete")

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
            metal_multiplicity: int = None,
            view_complex: bool = True,
            concatonate_xyz_name: str = None,
            write_gaussian_input_files: bool = False,
            output_directory: str =None,
            frames: int = None):

        #todo: in order to check for duplicates we may need to append a list here
        for complex_ in list_of_complexes_wih_isomers:  # We loop through all the created isomers
            if (complex_ is not None) and (complex_ != (None, None)):
                Assembled_complex = TMC(compl=complex_, ligands=ligands, metal=metal, metal_charge=metal_ox_state, spin=metal_multiplicity)
                #
                #
                # 1.
                if view_complex:
                    Assembled_complex.mol.view_3d()
                    input("Press Enter to Continue")
                #
                #
                # 2.
                if concatonate_xyz_name is not None:
                    concat_xyz_filename = Path(output_directory, concatonate_xyz_name)
                    self.touch_file(str(concat_xyz_filename))
                    Assembled_complex.mol.print_to_xyz(str(Path(output_directory, 'tmp_in_xyz.xyz')))  # Print to a temporary file
                    for i in range(frames):
                        pass
                        self.concatenate_files(file1_path=str(concat_xyz_filename), file2_path=str(Path(output_directory, 'tmp_in_xyz.xyz')), output_path=str(Path(output_directory, 'tmp_out_xyz.xyz')))
                        old_path = Path(output_directory, 'tmp_out_xyz.xyz')
                        old_path.rename(str(concat_xyz_filename))
                    Path(output_directory, 'tmp_in_xyz.xyz').unlink()
                #
                #
                # 3.
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
                warnings.warn("!!!Warning!!! -> None type complex detected in list -> skipping to next complex")

    @staticmethod
    def process_single_atom(bb_for_complex):
        # This is here to remove the hydrogen atoms for mono-atomic ligands
        stk_mol = bb_for_complex.to_rdkit_mol()
        edit_mol = Chem.RWMol(stk_mol)
        num_removed_atoms = []
        for atom in stk_mol.GetAtoms():
            pass
            if atom.GetSymbol() == "H":
                #print(atom.GetIdx())
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
            if ligand.denticity == 3:
                return ligand.planar_check()
            else:
                pass
        return True

    def Process_Monodentate(self, ligand: RCA_Ligand = None, coordinates: list = None):
        # coords specifies the position of the donor atom relative to the metal
        stk_e.Monodentate._ligand_vertex_prototypes[0]._position = np.array(coordinates)
        if len(ligand.atomic_props["atoms"]) == 1:
            #TODO: RDKIT really does not like H- as ligand. This will throw an error in this case.  !!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!
            # This bit of code is quite unstable to be brutally honest. errors are quite prone here
            # Todo
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
    def Process_Bidentate(ligand: RCA_Ligand = None, coordinates: list = None, bidentate_placed: bool = None, top_list: list = None, direction: str = None):
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
                                                                          bool_placed=bidentate_placed),
                                                        functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]',
                                                                                                            bonders=(0,),
                                                                                                            deleters=(), ), ], )
        return final_bb

    def convert_ligand_to_building_block_for_complex(self, ligands: dict[RCA_Ligand], topology, metal: str = None) -> tuple[dict[int, BuildingBlock], dict[int, Any]]:
        # Here we pick and choose are ligands and rotate and place them based on our topology
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
                if (topology_list == [4, 1, 0] and (topology_determining_ligand_planar is True)) or (topology_list == [3, 2, 0]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Top()

                elif (topology_list == [4, 1, 0]) and (topology_determining_ligand_planar is False):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Back_Left()
                    first_lig0_placed = True

                elif ((topology_list == [4, 1, 0]) and (topology_determining_ligand_planar is False) and (first_lig0_placed == True)) or (topology_list == [2, 1, 0]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Front_Left()


                elif topology_list == [5, 0]:
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Bottom()

                elif topology_list == [2, 0]:
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Middle_Left()

                else:
                    print(f"!!!Fatal Error!!! -> Your newly created topology {topology_list}, has not been accounted for in the assembly process (denticity = 0) -> Exiting Program ...")
                    raise ValueError
                bb_for_complex = self.Process_Monodentate(ligand=ligand, coordinates=coords)


            elif ligand.denticity == 1:

                if ((((topology_list == [4, 1, 1]) and (topology_determining_ligand_planar is False)) or (topology_list == [2, 1, 1])) and (first_lig1_placed == False)) or (
                        topology_list == [2, 1, 0]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Back_Left()
                    first_lig1_placed = True

                elif ((((topology_list == [4, 1, 1]) or (topology_list == [4, 1, 0])) and (topology_determining_ligand_planar is False)) or (topology_list == [2, 1, 1])) and (
                        first_lig1_placed == True):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Front_Left()

                elif ((topology_list == [4, 1, 1]) and (topology_determining_ligand_planar is True)) and (first_lig1_placed == False) or (topology_list == [3, 2, 1]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Top()
                    first_lig1_placed = True

                elif ((topology_list == [4, 1, 1] and (first_lig1_placed == True)) or (topology_list == [4, 1, 0])) and (topology_determining_ligand_planar is True) or (topology_list == [5, 1]):
                    coords = monodentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Bottom()

                else:
                    print(f"!!!Fatal Error!!! -> Your newly created topology {topology_list}, has not been accounted for in the assembly process (denticity = 1) -> Exiting Program ...")
                    raise ValueError

                bb_for_complex = self.Process_Monodentate(ligand=ligand, coordinates=coords)




            elif ligand.denticity == 4:
                #
                # If our ligand has denticity of 1 we enter this if statement
                #
                building_block = ligand.to_stk_bb()
                if topology_determining_ligand_planar is True:
                    # Then some rotation needs to be done
                    building_block = rotate_tetradentate_bb(building_block, ligand_=ligand)
                    tetra_topology_graph = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(), ligands=building_block)
                    bb_for_complex = stk.BuildingBlock.init_from_molecule(stk.ConstructedMolecule(topology_graph=tetra_topology_graph),
                                                                          functional_groups=[
                                                                              stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), )]
                                                                          )
                elif topology_determining_ligand_planar is False:
                    bb_for_complex = nonplanar_tetra_solver(stk_bb=building_block, ligand=ligand)
                    bb_for_complex = stk.BuildingBlock.init_from_molecule(bb_for_complex, functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg]', bonders=(0,), deleters=(), ), ], )
                else:
                    print("!!!Fatal Error!!! -> Program unable to determine if the tetradentate ligand is planar or not -> Exiting program")
                    raise ValueError


            elif ligand.denticity == 3:
                #
                # If our ligand has denticity of 3 we enter this if statement
                #
                building_block = ligand.to_stk_bb()

                if topology_determining_ligand_planar is True:
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
                        print(f"!!!Fatal Error!!! -> Your newly created topology {topology_list}, has not been accounted for in the assembly process (denticity = 1) -> Exiting Program ...")
                        raise ValueError
                else:
                    print("!!!Fatal Error!!! -> the geometry of the tridentate ligand is unresolved -> Exiting Program")
                    raise ValueError



            elif ligand.denticity == 2:
                #
                # If our ligand has denticity of 2 we enter this if statement
                #
                if topology_list == [3, 2, 0] or topology_list == [3, 2, 1]:
                    coord = Bidentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Bottom()
                    bb_for_complex = self.Process_Bidentate(ligand=ligand, coordinates=coord, direction="Bottom", bidentate_placed=first_lig2_placed, top_list=topology_list)

                elif (topology_list == [2, 2] or [2, 1, 1] or [2, 1, 0] or [2, 0]) and (first_lig2_placed == False):
                    coord = Bidentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Right()
                    bb_for_complex = self.Process_Bidentate(ligand=ligand, coordinates=coord, direction="Right", bidentate_placed=first_lig2_placed, top_list=topology_list)
                    first_lig2_placed = True

                elif (topology_list == [2, 2] or [2, 1, 1] or [2, 1, 0]) and (first_lig2_placed == True):
                    coord = Bidentate_coordinating_distance(metal=metal, ligand=ligand, offset=0).Left()
                    bb_for_complex = self.Process_Bidentate(ligand=ligand, coordinates=coord, direction="Left", bidentate_placed=first_lig2_placed, top_list=topology_list)

                else:
                    print("!!!Fatal_Error!!! -> Topology not accounted for in the context of bidentate ligands.py -> Exiting Program")
                    raise ValueError





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
                    print("!!!Fatal_Error!!! -> Error involving the orientation of the pentadenate ligand-> Exiting Program")
                    raise ValueError

                penta_topology = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(), ligands=tetra_bb_for_penta)

                bb_for_complex = stk.BuildingBlock.init_from_molecule(
                    stk.ConstructedMolecule(topology_graph=penta_topology),
                    functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ])

            else:
                print("!!!Fatal Error!!! -> Unknown Ligand Denticity -> Exiting Program")
                raise ValueError

            #
            # Here we store all the ligand building blocks and their denticities
            #
            ligand_buildingblocks[i] = bb_for_complex
            ligand_denticities[i] = ligand.denticity

        return ligand_buildingblocks, ligand_denticities




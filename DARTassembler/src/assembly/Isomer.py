import DARTassembler.src.assembly.stk_extension as stk_e
import numpy as np
import stk
from DARTassembler.src.assembly.building_block_utility import get_energy_stk, mercury_remover
import logging

class BuildIsomers:
    def __init__(self, topology, building_blocks_list, metal_input, charge_input, denticity_list, return_all_isomers, opt_choice: bool, ligand_list):
        self.topology = topology
        self.building_blocks_list = building_blocks_list
        self.metal_input = metal_input
        self.charge_input = charge_input
        self.denticity_list = denticity_list
        self.return_all_isomers = return_all_isomers
        self.opt_choice = opt_choice    # todo: not used atm, what should happen here?
        self.ligand_list = ligand_list  # todo: not used atm, what should happen here?
        self.ALL_building_blocks = []
        self.ALL_Assembled_complexes = []

    @staticmethod
    def create_metal_building_block(metal, charge):  # Creates the metal building block
        # build the metal block with the new metal atom
        # This may not work for neutral metals
        if str(charge) != "0":
            smiles_str = f"[{metal}{charge}]"
        elif str(charge) == "0":
            smiles_str = f"[{metal}]"
        else:
            smiles_str = None
        stk_metal_func = getattr(stk, metal)
        functional_groups = (stk.SingleAtom(stk_metal_func(0, charge=charge)) for i in range(6))
        final_metal_bb = stk.BuildingBlock(smiles=smiles_str,
                                           functional_groups=functional_groups,
                                           position_matrix=np.ndarray([0, 0, 0])
                                           )
        return final_metal_bb

    def Generate(self):
        #
        #
        # 1.
        if (str(self.topology) == '[2, 1, 0]--[1, 2, 3]') or (str(self.topology) == "[2, 1, 1]--[1, 2, 2]") or (str(self.topology) == "[2, 1, 1]--[1, 2, 3]") or (str(self.topology) == "[2, 2]--[1, 1]") or (
                str(self.topology) == "[2, 2]--[1, 2]"):
            # Expect only one isomer from this if statement (so a list of two constructed molecules). Only a single bidentate needs to be flipped 180
            bidentate_already_rotated = False  # This is to ensure we don't enter the same if statement twice
            building_blocks_rotated = {}  # This will contain all the ligands for the isomer of our complex
            for j in range(len(self.building_blocks_list)):
                # We iterate through all the complexes in our input complex
                if (self.denticity_list[j] == 2) and (bidentate_already_rotated == False):
                    # If we come across a bidentate ligand, and it's the first time we have come across one then ...
                    building_blocks_rotated[j] = self.building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)), )
                    bidentate_already_rotated = True
                else:
                    # Otherwise we just copy over the ligand from the previous complex
                    building_blocks_rotated[j] = self.building_blocks_list[j]
            self.ALL_building_blocks.extend([self.building_blocks_list, building_blocks_rotated])
        #
        #
        # 2.
        elif (str(self.topology) == '[3, 2, 0]--[1, 2, 3]') or (str(self.topology) == '[3, 2, 1]--[1, 2, 3]'):
            # These are the empty dictionaries that will contain all the building blocks for each isomer
            # This is by far the most 'involved' isomer generation processes
            building_blocks_rotated_bi = {}
            building_blocks_rotated_tri = {}
            building_blocks_rotated_bi_and_tri = {}
            for j in range(len(self.building_blocks_list)):
                if self.denticity_list[j] == 2:
                    building_blocks_rotated_bi[j] = self.building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((-1, 0, -1)), origin=np.array((0, 0, 0)), )
                    building_blocks_rotated_bi_and_tri[j] = self.building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((-1, 0, -1)), origin=np.array((0, 0, 0)), )
                    building_blocks_rotated_tri[j] = self.building_blocks_list[j]
                    pass
                elif self.denticity_list[j] == 3:
                    building_blocks_rotated_tri[j] = self.building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)), )
                    building_blocks_rotated_bi_and_tri[j] = self.building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)), )
                    building_blocks_rotated_bi[j] = self.building_blocks_list[j]
                else:
                    building_blocks_rotated_bi[j] = self.building_blocks_list[j]
                    building_blocks_rotated_tri[j] = self.building_blocks_list[j]
                    building_blocks_rotated_bi_and_tri[j] = self.building_blocks_list[j]
                    pass
            self.ALL_building_blocks.extend([self.building_blocks_list, building_blocks_rotated_bi, building_blocks_rotated_tri, building_blocks_rotated_bi_and_tri])

        #
        #
        # 3
        elif (str(self.topology) == "[4, 1, 1]--[1, 2, 3]") or (str(self.topology) == "[4, 1, 0]--[1, 2, 3]") or (str(self.topology) == "[4, 1, 1]--[1, 2, 2]"):
            #Here we are just flipping the tetradentate ligand
            building_blocks_rotated_tetra = {}
            for j in range(len(self.building_blocks_list)):
                if self.denticity_list[j] == 4:

                    # logging.debug(building_blocks_list[j].get_position_matrix())
                    building_blocks_rotated_tetra[j] = self.building_blocks_list[j].with_rotation_about_axis(angle=180.0 * (np.pi / 180.0), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)), )
                else:
                    # logging.debug(building_blocks_list[j].get_position_matrix())
                    # logging.debug(list(building_blocks_list[j].get_atoms()))

                    building_blocks_rotated_tetra[j] = self.building_blocks_list[j]
            self.ALL_building_blocks.extend([self.building_blocks_list, building_blocks_rotated_tetra])
        #
        #
        # 4
        else:
            self.ALL_building_blocks.append(self.building_blocks_list)


        for building_blocks_list in self.ALL_building_blocks:
            if len(building_blocks_list) == 2:
                # If there are only two ligands that comprise our complex then ...
                # This is our input complex
                complex = stk_e.complex_topology_two(metals=self.create_metal_building_block(self.metal_input, self.charge_input),
                                                     ligands={building_block: (i,) for i, building_block in building_blocks_list.items()})
            elif len(building_blocks_list) == 3:
                # If there are 3 ligands in our complex then ...
                complex = stk_e.complex_topology_three(metals=self.create_metal_building_block(self.metal_input, self.charge_input),
                                                       ligands={building_block: (i,) for i, building_block in building_blocks_list.items()})
            else:
                raise ValueError
            complex_built = stk.ConstructedMolecule(topology_graph=complex)
            complex_built = mercury_remover(complex_built)
            self.ALL_Assembled_complexes.append(complex_built)

        if self.return_all_isomers == "Generate All":
            return self.ALL_Assembled_complexes, self.ALL_building_blocks

        elif self.return_all_isomers == "Generate Lowest Energy":
            energy_list = []
            for complex in self.ALL_Assembled_complexes:
                energy_list.append(get_energy_stk(complex))
            index = energy_list.index(min(energy_list))
            return [self.ALL_Assembled_complexes[index]], [self.ALL_building_blocks[index]]
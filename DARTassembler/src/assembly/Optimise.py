from DARTassembler.src.assembly.building_block_utility import mercury_remover
import numpy as np
import re
from DARTassembler.src.ligand_extraction.utilities_Molecule import get_concatenated_xyz_string_from_coordinates
from DARTassembler.src.assembly.forcefields import ForceField
import logging

class OPTIMISE:
    def __init__(self, isomer, ligand_list, building_blocks, nsteps: int=50):
        self.isomer = isomer
        self.ligands = ligand_list
        self.building_blocks = building_blocks
        self.nsteps = nsteps

    def Optimise_STK_Constructed_Molecule(self, return_ff_movie: bool = False):
        # logging.debug("1 " + str(type(self.isomer)))

        logging.debug("Beginning Optimisation")
        # logging.debug("2 " + str(type(self.isomer)))
        # todo: Return the rotated stk building blocks as well, they may still be of use to someone
        # REFERENCE: https://github.com/hjkgrp/molSimplify/blob/c3776d0309b5757d5d593e42db411a251b558e59/molSimplify/Scripts/structgen.py#L658
        # REFERENCE: https://gist.github.com/andersx/7784817

        # stk to xyz string
        complex_mol = self.isomer.to_rdkit_mol()
        num_atoms = complex_mol.GetNumAtoms()

        # we constrain the all the coordinating atoms to the metal
        fixed_atom_indices = [0] # Lock the metal. For `fixed_atom_indices`, indexing is 0-based.
        sum_of_atoms = []
        for ligand in self.ligands.values():
            coord_indexes = ligand.ligand_to_metal
            for atom_index in coord_indexes:
                fixed_atom_idx = 1 + atom_index + sum(sum_of_atoms) # 0-based indexing, the +1 comes from somewhere else.
                fixed_atom_indices.append(fixed_atom_idx)  # The one is to account for open babel indexing starting at 1 and to account for the metal
                assert fixed_atom_idx < num_atoms
            # this is so we don't take into account any mercury that might be in the atomic props (really only an issue for tetradentate non-planar ligands as they make use of the add atom function)
            sum_of_atoms.append(len([i for i in ligand.atomic_props["atoms"] if i != "Hg"]))

        xyz_string_output, optimized_coords, optimized_elements = ForceField().optimize(complex_mol, fixed_atom_indices, self.nsteps)

        # UPDATE THE COORDINATES OF THE STK BUILDING BLOCK ISOMER WITH THE NEW COORDINATES
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
        # logging.debug("3 "+str(type(self.isomer)))
        #
        #
        # UPDATE THE COORDINATES OF THE STK BUILDING BLOCK WITH THE NEW COORDINATES
        i = 0
        num_of_lig_atoms = []
        for bb in self.building_blocks.values():
            bb = mercury_remover(bb)
            pos_matrix = bb.get_position_matrix()
            # logging.debug(pos_matrix)
            for atom in range(bb.get_num_atoms()):
                pos_matrix[atom] = new_position_matrix[atom + 1 + sum(num_of_lig_atoms)]
            num_of_lig_atoms.append(bb.get_num_atoms())
            bb = bb.with_position_matrix(pos_matrix)
            self.building_blocks[i] = bb
            i += 1

        if return_ff_movie:
            xzy_string = get_concatenated_xyz_string_from_coordinates(optimized_coords, optimized_elements)
            return self.isomer, self.building_blocks, xzy_string
        else:
            return self.isomer, self.building_blocks

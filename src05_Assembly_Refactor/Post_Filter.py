from building_block_utility import mercury_remover
from mendeleev import element
import numpy as np
import warnings
import os, stk



class PostFilter:
    def __init__(self, isomer, metal_centre, metal_offset, building_blocks, instruction):
        self.building_blocks = building_blocks
        self.metal = metal_centre
        self.isomer = isomer
        self.metal_offset = metal_offset
        self.failed_isomers = []
        self.threshold = 0.3  # Angstrom todo: needs to be tuned
        self.instruction = instruction

    @staticmethod
    def visualize(input_complex):
        print("initializing visualization")
        stk.MolWriter().write(input_complex, 'input_complex.mol')
        os.system('obabel .mol input_complex.mol .xyz -O  output_complex.xyz')
        os.system("ase gui output_complex.xyz")
        os.system("rm -f input_complex.mol")
        os.system("rm -f output_complex.xyz")
        print("visualization complete")

    def closest_distance(self):
        # This function will detect and filter out complexes that have clashing between atoms in different ligands.py but it WILL NOT detect clashing between atoms of the same
        # ligand or the metal centre
        # todo: detect clashing between atoms and the metal centre without taking into account the functional groups maybe an idea for the future
        for keys_1, values_1 in self.building_blocks.items():  # loop through the building blocks within each isomer
            values_1 = mercury_remover(values_1)  # Eliminate the temporary Mercury
            for keys_2, values_2 in self.building_blocks.items():  # loop through the building blocks for within each isomer
                values_2 = mercury_remover(values_2)  # Eliminate the temporary Mercury
                if keys_1 == keys_2:  # Don't compare anything if they are the same ligand
                    pass
                elif keys_1 != keys_2:  # Compare distance if the ligands.py are not the same ligand
                    for point_1, atom_1 in zip(values_1.get_position_matrix(), values_1.get_atoms()):  # we loop through all the positions of the atoms
                        atom_1_type = [str(atom_1).split("(")][0][0]
                        cov_1 = element(atom_1_type).covalent_radius / 100.0
                        cov_metal = element(self.metal).covalent_radius / 100.0
                        metal_position = [0, 0, 0]
                        distance_metal = np.linalg.norm(point_1 - metal_position)
                        if distance_metal < (cov_1 + cov_metal + self.metal_offset):
                            warnings.warn("!!!Warning!!! -> Pre-optimisation filter failed (1)-> Returning None")
                            #self.visualize(self.isomer)
                            return None
                        for point_2, atom_2 in zip(values_2.get_position_matrix(), values_2.get_atoms()):  # we loop through all the positions of the atoms
                            atom_2_type = [str(atom_2).split("(")][0][0]
                            cov_2 = element(atom_2_type).covalent_radius / 100.0
                            distance = np.linalg.norm(point_1 - point_2)  # Calculate distance
                            if distance < (cov_1 + cov_2 + self.metal_offset):  # This function shouldn't take into account ligand metal distances
                                warnings.warn("!!!Warning!!! -> Pre-optimisation filter failed (2)-> Returning None")
                                #self.visualize(self.isomer)
                                return None
                            else:
                                pass
        return self.isomer

    def post_optimisation_filter(self):
        if (self.building_blocks is None) and (self.isomer is None):
            warnings.warn("!!!Warning!!! -> None detect in post-optimisation filter -> Returning None")
            return None
        elif self.instruction == "False":
            return self.isomer
        else:
            pass
        for keys_1, values_1 in self.building_blocks.items():
            for bond in list(values_1.get_bonds()):
                atom_1_id = bond.get_atom1().get_id()
                atom_2_id = bond.get_atom2().get_id()
                atom_1_pos = list(values_1.get_atomic_positions(atom_1_id))
                atom_2_pos = list(values_1.get_atomic_positions(atom_2_id))
                atom_1_AN = bond.get_atom1().get_atomic_number()
                atom_2_AN = bond.get_atom2().get_atomic_number()
                distance = np.linalg.norm(atom_1_pos[0] - atom_2_pos[0])
                # if distance > ((Element.from_Z(atom_1_AN).atomic_radius + Element.from_Z(atom_2_AN).atomic_radius) + self.threshold):
                if distance > (((element(atom_1_AN).covalent_radius / 100) + (element(atom_2_AN).covalent_radius / 100)) + self.threshold):
                    warnings.warn("!!!Warning!!! -> Post-optimisation filter failed -> Returning None")
                    return None
                else:
                    pass
        return self.isomer

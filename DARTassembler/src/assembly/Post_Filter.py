from DARTassembler.src.assembly.building_block_utility import mercury_remover
import numpy as np
import warnings
import os
import stk
from DARTassembler.src.constants.Periodic_Table import DART_Element
import logging

elem_cov_radii = {'H': 0.32, 'He': 0.46, 'Li': 1.33, 'Be': 1.02, 'B': 0.85, 'C': 0.75, 'N': 0.71, 'O': 0.63, 'F': 0.64, 'Ne': 0.67, 'Na': 1.55, 'Mg': 1.39, 'Al': 1.26, 'Si': 1.16, 'P': 1.11,
                  'S': 1.03, 'Cl': 0.99, 'Ar': 0.96, 'K': 1.96, 'Ca': 1.71, 'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22, 'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11, 'Ni': 1.1, 'Cu': 1.12, 'Zn': 1.18,
                  'Ga': 1.24, 'Ge': 1.21, 'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17, 'Rb': 2.1, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25,
                  'Pd': 1.2, 'Ag': 1.28, 'Cd': 1.36, 'In': 1.42, 'Sn': 1.4, 'Sb': 1.4, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31, 'Cs': 2.32, 'Ba': 1.96, 'La': 1.8, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74,
                  'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68, 'Gd': 1.69, 'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.7, 'Lu': 1.62, 'Hf': 1.52, 'Ta': 1.46, 'W': 1.37, 'Re': 1.31,
                  'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33, 'Tl': 1.44, 'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01, 'Ac': 1.86, 'Th': 1.75,
                  'Pa': 1.69, 'U': 1.7, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66, 'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65, 'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57, 'Db': 1.49,
                  'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34, 'Mt': 1.29, 'Ds': 1.28, 'Rg': 1.21, 'Cn': 1.22, 'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75, 'Ts': 1.65, 'Og': 1.57}



class PostFilter:
    def __init__(self, isomer, metal_centre, metal_offset, ligand_atom_offset, building_blocks):
        self.building_blocks = building_blocks
        self.metal = metal_centre
        self.isomer = isomer
        self.metal_offset = metal_offset
        self.ligand_atom_offset = ligand_atom_offset
        self.failed_isomers = []
        self.threshold = 0.3  # Angstrom todo: needs to be tuned
        logging.debug("Post Filter Class initialized Succesfully")

    @staticmethod
    def visualize(input_complex):
        logging.debug("initializing visualization")
        stk.MolWriter().write(input_complex, 'input_complex.mol')
        os.system('obabel .mol input_complex.mol .xyz -O  output_complex.xyz')
        os.system("ase gui output_complex.xyz")
        os.system("rm -f input_complex.mol")
        os.system("rm -f output_complex.xyz")
        logging.debug("visualization complete")

    def closest_distance(self) -> bool:
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
                        # cov_1 = element(atom_1_type).covalent_radius / 100.0
                        cov_1 = elem_cov_radii[atom_1_type]
                        # cov_metal = element(self.metal).covalent_radius / 100.0
                        cov_metal = elem_cov_radii[self.metal]
                        metal_position = [0, 0, 0]
                        distance_metal = np.linalg.norm(point_1 - metal_position)
                        if distance_metal < (cov_1 + cov_metal + self.metal_offset):
                            logging.debug("!!!Warning!!! -> Pre-optimisation filter failed (1)-> Returning None")
                            # self.visualize(self.isomer)
                            return False
                        for point_2, atom_2 in zip(values_2.get_position_matrix(), values_2.get_atoms()):  # we loop through all the positions of the atoms
                            atom_2_type = [str(atom_2).split("(")][0][0]
                            # cov_2 = element(atom_2_type).covalent_radius / 100.0
                            cov_2 = elem_cov_radii[atom_2_type]
                            distance = np.linalg.norm(point_1 - point_2)  # Calculate distance
                            if distance < (cov_1 + cov_2 + self.ligand_atom_offset):  # This function shouldn't take into account ligand metal distances
                                logging.debug("!!!Warning!!! -> Pre-optimisation filter failed (2)-> Returning None")
                                # self.visualize(self.isomer)
                                return False
                            else:
                                pass
        return True

    def post_optimisation_filter(self) -> bool:
        for keys_1, values_1 in self.building_blocks.items():
            for bond in list(values_1.get_bonds()):
                atom_1_id = bond.get_atom1().get_id()
                atom_2_id = bond.get_atom2().get_id()
                atom_1_pos = list(values_1.get_atomic_positions(atom_1_id))
                atom_2_pos = list(values_1.get_atomic_positions(atom_2_id))
                atom_1_AN = DART_Element(bond.get_atom1()._atomic_number).symbol
                atom_2_AN = DART_Element(bond.get_atom2()._atomic_number).symbol
                distance = np.linalg.norm(atom_1_pos[0] - atom_2_pos[0])

                if distance > (((elem_cov_radii[atom_1_AN]) + (elem_cov_radii[atom_2_AN])) + self.threshold):
                    logging.debug("!!!Warning!!! -> Post-optimisation filter failed -> Returning None")
                    return False
                else:
                    pass

        return True

from rdkit.Chem import rdmolfiles
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from openbabel import openbabel as ob
import stk
import os
from openbabel import pybel
import itertools
from DARTassembler.src.assembly.stk_extension import *
from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand
from copy import deepcopy


class RandomComplexAssembler:
    """
    is kind of the hub to handle the random assembly.
    Stores the ligand database and stores the configuration for the random assembly
    """

    # def __init__(self, database: LigandDB, store_path: str = "../data/Assembled_Molecules"): If somethig breaks in the assembly it may be because I have replaced this line with the one below
    def __init__(self, database: LigandDB, store_path: str = "../data/Assembled_Molecules"):
        self.ligand_dict = database.get_lig_db_in_old_format()
        self.store_path = store_path

    def ligand_choice(self):
        pass


F = LigandDB.from_json(json_="/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/data/Filtered_Jsons/filteredLigDB_03032.json",
                       type_="Ligand")  # We initiate the database in the RandomComplexAssembler
RCA = RandomComplexAssembler(database=F)


def ligand_to_mol(ligand: RCA_Ligand):
    xyz_str = ligand.get_xyz_file_format_string()
    mol_b = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("xyz", "mol")
    conv.ReadString(mol_b, xyz_str)
    string = conv.WriteString(mol_b)
    return string


def get_energy(molecule):
    mol_string = ligand_to_mol(molecule)
    mol = pybel.readstring("mol", str(mol_string))
    obmol = mol.OBMol
    ff = ob._openbabel.OBForceField_FindType("uff")
    assert (ff.Setup(obmol))
    kj_to_kcal = 1.0 / 4.184
    ff.SetCoordinates(mol.OBMol)
    # false = don't calculate gradients
    uffE = ff.Energy(False) * kj_to_kcal
    return uffE


def visualize(input_complex):
    print("initializing visualization")
    stk.MolWriter().write(input_complex, 'input_complex.mol')
    os.system('obabel .mol input_complex.mol .xyz -O  output_complex.xyz')
    os.system("ase gui output_complex.xyz")
    os.system("rm -f input_complex.mol")
    os.system("rm -f output_complex.xyz")
    print("visualization complete")


# Notes: we need a function that takes a ligand and an stk building block as an input

def nonplanar_tetra_solver_modified(stk_bb, ligand):
    all_Energies, all_midpoints = [], []
    combo_list = list(itertools.combinations([0, 1, 2, 3], 2))
    assembly_dictionary = ligand.get_assembly_dict()
    location_list = assembly_dictionary["index"]  # location of every coord atom in file
    type_list = assembly_dictionary["type"]
    for combo in combo_list:  # this loop iterates through every combination of functional Groups
        ligand_copy = deepcopy(ligand)
        func_1_location = location_list[combo[0]]  # first atom location
        func_2_location = location_list[combo[1]]  # first atom type

        # we needed to initialise the building block solely to get the positions of the functional groups
        positions = list(stk_bb.get_atomic_positions(atom_ids=[int(func_1_location), int(func_2_location), ], ))
        x1 = positions[0][0]
        y1 = positions[0][1]
        z1 = positions[0][2]
        x2 = positions[1][0]
        y2 = positions[1][1]
        z2 = positions[1][2]
        x_mid = (x1 + x2) / 2  # Here we get the components of the mid points between coord groups
        y_mid = (y1 + y2) / 2
        z_mid = (z1 + z2) / 2
        mid_points = [x_mid, y_mid, z_mid]
        # The add_atom functionality does not update the coordinates attribute correctly in the molecule class

        ligand_copy.add_atom(symbol="Hg", coordinates=[x_mid, y_mid, z_mid])  # add dummy metal
        Energy = get_energy(ligand_copy)  # get energy
        # lig.remove_last_element_in_xyz()  # get rid of dummy metal
        all_Energies.append(Energy)  # list of all the energies of pacing dummy metal at each midpoint
        all_midpoints.append(mid_points)
        del ligand_copy
        print("iteration done")
    minimum_energy = min(all_Energies)
    minimum_energy_index = all_Energies.index(minimum_energy)
    ligand.add_atom(symbol="Hg", coordinates=[all_midpoints[minimum_energy_index][0], all_midpoints[minimum_energy_index][1],
                                              all_midpoints[minimum_energy_index][2]])  # paces Hg at midpoint with the smallest energy
    tetra_bb_2 = stk.BuildingBlock.init_from_rdkit_mol(rdmolfiles.MolFromMolBlock(ligand_to_mol(ligand=ligand), removeHs=False, sanitize=False, strictParsing=False), functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts="[Hg]", bonders=(0,), deleters=(), ), ], )

    tetra_bb_2 = tetra_bb_2.with_displacement(np.array(((-1) * all_midpoints[minimum_energy_index][0],
                                                        (-1) * all_midpoints[minimum_energy_index][1],
                                                        (-1) * all_midpoints[minimum_energy_index][2])))

    # This block of code allows me to remove all the coordinating groups that were previously used to coord the
    # temporary atom
    # Note location refers to the location within the file
    location_list = list(location_list)
    index1 = combo_list[minimum_energy_index][0]
    index2 = combo_list[minimum_energy_index][1]
    value1 = location_list[index1]
    value2 = location_list[index2]
    location_list.remove(value1)
    location_list.remove(value2)

    # this block removes the coord atoms used to coordinate to the temporary metal
    value1 = type_list[index1]
    value2 = type_list[index2]
    type_list.remove(value1)
    type_list.remove(value2)

    position_of_Hg_in_mol = ligand.atomic_props["atoms"].index("Hg")

    # The following Block of code ensures that the remaining coordinating groups exist in the xy plane
    complex_tetradentate = tetra_bb_2.with_rotation_to_minimize_angle(start=tetra_bb_2.get_plane_normal(
        atom_ids=[int(location_list[0]), int(location_list[1]),
                  int(position_of_Hg_in_mol), ]), target=np.array((0, 0, 1)), axis=np.array((0, 1, 0)),
        origin=np.array((0, 0, 0)), )

    complex_tetradentate = complex_tetradentate.with_rotation_to_minimize_angle(
        start=complex_tetradentate.get_plane_normal(
            atom_ids=[int(location_list[0]), int(location_list[0]),
                      int(position_of_Hg_in_mol), ]), target=np.array((0, 0, 1)), axis=np.array((1, 0, 0)),
        origin=np.array((0, 0, 0)), )

    for _ in np.arange(0, 361, 0.5):
        complex_tetradentate = complex_tetradentate.with_rotation_about_axis(angle=0.5 * (np.pi / 180.0),
                                                                             axis=np.array((0, 0, 1)),
                                                                             origin=np.array((0, 0, 0)),
                                                                             )

        position_of_coord_atom_1 = complex_tetradentate.get_atomic_positions(atom_ids=[int(location_list[0]), ])
        position_of_coord_atom_2 = complex_tetradentate.get_atomic_positions(atom_ids=[int(location_list[1]), ])

        # Here we are minimising the difference in distances between the two functional groups and a point far away on the -x axis
        distance1 = np.linalg.norm(list(position_of_coord_atom_1) - np.array((-10.0, 0, 0)))
        distance2 = np.linalg.norm(list(position_of_coord_atom_2) - np.array((-10.0, 0, 0)))

        mean = (distance1 + distance2) / 2.0
        deviation1 = (distance1 - mean) ** 2
        deviation2 = (distance2 - mean) ** 2
        variance = (deviation1 + deviation2) / 2.0
        if (variance < 0.001) and (distance1 > 10.0) and (distance2 > 10.0):
            return complex_tetradentate
        else:
            pass

def nonplanar_tetra_solver_modified_2(stk_bb, ligand):
    all_Energies, all_midpoints = [], []
    combo_list = list(itertools.combinations([0, 1, 2, 3], 2))
    assembly_dictionary = ligand.get_assembly_dict()
    location_list = assembly_dictionary["index"]  # location of every coord atom in file
    type_list = assembly_dictionary["type"]
    for combo in combo_list:  # this loop iterates through every combination of functional Groups
        ligand_copy = deepcopy(ligand)
        func_1_location = location_list[combo[0]]  # first atom loction
        func_2_location = location_list[combo[1]]  # first atom type

        # we needed to initialise the building block solely to get the positons of the functional groups
        positions = list(stk_bb.get_atomic_positions(atom_ids=[int(func_1_location), int(func_2_location), ], ))
        x1 = positions[0][0]
        y1 = positions[0][1]
        z1 = positions[0][2]
        x2 = positions[1][0]
        y2 = positions[1][1]
        z2 = positions[1][2]
        x_mid = (x1 + x2) / 2  # Here we get the components of the mid points between coord groups
        y_mid = (y1 + y2) / 2
        z_mid = (z1 + z2) / 2
        mid_points = [x_mid, y_mid, z_mid]
        # The add_atom functionality does not update the coordinates attribute correctly in the molecule class

        ligand_copy.add_atom(symbol="Hg", coordinates=[x_mid, y_mid, z_mid])  # add dummy metal
        Energy = get_energy(ligand_copy)  # get energy
        # lig.remove_last_element_in_xyz()  # get rid of dummy metal
        all_Energies.append(Energy)  # list of all the energies of pacing dummy metal at each midpoint
        all_midpoints.append(mid_points)
        del ligand_copy
        print("iteration done")
    minimum_energy = min(all_Energies)
    minimum_energy_index = all_Energies.index(minimum_energy)
    ligand.add_atom(symbol="Hg", coordinates=[all_midpoints[minimum_energy_index][0], all_midpoints[minimum_energy_index][1],
                                              all_midpoints[minimum_energy_index][2]])  # paces Hg at midpoint with the smallest energy
    tetra_bb_2 = stk.BuildingBlock.init_from_rdkit_mol(rdmolfiles.MolFromMolBlock(ligand_to_mol(ligand=ligand), removeHs=False, sanitize=False, strictParsing=False), functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts="[Hg]", bonders=(0,), deleters=(), ), ], )

    tetra_bb_2 = tetra_bb_2.with_displacement(np.array(((-1) * all_midpoints[minimum_energy_index][0],
                                                        (-1) * all_midpoints[minimum_energy_index][1],
                                                        (-1) * all_midpoints[minimum_energy_index][2])))

    # This block of code allows me to remove all the coordinating groups that were previously used to coord the
    # temporary atom
    # Note location refers to the location within the file
    location_list = list(location_list)
    index1 = combo_list[minimum_energy_index][0]
    index2 = combo_list[minimum_energy_index][1]
    value1 = location_list[index1]
    value2 = location_list[index2]
    location_list.remove(value1)
    location_list.remove(value2)

    # this block removes the coord atoms used to coordinate to the temporary metal
    value1 = type_list[index1]
    value2 = type_list[index2]
    type_list.remove(value1)
    type_list.remove(value2)

    position_of_Hg_in_mol = ligand.atomic_props["atoms"].index("Hg")

    # The following Block of code ensures that the remaining coordinating groups exist in the xy plane
    complex_tetradentate = tetra_bb_2.with_rotation_to_minimize_angle(start=tetra_bb_2.get_plane_normal(
        atom_ids=[int(location_list[0]), int(location_list[1]),
                  int(position_of_Hg_in_mol), ]), target=np.array((0, 0, 1)), axis=np.array((0, 1, 0)),
        origin=np.array((0, 0, 0)), )

    complex_tetradentate = complex_tetradentate.with_rotation_to_minimize_angle(
        start=complex_tetradentate.get_plane_normal(
            atom_ids=[int(location_list[0]), int(location_list[0]),
                      int(position_of_Hg_in_mol), ]), target=np.array((0, 0, 1)), axis=np.array((1, 0, 0)),
        origin=np.array((0, 0, 0)), )

    for _ in np.arange(0, 361, 0.5):
        complex_tetradentate = complex_tetradentate.with_rotation_about_axis(angle=0.5 * (np.pi / 180.0),
                                                                             axis=np.array((0, 0, 1)),
                                                                             origin=np.array((0, 0, 0)),
                                                                             )

        position_of_coord_atom_1 = complex_tetradentate.get_atomic_positions(atom_ids=[int(location_list[0]), ])
        position_of_coord_atom_2 = complex_tetradentate.get_atomic_positions(atom_ids=[int(location_list[1]), ])

        # Here we are minimising the difference in distances between the two functional groups and a point far away on the -x axis
        distance1 = np.linalg.norm(list(position_of_coord_atom_1) - np.array((-10.0, 0, 0)))
        distance2 = np.linalg.norm(list(position_of_coord_atom_2) - np.array((-10.0, 0, 0)))

        mean = (distance1 + distance2) / 2.0
        deviation1 = (distance1 - mean) ** 2
        deviation2 = (distance2 - mean) ** 2
        variance = (deviation1 + deviation2) / 2.0
        if (variance < 0.001) and (distance1 > 10.0) and (distance2 > 10.0):
            return complex_tetradentate
        else:
            pass


if __name__ == "__main__":

    denticity = 1
    for ligand in RCA.ligand_dict[denticity]:
        print(ligand.atomic_props["atoms"])

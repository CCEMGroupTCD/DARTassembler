import stk
from openbabel import pybel
from src03_Assembly.stk_extension import *
import os
from ase import io
from src01.Molecule import RCA_Molecule, RCA_Ligand
import rdkit
from src02_Pre_Ass_Filtering.constants import get_boxes, intensity, sharpness
from openbabel import openbabel as ob
import itertools


names_dict = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}


def get_topology_string(top):
    str_ = ""
    for el in top:
        str_ += f"{names_dict[el]}_"

    return f"{str_}assembly"



def optimize(input_file, option: str):  # option will be stk_to_xyz or stk_to_stk or file_to_file
    if option == "xyz_to_xyz":
        print("optimizing xyz_to_xyz")
        file_extension = input_file[input_file.find(".") + 1:].split()[0]
        mol = next(pybel.readfile(format=str(file_extension), filename=str(input_file)))
        mol.localopt(forcefield='uff', steps=10_000)
        mol.write("xyz", f"../tmp/{input_file}_optimization_output_if.xyz")

    elif option == "stk_to_xyz":
        print("optimizing stk_to_xyz")
        stk.MolWriter().write(input_file, '../tmp/input_complex.mol')
        mol = next(pybel.readfile("mol", "../tmp/input_complex.mol"))
        mol.localopt(forcefield='uff', steps=10_000)
        mol.write("xyz", "../tmp/Optimizer_output_stk_to_xyz.xyz", overwrite="True")
        os.system("rm -f ../tmp/input_complex.mol")

    elif option == "stk_to_mol":
        print("optimizing stk_to_mol")
        stk.MolWriter().write(input_file, '../tmp/input_complex.mol')
        mol = next(pybel.readfile("mol", "../tmp/input_complex.mol"))
        mol.localopt(forcefield='uff', steps=10_000)
        mol.write("mol", "../tmp/Optimizer_output_stk_to_mol.mol", overwrite="True")
        os.system("rm -f ../tmp/input_complex.mol")

    elif option == "stk_to_stk":
        print("optimizing stk_to_stk")
        stk.MolWriter().write(input_file, '../tmp/input_complex.mol')
        mol = next(pybel.readfile("mol", "../tmp/input_complex.mol"))
        mol.localopt(forcefield='uff', steps=3)
        mol.write("mol", "../tmp/Optimizer_output_else.mol", overwrite="True")
        os.system("rm -f ../tmp/input_complex.mol")
        return stk.BuildingBlock.init_from_file("../tmp/Optimizer_output_else.mol", functional_groups=[
            stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ], )
    else:
        print("you messed up the option bit of the optimise function")
    os.system("rm -f ../tmp/Optimizer_output_else.mol")





def complex_visualisation(input_complex):
    stk.XyzWriter().write(input_complex, '../tmp/input_complex.xyz')
    with open('../tmp/input_complex.xyz', "r+") as file:
        lines = file.readlines()
        counter = 0
        for i, line in enumerate(lines):
            if len(line.split()) > 0:
                if line.split()[0] == 'Hg':
                    del lines[i]
                    counter += 1
        lines[0] = f"{int(lines[0]) - counter}\n"

        with open('../tmp/input_complex.xyz', "w+") as file:
            file.write(''.join(lines))
        mol_ = io.read('../tmp/input_complex.xyz')
        ase_mol = RCA_Molecule(mol=mol_)
        ase_mol.view_3d()


def ligand_to_mol(ligand: RCA_Ligand, target_path="../tmp/tmp.mol", xyz_path="../tmp/tmp.xyz"):
    xyz_str = ligand.get_xyz_file_format_string()
    with open(xyz_path, "w+") as f:
        f.write(xyz_str)
    os.system(f'obabel .xyz {xyz_path} .mol -O  {target_path} ---errorlevel 1')
    return target_path


def tmp_clean_up(*args):
    for path in args:
        os.system(f"rm -f {path}")


def rotate_tridentate_ligand(tridentate_building_block, x, y, z, index_list) -> float:
    DictionaryA = {}
    DictionaryB = {}
    DictionaryC = {}
    minimum_angleA, minimum_angleB, minimum_angleC = 0, 0, 0
    for degree in np.arange(0, 361, 1.0):
        tridentate_building_block = tridentate_building_block.with_rotation_about_axis(angle=1.0 * (np.pi / 180.0),
                                                                                       axis=np.array((0, 0, 1)),
                                                                                       origin=np.array((0, 0, 0)), )
        A_position = tridentate_building_block.get_atomic_positions(atom_ids=[int(index_list[0] + 1), ])
        B_position = tridentate_building_block.get_atomic_positions(atom_ids=[int(index_list[1] + 1), ])
        C_position = tridentate_building_block.get_atomic_positions(atom_ids=[int(index_list[2] + 1), ])
        dist_A = np.linalg.norm(list(A_position) - np.array((x, y, z)))
        dist_B = np.linalg.norm(list(B_position) - np.array((x, y, z)))
        dist_C = np.linalg.norm(list(C_position) - np.array((x, y, z)))
        DictionaryA[float(degree)] = float(dist_A)
        DictionaryB[float(degree)] = float(dist_B)
        DictionaryC[float(degree)] = float(dist_C)
        minimum_angleA = min(DictionaryA, key=DictionaryA.get)
        minimum_angleB = min(DictionaryB, key=DictionaryB.get)
        minimum_angleC = min(DictionaryC, key=DictionaryC.get)

    tridentate_building_block = tridentate_building_block.with_rotation_about_axis(
        angle=float(minimum_angleA) * (np.pi / 180.0), axis=np.array((0, 0, 1)), origin=np.array((0, 0, 0)), )
    angle_A_position_matrix = tridentate_building_block.get_position_matrix()
    distance_with_origin_A = []
    for i in range(len(angle_A_position_matrix)):
        distance_with_origin_A.append(np.linalg.norm(angle_A_position_matrix[i] - np.array((-10.0, 0, 0))))
    min_distance_with_origin_A = min(distance_with_origin_A)
    tridentate_building_block = tridentate_building_block.with_rotation_about_axis(
        angle=(-1.0) * float(minimum_angleA) * (np.pi / 180.0), axis=np.array((0, 0, 1)), origin=np.array((0, 0, 0)), )

    tridentate_building_block = tridentate_building_block.with_rotation_about_axis(
        angle=float(minimum_angleB) * (np.pi / 180.0), axis=np.array((0, 0, 1)), origin=np.array((0, 0, 0)), )
    angle_B_position_matrix = tridentate_building_block.get_position_matrix()
    distance_with_origin_B = []
    for i in range(len(angle_B_position_matrix)):
        distance_with_origin_B.append(np.linalg.norm(angle_B_position_matrix[i] - np.array((-10.0, 0, 0))))
    min_distance_with_origin_B = min(distance_with_origin_B)
    tridentate_building_block = tridentate_building_block.with_rotation_about_axis(
        angle=(-1.0) * float(minimum_angleB) * (np.pi / 180.0), axis=np.array((0, 0, 1)), origin=np.array((0, 0, 0)), )

    tridentate_building_block = tridentate_building_block.with_rotation_about_axis(
        angle=float(minimum_angleC) * (np.pi / 180.0), axis=np.array((0, 0, 1)), origin=np.array((0, 0, 0)), )
    angle_C_position_matrix = tridentate_building_block.get_position_matrix()
    distance_with_origin_C = []
    for i in range(len(angle_C_position_matrix)):
        distance_with_origin_C.append(np.linalg.norm(angle_C_position_matrix[i] - np.array((-10.0, 0, 0))))
    min_distance_with_origin_C = min(distance_with_origin_C)
    tridentate_building_block = tridentate_building_block.with_rotation_about_axis(
        angle=(-1.0) * float(minimum_angleC) * (np.pi / 180.0), axis=np.array((0, 0, 1)), origin=np.array((0, 0, 0)), )

    if min_distance_with_origin_A > min_distance_with_origin_B and min_distance_with_origin_A > min_distance_with_origin_C:
        return float(minimum_angleA)

    elif min_distance_with_origin_B > min_distance_with_origin_A and min_distance_with_origin_B > min_distance_with_origin_C:
        return float(minimum_angleB)

    elif min_distance_with_origin_C > min_distance_with_origin_B and min_distance_with_origin_C > min_distance_with_origin_A:
        return float(minimum_angleC)

    else:
        print("somethings gone  wrong with the xy plane rotation :(")


def penta_as_tetra(ligand_bb, ligand):
    dict_ = ligand.get_assembly_dict()
    atom_func = dict_["index"]
    atom_type = dict_["type"]

    # translate it so centroid is placed at 0,0,0
    penta_bb_temp = ligand_bb.with_centroid(np.array((0, 0, 0)), atom_ids=atom_func)

    variance_list = []
    for i, _ in enumerate(atom_func):
        list_indices = atom_func.copy()
        del list_indices[i]
        penta_bb_temp = penta_bb_temp.with_rotation_to_minimize_angle(start=penta_bb_temp.get_plane_normal(
            atom_ids=[int(list_indices[0]), int(list_indices[1]), int(list_indices[2]), int(list_indices[3]), ]),
            target=np.array((0, 0, 1)),
            axis=np.array((0, 1, 0)),
            origin=np.array((0, 0, 0)), )

        penta_bb_temp_2 = penta_bb_temp.with_rotation_to_minimize_angle(start=penta_bb_temp.get_plane_normal(
            atom_ids=[int(list_indices[0]), int(list_indices[1]), int(list_indices[2]), int(list_indices[3]), ]),
            target=np.array((0, 0, 1)),
            axis=np.array((1, 0, 0)),
            origin=np.array((0, 0, 0)), )

        position_xyz = list(penta_bb_temp_2.get_atomic_positions(
            atom_ids=[int(list_indices[0]), int(list_indices[1]), int(list_indices[2]), int(list_indices[3]), ]))

        z_list = [position_xyz[0][2], position_xyz[1][2], position_xyz[2][2]]
        variance = np.var(z_list)
        variance_list.append(variance)

    index_ = variance_list.index(min(variance_list))

    #
    # depending on index
    position_index = atom_func[index_]
    atom_ids_ = [afunc_ for i, afunc_ in enumerate(atom_func) if i != index_]

    _mod_functional_groups = [stk.GenericFunctionalGroup(atoms=(getattr(stk_, atype_)(afunc_),),
                                                         bonders=(getattr(stk_, atype_)(afunc_),),
                                                         deleters=()
                                                         ) for k, (atype_, afunc_) in
                              enumerate(zip(atom_type, atom_func))
                              if k != index_]

    tmp_path = "../tmp/lig_mol.mol"

    ligand_to_mol(ligand=ligand, target_path=tmp_path)
    penta_bb_temp = stk.BuildingBlock.init_from_file(tmp_path, functional_groups=_mod_functional_groups)

    os.remove("../tmp/lig_mol.mol")

    return penta_bb_temp, position_index, atom_ids_





def mercury_remover(stk_Building_block):
    print("Removing Temporary Mercury")
    stk.MolWriter().write(stk_Building_block, '../tmp/stk_Building_block.mol')
    os.system('obabel .mol ../tmp/stk_Building_block.mol .xyz -O  ../tmp/stk_Building_block.xyz ---errorlevel 1')
    os.system("rm -f ../tmp/stk_Building_block.mol")
    path = "../tmp/stk_Building_block.xyz"
    with open(path, "r") as f:
        new_str = list()
        counter = 0
        for line in f.readlines():
            if len(line.split()) > 0:
                if line.split()[0] != "Hg":
                    new_str.append(line)
                else:
                    counter += 1
            else:
                new_str.append("\n\n")
        new_str[0] = str(int(new_str[0]) - counter)
    with open(path, "w+") as f:
        f.write(''.join([elem for elem in new_str]))
    os.system('obabel .xyz ../tmp/stk_Building_block.xyz .mol -O  ../tmp/stk_Building_block.mol ---errorlevel 1')
    os.system("rm -f stk_Building_block.xyz")
    stk_Building_block1 = stk.BuildingBlock.init_from_file('../tmp/stk_Building_block.mol')
    os.system("rm -f ../tmp/stk_Building_block.mol")
    return stk_Building_block1


def Bidentate_Rotator(ligand_bb, ligand):
    stk_Building_Block = mercury_remover(ligand_bb)

    index_list = ligand.get_assembly_dict()["index"]

    functional_group_2 = list(stk_Building_Block.get_atomic_positions(atom_ids=index_list[1]))

    vector = list(stk_Building_Block.get_direction(atom_ids=[int(index_list[0]), int(index_list[1]), ]))

    x2, y2, z2 = functional_group_2[0][0], functional_group_2[0][1], functional_group_2[0][2]
    x1, y1, z1 = vector[0], vector[1], vector[2]

    Boxes = get_boxes(denticity=ligand.denticity)

    rotation_increment = 1.0

    dict_ = {value: 0 for value in np.arange(0, 361, rotation_increment)}
    for angle in dict_:
        stk_Building_Block = stk_Building_Block.with_rotation_about_axis(angle=rotation_increment * (np.pi / 180.0),
                                                                         axis=np.array((x1, y1, z1)),
                                                                         origin=np.array((x2, y2, z2)), )
        # movie(stk_Building_Block)
        total_atoms_in_box = 0
        for counter, atom in enumerate(list(stk_Building_Block.get_atomic_positions())):
            point_ = [atom[i] for i in range(3)]
            for Box in Boxes:
                if Box.point_in_box(point=point_):
                    score_x = intensity / (1.0 + (sharpness * ((point_[0]) - ((Box.x2 - Box.x1) / 2.0) + Box.x1) ** 2))
                    score_y = intensity / (1.0 + (sharpness * ((point_[1]) - ((Box.y2 - Box.y1) / 2.0) + Box.y1) ** 2))
                    score_z = intensity / (1.0 + (sharpness * ((point_[2]) - ((Box.z2 - Box.z1) / 2.0) + Box.z1) ** 2))
                    total_atoms_in_box = total_atoms_in_box + score_x + score_y + score_z

        dict_[angle] = float(total_atoms_in_box)

    minimum_angle = min(dict_, key=dict_.get)

    #
    #
    stk_Building_Block = stk_Building_Block.with_rotation_about_axis(angle=minimum_angle * (np.pi / 180.0),
                                                                     axis=np.array((x1, y1, z1)),
                                                                     origin=np.array((x2, y2, z2)), )
    # visualize(stk_Building_Block)
    num_atoms = stk_Building_Block.get_num_atoms()
    stk_Building_Block.write("../tmp/temp_xyz.xyz")
    os.system("echo 'Hg 0.0       0.0        0.0' >> ../tmp/temp_xyz.xyz")
    exec(str(os.system("sed '1 s/.*/" + str(num_atoms + 1) + "/' ../tmp/temp_xyz.xyz >  ../tmp/temp_xyz_2.xyz")))
    os.system('obabel .xyz ../tmp/temp_xyz_2.xyz .mol -O  ../tmp/temp_mol.mol ---errorlevel 1')
    os.system(
        "sed 's/Hg[[:space:]][[:space:]]0[[:space:]][[:space:]]0/Hg  0  2/g' ../tmp/temp_mol.mol > ../tmp/temp_mol_2.mol")

    return "../tmp/temp_mol_2.mol"


def get_energy(molecule):
    path = ligand_to_mol(molecule)  # Here there is a dependency on the above ligand_to_mol function.
    print(path)
    mol = next(pybel.readfile("mol", str(path)))
    obmol = mol.OBMol
    ff = ob.OBForceField_FindType("uff")
    assert (ff.Setup(obmol))
    kj_to_kcal = 1.0 / 4.184
    ff.SetCoordinates(mol.OBMol)
    # false = don't calculate gradients
    uffE = ff.Energy(False) * kj_to_kcal
    return uffE


def nonplanar_tetra_solver(bb, lig):
    all_Energies, all_midpoints = [], []

    combo_list = list(itertools.combinations([0, 1, 2, 3], 2))

    assembly_dictionary = lig.get_assembly_dict()
    location_list = assembly_dictionary["index"]  # location of every coord atom in file
    type_list = assembly_dictionary["type"]

    for combo in combo_list:  # this loop iterates through every combination of functional Groups
        func_1_location = location_list[combo[0]]  # first atom loction
        func_2_location = location_list[combo[1]]  # first atom type

        # we needed to initialise the building block solely to get the positons of the functional groups
        positions = list(
            bb.get_atomic_positions(atom_ids=[int(func_1_location), int(func_2_location), ], ))
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
        lig.add_atom(symbol="Hg", coordinates=[x_mid, y_mid, z_mid])  # add dummy metal
        Energy = get_energy(lig)  # get energy
        lig.remove_last_element_in_xyz()  # get rid of dummy metal
        all_Energies.append(Energy)  # list of all the energies of pacing dummy metal at each midpoint
        all_midpoints.append(mid_points)

    minimum_energy = min(all_Energies)
    minimum_energy_index = all_Energies.index(minimum_energy)
    lig.add_atom(symbol="Hg",
                 coordinates=[all_midpoints[minimum_energy_index][0], all_midpoints[minimum_energy_index][1],
                              all_midpoints[minimum_energy_index][
                                  2]])  # paces Hg at midpoint with smallest energy

    file_path_2 = ligand_to_mol(lig)

    tetra_bb_2 = stk.BuildingBlock.init_from_file(str(file_path_2), functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts="[Hg]", bonders=(0,), deleters=(), ), ], )

    tetra_bb_2 = tetra_bb_2.with_displacement(np.array(((-1) * all_midpoints[minimum_energy_index][0],
                                                        (-1) * all_midpoints[minimum_energy_index][1],
                                                        (-1) * all_midpoints[minimum_energy_index][2])))
    # visualize(tetra_bb_2)

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

    position_of_Hg_in_mol = [key for key, item in lig.coordinates.items() if item[0] == "Hg"][0]

    # The following Block of code ensures that the remaning coordinating groups exist in the xy plane
    complex_tetradentate = tetra_bb_2.with_rotation_to_minimize_angle(start=tetra_bb_2.get_plane_normal(
        atom_ids=[int(location_list[0]), int(location_list[1]),
                  int(position_of_Hg_in_mol), ]), target=np.array((0, 0, 1)), axis=np.array((0, 1, 0)),
        origin=np.array((0, 0, 0)), )

    complex_tetradentate = complex_tetradentate.with_rotation_to_minimize_angle(
        start=complex_tetradentate.get_plane_normal(
            atom_ids=[int(location_list[0]), int(location_list[0]),
                      int(position_of_Hg_in_mol), ]), target=np.array((0, 0, 1)), axis=np.array((1, 0, 0)),
        origin=np.array((0, 0, 0)), )

    for degree in np.arange(0, 361, 0.5):
        complex_tetradentate = complex_tetradentate.with_rotation_about_axis(angle=0.5 * (np.pi / 180.0),
                                                                             axis=np.array((0, 0, 1)),
                                                                             origin=np.array((0, 0, 0)),
                                                                             )

        position_of_coord_atom_1 = complex_tetradentate.get_atomic_positions(atom_ids=[int(location_list[0]), ])
        position_of_coord_atom_2 = complex_tetradentate.get_atomic_positions(atom_ids=[int(location_list[1]), ])

        distance1 = np.linalg.norm(list(position_of_coord_atom_1) - np.array((-10.0, 0, 0)))
        distance2 = np.linalg.norm(list(position_of_coord_atom_2) - np.array((-10.0, 0, 0)))

        mean = (distance1 + distance2) / 2.0
        deviation1 = (distance1 - mean) ** 2
        deviation2 = (distance2 - mean) ** 2
        variance = (deviation1 + deviation2) / 2.0

        if (variance < 0.0001) and (distance1 > 10.0) and (distance2 > 10.0):
            complex_tetradentate.write("../tmp/tmp.mol")
        else:
            pass

    return "../tmp/tmp.mol"

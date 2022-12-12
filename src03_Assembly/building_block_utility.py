import stk
from openbabel import pybel

from openbabel import openbabel as ob
import itertools

import os

from src02_Pre_Ass_Filtering.constants import get_boxes, intensity, sharpness
from src03_Assembly.stk_extension import *
from src01.Molecule import RCA_Ligand


def rotate_tetradentate_bb(tetradentate_bb, ligand_):
    for axis_ in [np.array((0, 1, 0)), np.array((1, 0, 0))]:
        tetradentate_bb = tetradentate_bb.with_rotation_to_minimize_angle(
            start=tetradentate_bb.get_plane_normal(
                atom_ids=[ligand_.get_assembly_dict()["index"][i] for i in range(4)]),
            target=np.array((0, 0, 1)),
            axis=axis_,
            origin=np.array((0, 0, 0))
        )

    return tetradentate_bb


def rotate_tridentate_bb(tridentate_bb_, ligand_):
    index_list = ligand_.get_assembly_dict()["index"]
    for axis_ in [np.array((0, 1, 0)), np.array((1, 0, 0))]:
        tridentate_bb_ = tridentate_bb_.with_rotation_to_minimize_angle(
            start=tridentate_bb_.get_plane_normal(atom_ids=[index_list[0], index_list[1], index_list[2]]),
            target=np.array((0, 0, 1)),
            axis=axis_,
            origin=np.array((0, 0, 0))
        )

    return tridentate_bb_


def rotate_pentadentate_bb(pentadentate_bb,
                           indices,
                           rotations: [np.array, list[np.array]]
                           ):
    if not isinstance(rotations, list):
        rotations = [rotations]

    for rot in rotations:
        pentadentate_bb = pentadentate_bb.with_rotation_to_minimize_angle(
            start=pentadentate_bb.get_plane_normal(atom_ids=[int(ind) for ind in indices]),
            target=np.array((0, 0, 1)),
            axis=rot,
            origin=np.array((0, 0, 0))
        )

    return pentadentate_bb


def get_optimal_rotation_angle_tridentate(tridentate_building_block,
                                          x,
                                          y,
                                          z,
                                          ligand: RCA_Ligand
                                          ) -> float:
    index_list = ligand.get_assembly_dict()["index"]

    # init empty variables
    d = {
        "a": {},
        "b": {},
        "c": {}
    }

    for degree in np.arange(0, 360, 1.0):
        tridentate_building_block = tridentate_building_block.with_rotation_about_axis(angle=1.0 * (np.pi / 180.0),
                                                                                       axis=np.array((0, 0, 1)),
                                                                                       origin=np.array((0, 0, 0)), )
        for i, key in enumerate(d):
            pos = tridentate_building_block.get_atomic_positions(atom_ids=[int(index_list[i] + 1), ])
            dist = np.linalg.norm(list(pos) - np.array((x, y, z)))
            d[key][degree] = float(dist)

    minimum_angle = {key: min(dic, key=dic.get) for key, dic in d.items()}

    distance_with_origin = {key: [] for key in minimum_angle}

    for key in minimum_angle:
        # The tridentate shall be rotated further for each process
        tridentate_building_block = tridentate_building_block.with_rotation_about_axis(
            angle=float(minimum_angle[key]) * (np.pi / 180.0),
            axis=np.array((0, 0, 1)),
            origin=np.array((0, 0, 0))
        )

        angle_A_position_matrix = tridentate_building_block.get_position_matrix()
        for i in range(len(angle_A_position_matrix)):
            distance_with_origin[key].append(np.linalg.norm(angle_A_position_matrix[i] - np.array((-10.0, 0, 0))))

        tridentate_building_block = tridentate_building_block.with_rotation_about_axis(
            angle=(-1.0) * float(minimum_angle[key]) * (np.pi / 180.0), axis=np.array((0, 0, 1)),
            origin=np.array((0, 0, 0)), )

    min_distance_with_origin = {key: min(value) for key, value in distance_with_origin.items()}

    #
    # now we would like to return the angle with the maximal minimum corresponding distance
    max_distance_key = max(min_distance_with_origin, key=min_distance_with_origin.get)
    return minimum_angle[max_distance_key]


def penta_as_tetra(ligand: RCA_Ligand):
    """
    Here we convert a pentadentate ligand to a tetradentate one
    """
    ligand_bb = ligand.to_stk_bb()

    # translate it so centroid is placed at 0,0,0
    penta_building_block = ligand_bb.with_centroid(np.array((0, 0, 0)), atom_ids=ligand.get_assembly_dict()["index"])

    variances = {}
    for i, _ in enumerate(ligand.get_assembly_dict()["index"]):
        list_indices = [ind for k, ind in enumerate(ligand.get_assembly_dict()["index"]) if k != i]

        penta_building_block = rotate_pentadentate_bb(penta_building_block,
                                                      list_indices,
                                                      rotations=[np.array((0, 1, 0)), np.array((1, 0, 0))]
                                                      )

        positions = penta_building_block.get_atomic_positions(atom_ids=[int(ind) for ind in list_indices])

        # now we are interested in the variance of the z-coordinates of the above positions
        variances[i] = np.var([pos[2] for pos in positions])

    index_with_min_variance = max(variances, key=variances.get)

    # now we need the modified functional groups, because we leave out the functional_atom corresponding to the index

    # first we modifiy the atom ids:
    modified_atom_ids = [ind for i, ind in enumerate(ligand.get_assembly_dict()["index"]) if
                         i != index_with_min_variance]

    modified_atom_types = [ind for i, ind in enumerate(ligand.get_assembly_dict()["type"]) if
                           i != index_with_min_variance]

    _mod_functional_groups = [
        stk.GenericFunctionalGroup(atoms=[getattr(stk, a)(i)], bonders=[getattr(stk, a)(i)], deleters=())
        for (a, i) in zip(modified_atom_types, modified_atom_ids)
    ]

    penta_bb_temp = stk.BuildingBlock.init_from_molecule(molecule=ligand.to_stk_mol(),
                                                         functional_groups=_mod_functional_groups
                                                         )

    return penta_bb_temp.with_centroid(np.array((0, 0, 0)), atom_ids=modified_atom_ids), ligand.get_assembly_dict()["index"][index_with_min_variance]


# todo: From here on
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
        #lig.remove_last_element_in_xyz()  # get rid of dummy metal
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

    for _ in np.arange(0, 361, 0.5):
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


def ligand_to_mol(ligand: RCA_Ligand, target_path="../tmp/tmp.mol", xyz_path="../tmp/tmp.xyz"):
    xyz_str = ligand.get_xyz_file_format_string()
    with open(xyz_path, "w+") as f:
        f.write(xyz_str)
    os.system(f'obabel .xyz {xyz_path} .mol -O  {target_path} ---errorlevel 1')
    return target_path


def tmp_clean_up(*args):
    for path in args:
        os.system(f"rm -f {path}")


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

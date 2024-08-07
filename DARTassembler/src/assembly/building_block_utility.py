import os
import stk
import itertools
from rdkit import Chem
from rdkit.Chem import rdmolfiles

from DARTassembler.src.assembly.forcefields import ForceField
from DARTassembler.src.ligand_filters.constants_BoxExcluder import get_boxes, intensity, sharpness
from DARTassembler.src.assembly.stk_extension import *
from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand
from copy import deepcopy
import logging


def rotate_tetradentate_bb(tetradentate_bb, ligand_):
    # This function rotates the tetradentate building block such that the normal axis is parallel with the z-axis.
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
    # This function rotates the tridentate building block such that the normal axis is parallel with the z-axis.
    index_list = ligand_.get_assembly_dict()["index"]
    for axis_ in [np.array((0, 1, 0)), np.array((1, 0, 0))]:
        tridentate_bb_ = tridentate_bb_.with_rotation_to_minimize_angle(
            start=tridentate_bb_.get_plane_normal(atom_ids=[index_list[0], index_list[1], index_list[2]]),
            target=np.array((0, 0, 1)),
            axis=axis_,
            origin=np.array((0, 0, 0))
        )
    return tridentate_bb_


def rotate_pentadentate_bb(pentadentate_bb, indices, rotations: [np.array, list[np.array]]):
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


def Bidentate_Rotator(ligand_bb, ligand, top_list=None, bool_placed=None, build_options: str = None):
    # TODO: rewrite this so it doesn't write to disk
    stk_Building_Block = mercury_remover(ligand_bb)

    index_list = ligand.get_assembly_dict()["index"]

    functional_group_2 = list(stk_Building_Block.get_atomic_positions(atom_ids=index_list[1]))

    vector = list(stk_Building_Block.get_direction(atom_ids=[int(index_list[0]), int(index_list[1]), ]))

    x2, y2, z2 = functional_group_2[0][0], functional_group_2[0][1], functional_group_2[0][2]
    x1, y1, z1 = vector[0], vector[1], vector[2]

    Boxes = get_boxes(denticity=ligand.denticity, input_topology=top_list, bool_placed_boxes=bool_placed, build_options=build_options)

    rotation_increment = 1.0
    dict_box = {}
    dict_ = {value: 0 for value in np.arange(0, 361, rotation_increment)}
    for angle in dict_:
        stk_Building_Block = stk_Building_Block.with_rotation_about_axis(angle=rotation_increment * (np.pi / 180.0),
                                                                         axis=np.array((x1, y1, z1)),
                                                                         origin=np.array((x2, y2, z2)), )
        # movie(stk_Building_Block)
        total_atoms_in_box = 0
        box_entered = []
        for counter, atom in enumerate(list(stk_Building_Block.get_atomic_positions())):
            point_ = [atom[i] for i in range(3)]
            k = 1
            for Box in Boxes:
                if Box.point_in_box(point=point_):
                    score_x = intensity / (1.0 + (sharpness * ((point_[0]) - ((Box.x2 - Box.x1) / 2.0) + Box.x1) ** 2))
                    score_y = intensity / (1.0 + (sharpness * ((point_[1]) - ((Box.y2 - Box.y1) / 2.0) + Box.y1) ** 2))
                    score_z = intensity / (1.0 + (sharpness * ((point_[2]) - ((Box.z2 - Box.z1) / 2.0) + Box.z1) ** 2))
                    total_atoms_in_box = total_atoms_in_box + score_x + score_y + score_z
                    box_entered.append(k)
                k = k + 1
        dict_box.update({str(angle): box_entered})
        dict_[angle] = float(total_atoms_in_box)

    minimum_angle = min(dict_, key=dict_.get)
    logging.debug("minimum angle = " + str(minimum_angle) + " boxes entered " + str(dict_box[str(minimum_angle)]))
    # logging.debug(dict_box)

    #
    #
    stk_Building_Block = stk_Building_Block.with_rotation_about_axis(angle=minimum_angle * (np.pi / 180.0),
                                                                     axis=np.array((x1, y1, z1)),
                                                                     origin=np.array((x2, y2, z2)), )


    """
    centroid = stk_Building_Block.get_centroid()
    # if the angle between the centroid of a bidentate ligand and the xy plane is less than 10 degree then the ligand probably needs to sit planar so we rotate it like that
    if (np.arcsin(list(centroid)[2] / np.linalg.norm(centroid)) * (180 / np.pi)) < 2 and (top_list != [3, 2, 0]):
        logging.debug("forcing bidentate planar")
        if not bool_placed:
            stk_Building_Block = stk_Building_Block.with_rotation_to_minimize_angle(start=centroid,
                                                                                    target=np.array([10, 0, 0]),
                                                                                    axis=np.array([0, 1, 0]),
                                                                                    origin=np.array((x2, y2, z2)))
        elif bool_placed:
            stk_Building_Block = stk_Building_Block.with_rotation_to_minimize_angle(start=centroid,
                                                                                    target=np.array([-10, 0, 0]),
                                                                                    axis=np.array([0, 1, 0]),
                                                                                    origin=np.array((x2, y2, z2)))
    else:
        logging.debug("Not forcing bidentate planar")
        pass
    """

    # visualize(stk_Building_Block)
    # I think after here we need to add a mercury
    tmp_mol = stk_Building_Block.to_rdkit_mol()
    metal_string_output = """
 OpenBabel08012417113D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Hg  0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""
    mol_Hg = rdmolfiles.MolFromMolBlock(metal_string_output, removeHs=False, sanitize=False, strictParsing=False)  # Created rdkit object of just metal atom
    mol_Hg.GetAtomWithIdx(0).SetFormalCharge(2)
    combined_mol = Chem.CombineMols(tmp_mol, mol_Hg)
    stk.BuildingBlock.init_from_rdkit_mol(combined_mol)
    return stk.BuildingBlock.init_from_rdkit_mol(combined_mol)

# Used by the nonplanar tetra solver
# def get_energy(molecule):
#     string = ligand_to_mol(molecule)
#     uffE = UFF().singlepoint(string)
#     return uffE

# The nonplanar tetra solver is not used at the moment and uses some openbabel calls. It is kept here for reference.
# def nonplanar_tetra_solver(stk_bb, ligand):
#     all_Energies, all_midpoints = [], []
#     combo_list = list(itertools.combinations([0, 1, 2, 3], 2))
#     assembly_dictionary = ligand.get_assembly_dict()
#     location_list = assembly_dictionary["index"]  # location of every coord atom in file
#     type_list = assembly_dictionary["type"]
#     for combo in combo_list:  # this loop iterates through every combination of functional Groups
#         ligand_copy = deepcopy(ligand)
#         func_1_location = location_list[combo[0]]  # first atom loction
#         func_2_location = location_list[combo[1]]  # first atom type
#
#         # we needed to initialise the building block solely to get the positons of the functional groups
#         positions = list(stk_bb.get_atomic_positions(atom_ids=[int(func_1_location), int(func_2_location), ], ))
#         x1 = positions[0][0]
#         y1 = positions[0][1]
#         z1 = positions[0][2]
#         x2 = positions[1][0]
#         y2 = positions[1][1]
#         z2 = positions[1][2]
#         x_mid = (x1 + x2) / 2  # Here we get the components of the mid points between coord groups
#         y_mid = (y1 + y2) / 2
#         z_mid = (z1 + z2) / 2
#         mid_points = [x_mid, y_mid, z_mid]
#         # The add_atom functionality does not update the coordinates attribute correctly in the molecule class
#
#         ligand_copy.add_atom(symbol="Hg", coordinates=[x_mid, y_mid, z_mid])  # add dummy metal
#         Energy = get_energy(ligand_copy)  # get energy
#         # lig.remove_last_element_in_xyz()  # get rid of dummy metal
#         all_Energies.append(Energy)  # list of all the energies of pacing dummy metal at each midpoint
#         all_midpoints.append(mid_points)
#         del ligand_copy
#         logging.debug("iteration done")
#     minimum_energy = min(all_Energies)
#     minimum_energy_index = all_Energies.index(minimum_energy)
#     ligand.add_atom(symbol="Hg", coordinates=[all_midpoints[minimum_energy_index][0], all_midpoints[minimum_energy_index][1],
#                                               all_midpoints[minimum_energy_index][2]])  # paces Hg at midpoint with the smallest energy
#     tetra_bb_2 = stk.BuildingBlock.init_from_rdkit_mol(rdmolfiles.MolFromMolBlock(ligand_to_mol(ligand=ligand), removeHs=False, sanitize=False, strictParsing=False), functional_groups=[
#         stk.SmartsFunctionalGroupFactory(smarts="[Hg]", bonders=(0,), deleters=(), ), ], )
#
#     tetra_bb_2 = tetra_bb_2.with_displacement(np.array(((-1) * all_midpoints[minimum_energy_index][0],
#                                                         (-1) * all_midpoints[minimum_energy_index][1],
#                                                         (-1) * all_midpoints[minimum_energy_index][2])))
#
#     # This block of code allows me to remove all the coordinating groups that were previously used to coord the
#     # temporary atom
#     # Note location refers to the location within the file
#     location_list = list(location_list)
#     index1 = combo_list[minimum_energy_index][0]
#     index2 = combo_list[minimum_energy_index][1]
#     value1 = location_list[index1]
#     value2 = location_list[index2]
#     location_list.remove(value1)
#     location_list.remove(value2)
#
#     # this block removes the coord atoms used to coordinate to the temporary metal
#     value1 = type_list[index1]
#     value2 = type_list[index2]
#     type_list.remove(value1)
#     type_list.remove(value2)
#
#     position_of_Hg_in_mol = ligand.atomic_props["atoms"].index("Hg")
#     ligand.del_Hg_atom(coordinates=[all_midpoints[minimum_energy_index][0], all_midpoints[minimum_energy_index][1], all_midpoints[minimum_energy_index][2]])
#
#     # The following Block of code ensures that the remaining coordinating groups exist in the xy plane
#     complex_tetradentate = tetra_bb_2.with_rotation_to_minimize_angle(start=tetra_bb_2.get_plane_normal(
#         atom_ids=[int(location_list[0]), int(location_list[1]),
#                   int(position_of_Hg_in_mol), ]), target=np.array((0, 0, 1)), axis=np.array((0, 1, 0)),
#         origin=np.array((0, 0, 0)), )
#
#     complex_tetradentate = complex_tetradentate.with_rotation_to_minimize_angle(
#         start=complex_tetradentate.get_plane_normal(
#             atom_ids=[int(location_list[0]), int(location_list[0]),
#                       int(position_of_Hg_in_mol), ]), target=np.array((0, 0, 1)), axis=np.array((1, 0, 0)),
#         origin=np.array((0, 0, 0)), )
#
#     for _ in np.arange(0, 361, 0.5):
#         complex_tetradentate = complex_tetradentate.with_rotation_about_axis(angle=0.5 * (np.pi / 180.0),
#                                                                              axis=np.array((0, 0, 1)),
#                                                                              origin=np.array((0, 0, 0)),
#                                                                              )
#
#         position_of_coord_atom_1 = complex_tetradentate.get_atomic_positions(atom_ids=[int(location_list[0]), ])
#         position_of_coord_atom_2 = complex_tetradentate.get_atomic_positions(atom_ids=[int(location_list[1]), ])
#
#         # Here we are minimising the difference in distances between the two functional groups and a point far away on the -x axis
#         distance1 = np.linalg.norm(list(position_of_coord_atom_1) - np.array((-10.0, 0, 0)))
#         distance2 = np.linalg.norm(list(position_of_coord_atom_2) - np.array((-10.0, 0, 0)))
#
#         mean = (distance1 + distance2) / 2.0
#         deviation1 = (distance1 - mean) ** 2
#         deviation2 = (distance2 - mean) ** 2
#         variance = (deviation1 + deviation2) / 2.0
#         if (variance < 0.001) and (distance1 > 10.0) and (distance2 > 10.0):
#             return complex_tetradentate
#         else:
#             pass

# Used by the nonplanar tetra solver
# def ligand_to_mol(ligand: RCA_Ligand):
#     xyz_str = ligand.get_xyz_file_format_string()
#     mol_b = ob.OBMol()
#     conv = ob.OBConversion()
#     conv.SetInAndOutFormats("xyz", "mol")
#     conv.ReadString(mol_b, xyz_str)
#     string = conv.WriteString(mol_b)
#     return string


def tmp_clean_up(*args):
    for path in args:
        os.system(f"rm -f {path}")


def mercury_remover(stk_Building_block):
    mol = stk_Building_block.to_rdkit_mol()
    Hg_index = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 80:
            Hg_index.append(int(atom.GetIdx()))
        else:
            pass
    mol_RW = Chem.EditableMol(mol)
    for index in Hg_index:
        mol_RW.RemoveAtom(index - Hg_index.index(index))
    non_editable_mol = mol_RW.GetMol()
    output_stk_bb = stk.BuildingBlock.init_from_rdkit_mol(non_editable_mol)
    return output_stk_bb


def building_block_to_mol(bb):
    string = stk.MolWriter().to_string(bb)
    return string


def get_energy_stk(building_block):
    if building_block is None:
        logging.debug("in get energy function returning none")
        return None
    else:
        uffE = ForceField().singlepoint(building_block)

        return uffE

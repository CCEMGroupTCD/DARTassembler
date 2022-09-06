# Module-Cian.2

import stk
from stk import *
import numpy as np
from openbabel import pybel
from stk_extension import *
import os
from ase import io
from ASE_Molecule import ASE_Molecule, ASE_Ligand


def build_ligand(type_list, index_list, path_):
    func_dict = {type_: globals()[f"{type_}"] for type_ in type_list}
    atoms_ = [func_dict[type_](index_list[i]) for i, type_ in enumerate(type_list)]

    functional_groups_ = [stk.GenericFunctionalGroup(atoms=(a,), bonders=(a,), deleters=()) for a in atoms_]
    return stk.BuildingBlock.init_from_file(path_, functional_groups=functional_groups_)


def post_process_complex(input_complex, name, visualize_=True, print_to_xyz=True, return_ase=False, path="../data/Assembled_Molecules"):

    ## müssen aus dem complex noch alle Hg Atome entfernen
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

    if print_to_xyz is True:
        with open(f'{path}/{name}.xyz', "w+") as file:
            file.write(''.join(lines))

    if visualize_ is True:
        with open('../tmp/input_complex.xyz', "w+") as file:
            file.write(''.join(lines))
        mol_ = io.read('../tmp/input_complex.xyz')
        ase_mol = ASE_Molecule(mol=mol_)
        ase_mol.view_3d()

        if return_ase is True:
            return ase_mol


def optimize(input_file, option: str):  # option will be stk_to_xyz or stk_to_stk or file_to_file
    if option == "xyz_to_xyz":
        print("optimizing xyz_to_xyz")
        file_extension = input_file[input_file.find(".") + 1:].split()[0]
        mol = next(pybel.readfile(format=str(file_extension), filename=str(input_file)))
        mol.localopt(forcefield='uff', steps=10_000)
        mol.write("xyz", str(input_file) + "_optimization_output_if.xyz")

    elif option == "stk_to_xyz":
        print("optimizing stk_to_xyz")
        stk.MolWriter().write(input_file, 'input_complex.mol')
        mol = next(pybel.readfile("mol", "input_complex.mol"))
        mol.localopt(forcefield='uff', steps=10_000)
        mol.write("xyz", "Optimizer_output_stk_to_xyz.xyz", overwrite="True")
        os.system("rm -f input_complex.mol")

    elif option == "stk_to_mol":
        print("optimizing stk_to_mol")
        stk.MolWriter().write(input_file, 'input_complex.mol')
        mol = next(pybel.readfile("mol", "input_complex.mol"))
        mol.localopt(forcefield='uff', steps=10_000)
        mol.write("mol", "Optimizer_output_stk_to_mol.mol", overwrite="True")
        os.system("rm -f input_complex.mol")

    elif option == "stk_to_stk":
        print("optimizing stk_to_stk")
        stk.MolWriter().write(input_file, 'input_complex.mol')
        mol = next(pybel.readfile("mol", "input_complex.mol"))
        mol.localopt(forcefield='uff', steps=3)
        mol.write("mol", "Optimizer_output_else.mol", overwrite="True")
        os.system("rm -f input_complex.mol")
        return stk.BuildingBlock.init_from_file("Optimizer_output_else.mol", functional_groups=[
            stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ], )
    else:
        print("you messed up the option bit of the optimise function")
    os.system("rm -f Optimizer_output_else.mol")


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


def post_process_tridentate(_metal_bb, _tridentate_bb, index_list, optimize_=False):

    # Here we rotate our tridentate ligand to ensure it sits in the xy plane
    compl_tri = stk.ConstructedMolecule(topology_graph=tridentate(metals=_metal_bb, ligands=_tridentate_bb,),)

    compl_tri = compl_tri.with_rotation_to_minimize_angle(
        start=compl_tri.get_plane_normal(atom_ids=[index_list[0], index_list[1], index_list[2]]),
        target=np.array((0, 0, 1)),
        axis=np.array((0, 1, 0)),
        origin=np.array((0, 0, 0))
    )

    compl_tri = compl_tri.with_rotation_to_minimize_angle(
        start=compl_tri.get_plane_normal(atom_ids=[index_list[0], index_list[1], index_list[2], ]),
        target=np.array((0, 0, 1)),
        axis=np.array((1, 0, 0)),
        origin=np.array((0, 0, 0))
    )

    compl_tri = compl_tri.with_rotation_about_axis(axis=np.array((0, 0, 1)),
                                                   angle=float(np.radians(rotate_tridentate_ligand(compl_tri, 10.0, 0.0,0.0, index_list))),
                                                   origin=np.array((0, 0, 0)))

    complex_tridentate_bb_ = stk.BuildingBlock.init_from_molecule(
        compl_tri,
        functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())]
    )

    if optimize_ is True:
        complex_tridentate_bb_ = optimize(complex_tridentate_bb_, "stk_to_stk")

    return complex_tridentate_bb_


def post_process_bidentate(metal_bb_, bidentate_bb_, optimize_=False):

    complex_bidentate = stk.ConstructedMolecule(topology_graph=bidentate(metals=metal_bb_, ligands=bidentate_bb_))
    complex_bidentate_bb_ = stk.BuildingBlock.init_from_molecule(complex_bidentate, functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])

    if optimize_ is True:
        complex_bidentate_bb_ = optimize(complex_bidentate_bb_, "stk_to_stk")

    return complex_bidentate_bb_


def post_process_monodentate(metal_bb_, monodentate_bb_, optimize_=False):

    complex_monodentate = stk.ConstructedMolecule(topology_graph=monodentate(metals=metal_bb_, ligands=monodentate_bb_))
    complex_monodentate = stk.BuildingBlock.init_from_molecule(complex_monodentate, functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])

    if optimize_ is True:
        complex_monodentate = optimize(complex_monodentate, "stk_to_stk")

    return complex_monodentate


def post_process_tetradentate(metal_bb, tetradentate_bb, ligand_, **kwargs):
    index_list = ligand_.get_assembly_dict()["index"]

    complex_tetradentate = tetradentate_bb.with_rotation_to_minimize_angle(
        start=tetradentate_bb.get_plane_normal(
            atom_ids=[index_list[0], index_list[1], index_list[2], index_list[3]]),
            target=np.array((0, 0, 1)), axis=np.array((0, 1, 0)), origin=np.array((0, 0, 0))
            )
    complex_tetradentate = complex_tetradentate.with_rotation_to_minimize_angle(
        start=complex_tetradentate.get_plane_normal(
            atom_ids=[index_list[0], index_list[1], index_list[2], index_list[3]]),
            target=np.array((0, 0, 1)), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0))
                )

    topo_graph = stk.metal_complex.Porphyrin(metals=metal_bb, ligands=complex_tetradentate)
    complex_tetradentate = stk.ConstructedMolecule(topology_graph=topo_graph)

    complex_tetradentate_bb = stk.BuildingBlock.init_from_molecule(complex_tetradentate, functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), )])

    return complex_tetradentate_bb


def post_process_two_monodentates(metal_bb, ligand_bb_dict, optimize_=False):

    (_, bb_top) = ligand_bb_dict[list(ligand_bb_dict.keys())[0]]
    top_ = stk.ConstructedMolecule(
        topology_graph=monodentate(metals=metal_bb, ligands=bb_top))

    top_ = stk.BuildingBlock.init_from_molecule(top_, functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])

    if optimize_ is True:
        top_ = optimize(top_, "stk_to_stk")

    #
    #
    (_, bb_bottom) = ligand_bb_dict[list(ligand_bb_dict.keys())[1]]
    bottom_ = stk.ConstructedMolecule(
        topology_graph=monodentate_flipped(metals=metal_bb, ligands=bb_bottom))
    bottom_ = stk.BuildingBlock.init_from_molecule(bottom_, functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])

    if optimize_ is True:
        bottom_ = optimize(bottom_, "stk_to_stk")

    return top_, bottom_


def pentadentate_Solver(ligand: ASE_Ligand):

    dict_ = ligand.get_assembly_dict()
    atom_func = dict_["index"]
    atom_type = dict_["type"]

    tmp_path = "../tmp/lig_mol.mol"

    xyz_str = ligand.get_assembly_dict()["str"]
    with open("../tmp/lig_xyz.xyz", "w+") as f:
        f.write(xyz_str)

    os.system('obabel .xyz ../tmp/lig_xyz.xyz .mol -O  ../tmp/lig_mol.mol')

    stk_ = __import__("stk")

    _functional_groups = [stk.GenericFunctionalGroup(atoms=(getattr(stk_, atype_)(afunc_), ),
                                                     bonders=(getattr(stk_, atype_)(afunc_), ),
                                                     deleters=()
                                                    ) for (atype_, afunc_) in zip(atom_type, atom_func)]

    penta_bb_temp_1 = stk.BuildingBlock.init_from_file(tmp_path, functional_groups=_functional_groups)

    # translate it so centroid is placed at 0,0,0
    penta_bb_temp_2 = penta_bb_temp_1.with_centroid(np.array((0, 0, 0)), atom_ids=atom_func)

    variance_list = []
    for i, _ in enumerate(atom_func):
        list_indices = atom_func.copy()
        del list_indices[i]
        penta_bb_temp_2 = penta_bb_temp_2.with_rotation_to_minimize_angle(start=penta_bb_temp_2.get_plane_normal(
            atom_ids=[int(list_indices[0]), int(list_indices[1]), int(list_indices[2]), int(list_indices[3]), ]),
                                                                          target=np.array((0, 0, 1)),
                                                                          axis=np.array((0, 1, 0)),
                                                                          origin=np.array((0, 0, 0)), )

        penta_bb_temp_2 = penta_bb_temp_2.with_rotation_to_minimize_angle(start=penta_bb_temp_2.get_plane_normal(
            atom_ids=[int(list_indices[0]), int(list_indices[1]), int(list_indices[2]), int(list_indices[3]), ]),
                                                                          target=np.array((0, 0, 1)),
                                                                          axis=np.array((1, 0, 0)),
                                                                          origin=np.array((0, 0, 0)), )

        position_xyz = list(penta_bb_temp_2.get_atomic_positions(
            atom_ids=[int(list_indices[0]), int(list_indices[1]), int(list_indices[2]), int(list_indices[3]), ]))
        # print(position_xyz)
        z1 = position_xyz[0][2]
        z2 = position_xyz[1][2]
        z3 = position_xyz[2][2]
        z_list = [z1, z2, z3]
        variance = np.var(z_list)
        variance_list.append(variance)
        i = i + 1

    index_ = variance_list.index(min(variance_list))

    #
    # depending on index
    position_index = atom_func[index_]
    atom_ids_ = [afunc_ for i, afunc_ in enumerate(atom_func) if i != index_]

    _mod_functional_groups = [stk.GenericFunctionalGroup(atoms=(getattr(stk_, atype_)(afunc_), ),
                                                         bonders=(getattr(stk_, atype_)(afunc_), ),
                                                         deleters=()
                                                        ) for k, (atype_, afunc_) in enumerate(zip(atom_type, atom_func))
                                                        if k != index_]

    penta_bb_temp = stk.BuildingBlock.init_from_file(tmp_path, functional_groups=_mod_functional_groups)

    penta_bb_temp = penta_bb_temp.with_centroid(np.array((0, 0, 0)), atom_ids=atom_ids_)

    penta_bb_temp = penta_bb_temp.with_rotation_to_minimize_angle(start=penta_bb_temp.get_plane_normal(atom_ids=atom_ids_),
                                                                  target=np.array((0, 0, 1)),
                                                                  axis=np.array((0, 1, 0)),
                                                                  origin=np.array((0, 0, 0))
                                                                  )

    penta_bb_temp = penta_bb_temp.with_rotation_to_minimize_angle(start=penta_bb_temp.get_plane_normal(atom_ids=atom_ids_),
                                                                  target=np.array((0, 0, 1)),
                                                                  axis=np.array((1, 0, 0)),
                                                                  origin=np.array((0, 0, 0))
                                                                  )
    #
    #
    #
    tip_position = list(penta_bb_temp.get_atomic_positions(atom_ids=[int(position_index), ]))
    if float(tip_position[0][2]) > 0:
        penta_bb_temp = penta_bb_temp.with_rotation_about_axis(angle=np.radians(180),
                                                               axis=np.array((1, 0, 0)),
                                                               origin=np.array((0, 0, 0))
                                                               )

    elif float(tip_position[0][2]) == 0:
        print("SOmething went wrong")
        return 0

    os.remove("../tmp/lig_mol.mol")

    return penta_bb_temp


def post_process_pentadentate(ligand, _metal_bb):

    pentadentate_bb = pentadentate_Solver(ligand)

    complex_pentadentate = stk.ConstructedMolecule(topology_graph=stk.metal_complex.Porphyrin(metals=_metal_bb,
                                                                                              ligands=pentadentate_bb,
                                                                                              ),
                                                   )

    functional_grps = [stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ]

    penta = stk.BuildingBlock.init_from_molecule(complex_pentadentate, functional_groups=functional_grps)

    return penta



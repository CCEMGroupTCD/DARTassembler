from src03_Assembly.utilities_assembly import *
import stk

def convert_raw_monodentate_bb(metal_bb_, monodentate_bb_, optimize_=False, coordinates=None):
    """
    here we convert a raw monodentate building block into the actual format we need
    for the complex assembly
    :metal_bb: raw metal building block
    :monodentate_bb_: raw monodentate building block
    :coordinates: If coordinate shift is needed provide coordinates as an np array
    """

    mono_topology_graph = Monodentate(metals=metal_bb_, ligands=monodentate_bb_)           # build the final Building block

    if coordinates is not None:
        mono_topology_graph.set_ligand_coordinates(coordinates=coordinates)

    intermediate_bb = stk.ConstructedMolecule(topology_graph=mono_topology_graph)
    mono_for_complex = stk.BuildingBlock.init_from_molecule(intermediate_bb, functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])

    if optimize_ is True:
        mono_for_complex = optimize(mono_for_complex, "stk_to_stk")

    return mono_for_complex


def convert_raw_bidendate_bb(metal_bb_, bidentate_bb_, optimize_=False):

    bi_topology_graph = Bidentate(metals=metal_bb_, ligands=bidentate_bb_)

    complex_bidentate = stk.ConstructedMolecule(topology_graph=bi_topology_graph)
    complex_bidentate_bb_ = stk.BuildingBlock.init_from_molecule(complex_bidentate, functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])

    if optimize_ is True:
        complex_bidentate_bb_ = optimize(complex_bidentate_bb_, "stk_to_stk")

    return complex_bidentate_bb_


def convert_raw_planar_tridentate_bb(metal_bb_, tridentate_bb_, index_list, optimize_=False):

    # We first rotate the tridentate building block around two axis

    for axis_ in [np.array((0, 1, 0)), np.array((1, 0, 0))]:
        tridentate_bb_ = tridentate_bb_.with_rotation_to_minimize_angle(
            start=tridentate_bb_.get_plane_normal(atom_ids=[index_list[0], index_list[1], index_list[2]]),
            target=np.array((0, 0, 1)),
            axis=axis_,
            origin=np.array((0, 0, 0))
        )

    tri_topology_graph = Tridentate(metals=metal_bb_, ligands=tridentate_bb_, )

    compl_tri = stk.ConstructedMolecule(topology_graph=tri_topology_graph, )

    compl_tri = compl_tri.with_rotation_about_axis(axis=np.array((0, 0, 1)),
                                                   angle=float(np.radians(rotate_tridentate_ligand(compl_tri, 10.0, 0.0,0.0, index_list))),
                                                   origin=np.array((0, 0, 0)))

    complex_tridentate_bb_ = stk.BuildingBlock.init_from_molecule(compl_tri,
                                                                  functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())]
                                                                  )

    if optimize_ is True:
        complex_tridentate_bb_ = optimize(complex_tridentate_bb_, "stk_to_stk")

    return complex_tridentate_bb_


def convert_raw_planaer_tetradentate_bb(metal_bb, tetradentate_bb, ligand_, optmize=False,  **kwargs):

    # again we need a twofold rotation for the building block
    for axis_ in [np.array((0, 1, 0)), np.array((1, 0, 0))]:

        tetradentate_bb = tetradentate_bb.with_rotation_to_minimize_angle(
                                start=tetradentate_bb.get_plane_normal(atom_ids=[ligand_.get_assembly_dict()["index"][i] for i in range(4)]),
                                target=np.array((0, 0, 1)),
                                axis=axis_,
                                origin=np.array((0, 0, 0))
                )

    tetra_topology_graph = stk.metal_complex.Porphyrin(metals=metal_bb, ligands=tetradentate_bb)

    complex_tetradentate = stk.ConstructedMolecule(topology_graph=tetra_topology_graph)

    complex_tetradentate_bb = stk.BuildingBlock.init_from_molecule(complex_tetradentate, functional_groups=[
        stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), )])

    return complex_tetradentate_bb


def post_process_pentadentate(ligand, _metal_bb):

    pentadentate_bb = pentadentate_Solver(ligand)

    penta_topology_graph = stk.metal_complex.Porphyrin(metals=_metal_bb, ligands=pentadentate_bb, )

    complex_pentadentate = stk.ConstructedMolecule(topology_graph = penta_topology_graph, )

    penta_for_complex = stk.BuildingBlock.init_from_molecule(complex_pentadentate,
                                                 functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ])

    return penta_for_complex



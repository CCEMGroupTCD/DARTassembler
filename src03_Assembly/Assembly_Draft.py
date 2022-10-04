from src03_Assembly.Convert_Buildingblocks import *


def four_one_one_assembly(metal_bb, final_metal_bb, ligand_bb_dict, optimize_, top_done=False, planar_=True):
    """
    Assembly of a TMC with topology 4-1-1;
    either for a planar or a non-planar tetradentate ligand
    """
    mono_one_bb_for_comp, mono_two_bb_for_comp, tetra_bb_for_comp = None, None, None

    if planar_ is True:
        for key, (lig, lig_bb) in ligand_bb_dict.items():
            if lig.denticity == 4:
                tetra_bb_for_comp = convert_raw_planaer_tetradentate_bb(metal_bb=metal_bb,
                                                                        tetradentate_bb=lig_bb,
                                                                        ligand_=lig
                                                                        )
            elif lig.denticity == 1 and top_done is True:
                mono_one_bb_for_comp = convert_raw_monodentate_bb(metal_bb_=metal_bb,
                                                                  monodentate_bb_=lig_bb,
                                                                  optimize_=optimize_,
                                                                  coordinates=np.array([0.0, 0.0, -1.9])
                                                                  )
            elif lig.denticity == 1 and top_done is False:
                top_done = True
                mono_two_bb_for_comp = convert_raw_monodentate_bb(metal_bb_=metal_bb,
                                                                  monodentate_bb_=lig_bb,
                                                                  optimize_=optimize_
                                                                  )

        complex_top = complex_topology_three(metals=final_metal_bb,
                                             ligands={tetra_bb_for_comp: (0,),
                                                      mono_one_bb_for_comp: (1,),
                                                      mono_two_bb_for_comp: (2,)}
                                             )

        complex_ = stk.ConstructedMolecule(topology_graph=complex_top)

        return complex_

    elif planar_ is False:
        print("Non-Planar Tetradentate Ligand")
        return None


def three_two_one_assembly(metal_bb, final_metal_bb, ligand_bb_dict, optimize_, planar_):
    """
    Assembly for a 3-2-1 topoology TMC
    either for a planar or non-planar tridentate ligand
    """
    mono_bb_for_comp, bi_bb_for_comp, tri_bb_for_comp = None, None, None

    if planar_ is True:

        for (lig, lig_bb) in ligand_bb_dict.values():
            if lig.denticity == 3:
                tri_bb_for_comp = convert_raw_planar_tridentate_bb(metal_bb_=metal_bb,
                                                                   tridentate_bb_=lig_bb,
                                                                   index_list=lig.get_assembly_dict()["index"],
                                                                   optimize_=optimize_
                                                                   )
            elif lig.denticity == 2:
                bi_bb_for_comp = convert_raw_bidendate_bb(metal_bb_=metal_bb,
                                                          bidentate_bb_=lig_bb,
                                                          optimize_=optimize_,
                                                          ligand_=lig
                                                          )
            elif lig.denticity == 1:
                mono_bb_for_comp = convert_raw_monodentate_bb(metal_bb_=metal_bb,
                                                              monodentate_bb_=lig_bb,
                                                              optimize_=optimize_)

        complex_top = complex_topology_three(metals=final_metal_bb,
                                             ligands={bi_bb_for_comp: (0,),
                                                      mono_bb_for_comp: (1,),
                                                      tri_bb_for_comp: (2,), }
                                             )

        complex_ = stk.ConstructedMolecule(topology_graph=complex_top)

        return complex_

    else:
        print("Non-Planar Tridentate Ligand")
        return None


def five_one_assembly(metal_bb, final_metal_bb, ligand_bb_dict, optimize_, **kwargs):
    mono_bb_for_comp, penta_bb_for_comp = None, None

    for (lig, lig_bb) in ligand_bb_dict.values():
        if lig.denticity == 5:
            penta_bb_for_comp = convert_raw_pentadentate_bb(penta_ligand=lig,
                                                            penta_bb=lig_bb,
                                                            metal_bb=metal_bb,
                                                            optimize=optimize_
                                                            )
            if penta_bb_for_comp is None:
                return None
        elif lig.denticity == 1:
            mono_bb_for_comp = convert_raw_monodentate_bb(metal_bb_=metal_bb,
                                                          monodentate_bb_=lig_bb,
                                                          optimize_=optimize_
                                                          )

    _complex_top_ = complex_topology_two(metals=final_metal_bb,
                                         ligands={penta_bb_for_comp: (0,), mono_bb_for_comp: (1,)}
                                         )

    complex_ = stk.ConstructedMolecule(topology_graph=_complex_top_)

    return complex_








import rdkit
from src03_Assembly.Convert_Buildingblocks import *
import random


def post_process_complex(input_complex, name, visualize_=True, print_to_xyz=True, return_ase=False,
                         path="../data/Assembled_Molecules"):
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
        ase_mol = RCA_Molecule(mol=mol_)
        ase_mol.view_3d()

        if return_ase is True:
            return ase_mol


def four_one_one_assembly(metal_bb, final_metal_bb, ligand_bb_dict, _optimize, top_done=False):
    mono_one_bb_for_comp, mono_two_bb_for_comp, tetra_bb_for_comp = None, None, None

    planar_ = planar_ceck(ligand_bb_dict)

    if planar_ is True:
        for key, (lig, lig_bb) in ligand_bb_dict.items():
            if lig.denticity == 4:
                tetra_bb_for_comp = convert_raw_planaer_tetradentate_bb(
                    metal_bb=metal_bb,
                    tetradentate_bb=lig_bb,
                    ligand_=lig
                )
            elif lig.denticity == 1 and top_done is True:
                mono_one_bb_for_comp = convert_raw_monodentate_bb(metal_bb,
                                                                  lig_bb,
                                                                  optimize_=_optimize,
                                                                  coordinates=np.array([0.0, 0.0, -1.9])
                                                                  )
            elif lig.denticity == 1 and top_done is False:
                top_done = True
                mono_two_bb_for_comp = convert_raw_monodentate_bb(metal_bb, lig_bb, optimize_=_optimize)

        complex_top = complex_topology_three(metals=final_metal_bb,
                                             ligands={tetra_bb_for_comp: (0,),
                                                      mono_one_bb_for_comp: (1,),
                                                      mono_two_bb_for_comp: (2,)}
                                             )

        complex_ = stk.ConstructedMolecule(topology_graph=complex_top)

        return complex_

    elif planar_ is False:
        pass


def three_two_one_assembly(metal_bb, final_metal_bb, ligand_bb_dict, _optimize):
    mono_bb_for_comp, bi_bb_for_comp, tri_bb_for_comp = None, None, None

    for (lig, lig_bb) in ligand_bb_dict.values():
        if lig.denticity == 3:
            tri_bb_for_comp = convert_raw_planar_tridentate_bb(metal_bb_=metal_bb,
                                                               tridentate_bb_=lig_bb,
                                                               index_list=lig.get_assembly_dict()["index"],
                                                               optimize_=_optimize)
        elif lig.denticity == 2:
            bi_bb_for_comp = convert_raw_bidendate_bb(metal_bb, lig_bb, optimize_=_optimize)
        elif lig.denticity == 1:
            mono_bb_for_comp = convert_raw_monodentate_bb(metal_bb, lig_bb, optimize_=_optimize)

    complex_top = complex_topology_three(metals=final_metal_bb,
                                         ligands={bi_bb_for_comp: (0,),
                                                  mono_bb_for_comp: (1,),
                                                  tri_bb_for_comp: (2,), }
                                         )

    complex_ = stk.ConstructedMolecule(topology_graph=complex_top)

    return complex_


def five_one_assembly(metal_bb, final_metal_bb, ligand_bb_dict, _optimize):
    mono_bb_for_comp, penta_bb_for_comp = None, None

    for (lig, lig_bb) in ligand_bb_dict.values():
        if lig.denticity == 5:
            penta_bb_for_comp = post_process_pentadentate(ligand=lig, _metal_bb=metal_bb)
        elif lig.denticity == 1:
            mono_bb_for_comp = convert_raw_monodentate_bb(metal_bb, lig_bb, optimize_=_optimize)

    _complex_top_ = complex_topology_two(metals=final_metal_bb,
                                         ligands={penta_bb_for_comp: (0,), mono_bb_for_comp: (1,)}
                                         )

    complex_ = stk.ConstructedMolecule(topology_graph=_complex_top_)

    return complex_


# todo das in die assmebler Klasse
# mit global oder local dict von diesem module hier
def assembly(metal_bb, final_metal_bb, ligand_bb_dict, comp, _optimize=True):
    if set(comp) == {3, 2, 1}:
        complex_ = three_two_one_assembly(metal_bb, final_metal_bb, ligand_bb_dict, _optimize)

    elif set(comp) == {4, 1, 1}:
        # we assume always planar
        complex_ = four_one_one_assembly(metal_bb, final_metal_bb, ligand_bb_dict, _optimize)

    elif set(comp) == {5, 1}:
        complex_ = five_one_assembly(metal_bb, final_metal_bb, ligand_bb_dict, _optimize)
    else:
        print("No valid composition choosen!")
        complex_ = None

    return complex_


def random_assembly(ligand_dict: dict, comps: list, safe_path: str, metals, visualize_, _optimize=False):
    try:
        # choose random metal center
        (metal, charge) = random.choice(metals)

        # build the center atom with 6 connections
        metal_bb = stk.BuildingBlock(smiles='[Hg+2]',
                                     functional_groups=(stk.SingleAtom(stk.Hg(0, charge=2)) for i in range(6)),
                                     position_matrix=np.ndarray([0, 0, 0])
                                     )

        # build the metal block with the new metal atom
        smiles_str = f"[{metal}{charge}]"
        stk_metal_func = globals()[f"{metal}"]
        # stk_metal_func = getattr(__import__("stk"), metal)
        functional_groups = (stk.SingleAtom(stk_metal_func(0, charge=charge)) for i in range(6))
        final_metal_bb = stk.BuildingBlock(smiles=smiles_str,
                                           functional_groups=functional_groups,
                                           position_matrix=np.ndarray([0, 0, 0])
                                           )
        #
        #
        # chose a random composition of the possible comps
        comp = random.choice(comps)

        # and choose a random set of ligands
        ligands = {i: random.choice(ligand_dict[index]) for i, index in enumerate(comp)}

        # dict: ligand: ligand_bb
        # nr : (ligand, ligand_bb)
        ligand_bb_dict = {}

        for k, lig in enumerate(ligands.values()):
            lig_assembly_dict = lig.get_assembly_dict()

            xyz_str = lig_assembly_dict["str"]
            with open("../tmp/lig_xyz.xyz", "w+") as f:
                f.write(xyz_str)

            os.system('obabel .xyz ../tmp/lig_xyz.xyz .mol -O  ../tmp/lig_mol.mol')

            ligand_bb_dict[k] = (lig, build_ligand(type_list=lig_assembly_dict["type"],
                                                   index_list=lig_assembly_dict["index"],
                                                   path_="../tmp/lig_mol.mol"))
            os.remove("../tmp/lig_mol.mol")

        # Ligand_bb_dict contains the RAW building blocks, which are getting converted into actual buildings block later down the road

        complex_ = assembly(metal_bb, final_metal_bb, ligand_bb_dict, comp, _optimize=_optimize)

        post_process_complex(complex_, name={f"{ligands[0].csd_code}"},
                             visualize_=visualize_, print_to_xyz=True, path=safe_path)

        return complex_

    except rdkit.Chem.rdchem.AtomValenceException as ex:
        print(f"The standard Error: {ex} has occured. Solveable by turning sanitize in stk off")
        return -1

    except Exception as e:
        print(f"Oh. A new error: {e}!!")
        return -1

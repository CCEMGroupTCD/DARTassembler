
from src03_Assembly.building_block_utility import *
from src02_Pre_Ass_Filtering.constants import get_boxes
from ase import io
from stk import *
import logging


def delete_Hg(complex_):
    stk.MolWriter().write(complex_, '../tmp/complex.mol')
    os.system('obabel .mol ../tmp/complex.mol .xyz -O  ../tmp/complex.xyz ---errorlevel 1')
    path = "../tmp/complex.xyz"
    os.system("rm -f ../tmp/complex.mol")
    with open(path, "r") as f:
        new_str = []
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

    os.system('obabel .xyz ../tmp/complex.xyz .mol -O  ../tmp/complex.mol ---errorlevel 1')

    return stk.BuildingBlock.init_from_file('../tmp/complex.mol')


def visualize_complex():
    mol_ = io.read('../tmp/complex.xyz')
    ase_mol = RCA_Molecule(mol=mol_)
    ase_mol.view_3d()
    os.system("rm -f ../tmp/complex.xyz")


def box_excluder_descision(stk_Building_Block, denticity_, planar=True, threshhold=0.):

    boxes = get_boxes(denticity=denticity_, planar=planar)
    atom_positions = list(stk_Building_Block.get_atomic_positions())

    score = 0

    for atom in atom_positions:  # loop through atoms
        point = [atom[j] for j in range(3)]

        for Box in boxes:
            if Box.point_in_box(point):
                score += Box.intensity / (1.0 + (
                            Box.sharpness * ((point[0]) - ((Box.x2 - Box.x1) / 2.0) + Box.x1) ** 2))
                score += Box.intensity / (1.0 + (
                            Box.sharpness * ((point[1]) - ((Box.y2 - Box.y1) / 2.0) + Box.y1) ** 2))
                score += Box.intensity / (1.0 + (
                            Box.sharpness * ((point[2]) - ((Box.z2 - Box.z1) / 2.0) + Box.z1) ** 2))
                break

    if score > threshhold:
        # did not pass
        return False
    else:
        return True


def box_filter(ligand: RCA_Ligand, optimize_=True, box_default_descicion=False) -> bool:
    """
    Returns True if a ligand looks good and can pass
    and false if not
    """
    try:
        metal, charge = "Fe", "+2"

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

        lig_assembly_dict = ligand.get_assembly_dict()

        xyz_str = lig_assembly_dict["str"]
        with open("../tmp/lig_xyz.xyz", "w+") as f:
            f.write(xyz_str)

        os.system('obabel .xyz ../tmp/lig_xyz.xyz .mol -O  ../tmp/lig_mol.mol ---errorlevel 1')

        ligand_bb = build_ligand(type_list=lig_assembly_dict["type"],
                                 index_list=lig_assembly_dict["index"],
                                 path_="../tmp/lig_mol.mol"
                                 )
        os.remove("../tmp/lig_mol.mol")

        # build the building blocks
        if ligand.denticity == 3 and ligand.planar_check() is True:
            bb_for_comp = convert_raw_planar_tridentate_bb(metal_bb_=metal_bb,
                                                           tridentate_bb_=ligand_bb,
                                                           index_list=ligand.get_assembly_dict()["index"],
                                                           optimize_=True
                                                           )

        elif ligand.denticity == 3 and ligand.planar_check() is False:
            logging.info("Not implemented yet")
            print("Not implemented yet")
            return True                     # will always let it through
        elif ligand.denticity == 4 and ligand.planar_check() is True:
            bb_for_comp = convert_raw_planar_tetradentate_bb(metal_bb=metal_bb,
                                                             tetradentate_bb=ligand_bb,
                                                             ligand_=ligand,
                                                             optmize=optimize_
                                                             )
        elif ligand.denticity == 4 and ligand.planar_check() is False:
            print("Not implemented yet")
            return True                     # will always let it through
        elif ligand.denticity == 5:
            bb_for_comp = convert_raw_pentadentate_bb(metal_bb=metal_bb,
                                                      penta_bb=ligand_bb,
                                                      penta_ligand=ligand,
                                                      optimize=optimize_
                                                      )
        elif ligand.denticity == 2:
            bb_for_comp = convert_raw_bidendate_bb(metal_bb_=metal_bb,
                                                   bidentate_bb_=ligand_bb,
                                                   ligand_=ligand,
                                                   optimize_=optimize_,
                                                   )
        else:
            print("Something went wrong")
            raise ValueError

        #
        #
        #
        # convert the building blocks to topologies
        if ligand.denticity == 5:
            complex_ = stk.ConstructedMolecule(topology_graph=complex_topology_two(metals=final_metal_bb,
                                                                                   ligands=bb_for_comp
                                                                                   )
                                               )

            complex_ = delete_Hg(complex_)

            return box_excluder_descision(complex_, denticity_=ligand.denticity)
        else:
            complex_ = stk.ConstructedMolecule(topology_graph=complex_topology_three(metals=final_metal_bb,
                                                                                     ligands=bb_for_comp
                                                                                     )
                                               )

            complex_ = delete_Hg(complex_)

            return box_excluder_descision(complex_, denticity_=ligand.denticity, planar=ligand.planar_check())

    except Exception as e:
        print(f"Box Excluder filter failed, will return default ({box_default_descicion}) (I.e. dont let it pass)"
              f"\n Reason for failing: {e}")
        return box_default_descicion



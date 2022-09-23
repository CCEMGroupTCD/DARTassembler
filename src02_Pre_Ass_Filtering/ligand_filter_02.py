# Box Excluder coming in here

import pickle
from src03_Assembly.utilities_assembly import *
from src02_Pre_Ass_Filtering.constants import get_boxes
from ase import io


def delete_Hg(complex_):
    stk.MolWriter().write(complex_, '../tmp/complex.xyz')
    os.system('obabel ../tmp/complex.xyz .xyz -O  ../tmp/complex.xyz')
    os.system("rm -f complex.mol")

    path = "../tmp/complex.xyz"
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

    os.system('obabel .xyz ../tmp/complex.xyz .mol -O  ../tmp/complex.mol')

    return stk.BuildingBlock.init_from_file('../tmp/complex.mol')


def visualize_complex():
    mol_ = io.read('../tmp/complex.xyz')
    ase_mol = RCA_Molecule(mol=mol_)
    ase_mol.view_3d()
    os.system("rm -f ../tmp/complex.xyz")


def box_excluder(stk_Building_Block, denticity, threshhold=0.):

    boxes = get_boxes(denticity=denticity)
    atom_positions = list(stk_Building_Block.get_atomic_positions())

    score = 0

    for atom in atom_positions:  # loop through atoms
        point = [atom[j] for j in range(3)]

        for Box in boxes:
            if Box.point_in_box(point):
                score += (Box.intensity) / (1.0 + ((Box.sharpness) * ((point[0]) - ((Box.x2 - Box.x1) / 2.0) + Box.x1) ** 2))
                score += (Box.intensity) / (1.0 + ((Box.sharpness) * ((point[1]) - ((Box.y2 - Box.y1) / 2.0) + Box.y1) ** 2))
                score += (Box.intensity) / (1.0 + ((Box.sharpness) * ((point[2]) - ((Box.z2 - Box.z1) / 2.0) + Box.z1) ** 2))
                break

    if score > threshhold:
        # did not pass
        return False
    else:
        return True


if __name__ == "__main__":

    with open("../data/ligand_dict.pickle", "rb") as handle:  # this command initialises the Dictionary of all ligands
        ligand_dict = pickle.load(handle)

    #
    #
    #
    # Test-setup
    # todo: hier spaeter einen for loop einbauen
    denticity = 3
    ligand_list = ligand_dict[denticity][:10]

    evaluation = {"pt": 0, "pf": 0, "ft": 0, "ff": 0}
    #
    #
    tmp_clean_up("../tmp/tmp.xyz", "../tmp/tmp.mol")

    for ligand in ligand_list:
        (metal, charge) = ("Fe", "+2")

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

        os.system('obabel .xyz ../tmp/lig_xyz.xyz .mol -O  ../tmp/lig_mol.mol')

        ligand_bb = build_ligand(type_list=lig_assembly_dict["type"],
                                 index_list=lig_assembly_dict["index"],
                                 path_="../tmp/lig_mol.mol"
                                 )
        os.remove("../tmp/lig_mol.mol")

        tri_bb_for_comp = post_process_tridentate(_metal_bb=metal_bb, _tridentate_bb=ligand_bb,
                                                  index_list=ligand.get_assembly_dict()["index"],
                                                  optimize_=True)

        complex = stk.ConstructedMolecule(topology_graph=complex_topology(metals=final_metal_bb,
                                                                          ligands=tri_bb_for_comp
                                                                          )
                                          )

        complex_ = delete_Hg(complex)

        decision = box_excluder(complex_, denticity=denticity)

        visualize_complex()

        var = input(f"Descicion of BoxExcluder was: {decision} (False=no pass)\n What is your decision? (p for pass, n for no pass)")

        if var == "p" and decision is True:
            evaluation["pt"] += 1
        elif var == "n" and decision is True:
            evaluation["pf"] += 1
        elif var == "p" and decision is False:
            evaluation["ft"] += 1
        else:
            evaluation["ff"] += 1


    for key, item in evaluation.items():
        print(f"number of {key} is {item}\n")
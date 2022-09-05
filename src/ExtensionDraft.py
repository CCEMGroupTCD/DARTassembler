# todo: später durch richtiges modul ersetzen
from utilities_assembly import *


def visualize(input_complex):
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
        ase_mol = ASE_Molecule(mol=mol_)
        ase_mol.view_3d()


def get_ligand_bb(ligands):
    # dict: ligand: ligand_bb
    ligand_bb_dict = {}

    for lig in ligands.values():
        lig_assembly_dict = lig.get_assembly_dict()

        xyz_str = lig_assembly_dict["str"]
        with open("../tmp/lig_xyz.xyz", "w+") as f:
            f.write(xyz_str)

        # todo: convert .xyz to .mol
        os.system('obabel .xyz tmp/lig_xyz.xyz .mol -O  tmp/lig_mol.mol')

        ligand_bb_dict[lig] = build_ligand(type_list=lig_assembly_dict["type"],
                                           index_list=lig_assembly_dict["index"],
                                           path_="../tmp/lig_mol.mol")

    return ligand_bb_dict


def assembly(metal_bb, final_metal_bb, ligands, comp, _optimize=True):
    # todo: das müsste man noch umbauen, weil hier spiel glaube ich tatsächlich dann die Topologie mit rein
    # todo: später in das Stuff_Cian_new modul verschieben oder das hier und das Stuff_Cian_new modul mergen
    if set(comp) == {1, 2, 3}:

        ligand_bb_dict = get_ligand_bb(ligands)

        mono_bb_for_comp, bi_bb_for_comp, tri_bb_for_comp = None, None, None

        for lig, lig_bb in ligand_bb_dict.items():
            if lig.denticity == 3:
                tri_bb_for_comp = post_process_tridentate(_metal_bb=metal_bb, _tridentate_bb=lig_bb,
                                                        index_list=lig.get_assembly_dict()["index"],
                                                        optimize_=_optimize)
            elif lig.denticity == 2:
                bi_bb_for_comp = post_process_bidentate(metal_bb, lig_bb, optimize_=_optimize)
            elif lig.denticity == 1:
                mono_bb_for_comp = post_process_monodentate(metal_bb, lig_bb, optimize_=_optimize)

        complex_top = complex_topology(metals=final_metal_bb,
                                       ligands={bi_bb_for_comp: (0,),
                                                mono_bb_for_comp: (1,),
                                                tri_bb_for_comp: (2,), }
                                       )

        complex_ = stk.ConstructedMolecule(topology_graph=complex_top)

    elif set(comp) == {1, 1, 2, 2}:
        # kann zu 1,2,3-Fall vereinfacht werden?!
        if ligands[2].denticity != 2 or ligands[3].denticity != 1:
            print("Wrong format: composition needs to be (2,1,2,1)!!!")
            return 0

        new_ligand = merge(ligands[3], ligands[2])

        new_ligands = {0: ligands[0], 1: ligands[1], 2: new_ligand}

        return assembly(metal_bb, final_metal_bb, new_ligands, (3, 2, 1), _optimize)

    else:
        complex_ = None

    return complex_


def random_assembly(num, ligand_dict: dict, comps: list, safe_path: str, metals, visualize_):

    generated_complexes = list()

    # for all ligands we want to create
    for i in range(num):
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

            complex_ = assembly(metal_bb, final_metal_bb, ligands, comp)

            post_process_complex(complex_, name={f"{ligands[0].csd_code}"},
                                 visualize_=visualize_, print_to_xyz=True, path=safe_path)

            generated_complexes.append(complex_)

        except Exception as ex:
            print(f"An expected error has occured: {ex}")
            continue

    return generated_complexes


if __name__ == '__main__':
    # Input setting for generating ligands
    #
    # number of ligands to generate
    number_of_ligands = 1

    #
    # the ligand dict
    with open("../data/ligand_dict.pickle", "rb") as handle:
        ligand_dict = pickle.load(handle)

    #
    # all denticity combinations
    possible_compositions = [(2,1,2,1)
                             #(3, 2, 1)
                             # (3, 3)
                             # ...
                             ]

    #
    # where to store the resulting ligands
    store_path = "../data/Assembled_Molecules"

    #
    # list of possible metal centers with respective charge
    list_of_metals = [("Fe", "+2")
                      #("Fe", 3)
                      # ...
                      ]

    #
    # visualize the constructed molecules during the process
    visualize_ = True

    complexes = random_assembly(num=number_of_ligands,
                    ligand_dict=ligand_dict,
                    comps=possible_compositions,
                    safe_path=store_path,
                    metals=list_of_metals,
                    visualize_=visualize_
                    )


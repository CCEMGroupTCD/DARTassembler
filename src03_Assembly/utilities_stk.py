import numpy as np
from ase import io

import stk
stk_ = __import__("stk")

from src01.Molecule import RCA_Molecule


def create_placeholder_Hg_bb() -> stk.BuildingBlock:
    """
    We use that frequently throughout so we shall return that special building block
    center atom with 6 connections
    """
    return stk.BuildingBlock(smiles='[Hg+2]',
                             functional_groups=(stk.SingleAtom(stk.Hg(0, charge=2)) for i in range(6)),
                             position_matrix=np.ndarray([0, 0, 0])
                             )


def build_ligand(type_list, index_list, path_):
    func_dict = {type_: getattr(stk_, type_) for type_ in type_list}
    atoms_ = [func_dict[type_](index_list[i]) for i, type_ in enumerate(type_list)]

    functional_groups_ = [stk.GenericFunctionalGroup(atoms=(a,), bonders=(a,), deleters=()) for a in atoms_]
    return stk.BuildingBlock.init_from_file(path_, functional_groups=functional_groups_)


def remove_Hg(input_complex, visualize_: bool = True):
    stk.XyzWriter().write(input_complex, '../tmp/input_complex.xyz')
    with open('../tmp/input_complex.xyz', "r+") as file:
        lines = file.readlines()
        counter = 0
        new_lines = ["0", ""]
        for i, line in enumerate(lines):
            if len(line.split()) > 0:
                if line.split()[0] == 'Hg':
                    counter += 1
                else:
                    new_lines.append(line)
        new_lines[0] = f"{int(lines[0]) - counter}\n"

    with open('../tmp/input_complex.xyz', "w+") as file:
        file.write(''.join(new_lines))
    mol_ = io.read('../tmp/input_complex.xyz')
    rca_mol = RCA_Molecule(mol=mol_)
    if visualize_ is True:
        rca_mol.view_3d()

    return rca_mol, ''.join(new_lines)

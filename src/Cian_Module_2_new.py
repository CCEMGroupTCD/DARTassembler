## Module-Cian.2

import stk
from stk import *
import numpy as np
from openbabel import pybel
from Stk_Extension import *
import os
from ase import io
from ASE_Molecule import ASE_Molecule
import pickle
import random


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


def build_ligand(type_list, index_list, path_):
    func_dict = {type_: globals()[f"{type_}"] for type_ in type_list}
    atoms_ = [func_dict[type_](index_list[i]) for i, type_ in enumerate(type_list)]

    functional_groups_ = [stk.GenericFunctionalGroup(atoms=(a,), bonders=(a,), deleters=()) for a in atoms_]
    return stk.BuildingBlock.init_from_file(path_, functional_groups=functional_groups_)



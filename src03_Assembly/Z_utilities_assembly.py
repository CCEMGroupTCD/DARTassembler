import stk
from openbabel import pybel
from src03_Assembly.stk_extension import *
import os
from ase import io
from src01.Molecule import RCA_Molecule, RCA_Ligand
import rdkit
from src02_Pre_Ass_Filtering.constants import get_boxes, intensity, sharpness
from openbabel import openbabel as ob
import itertools
stk_ = __import__("stk")


def get_topology_string(top):
    """
    Issue is that we need to convert our topology_list into a string
    [1,2,3] -> "one_two_three"
    """
    names_dict = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}
    str_ = ""
    for el in top:
        str_ += f"{names_dict[el]}_"

    return f"{str_}assembly"


def complex_visualisation(input_complex):
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
        ase_mol = RCA_Molecule(mol=mol_)
        ase_mol.view_3d()




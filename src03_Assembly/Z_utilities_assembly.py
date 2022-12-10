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


def ligand_to_mol(ligand: RCA_Ligand, target_path="../tmp/tmp.mol", xyz_path="../tmp/tmp.xyz"):
    xyz_str = ligand.get_xyz_file_format_string()
    with open(xyz_path, "w+") as f:
        f.write(xyz_str)
    os.system(f'obabel .xyz {xyz_path} .mol -O  {target_path} ---errorlevel 1')
    return target_path


def tmp_clean_up(*args):
    for path in args:
        os.system(f"rm -f {path}")







def mercury_remover(stk_Building_block):
    print("Removing Temporary Mercury")
    stk.MolWriter().write(stk_Building_block, '../tmp/stk_Building_block.mol')
    os.system('obabel .mol ../tmp/stk_Building_block.mol .xyz -O  ../tmp/stk_Building_block.xyz ---errorlevel 1')
    os.system("rm -f ../tmp/stk_Building_block.mol")
    path = "../tmp/stk_Building_block.xyz"
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
    os.system('obabel .xyz ../tmp/stk_Building_block.xyz .mol -O  ../tmp/stk_Building_block.mol ---errorlevel 1')
    os.system("rm -f stk_Building_block.xyz")
    stk_Building_block1 = stk.BuildingBlock.init_from_file('../tmp/stk_Building_block.mol')
    os.system("rm -f ../tmp/stk_Building_block.mol")
    return stk_Building_block1



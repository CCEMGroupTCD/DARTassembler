import stk
import os
import collections
from ase import io
from RCA_Molecule import *


# some global coordinates
metals_in_pse = [el for a in [[21, 31], [39, 49], [57, 81], [89, 113]] for el in range(a[0], a[1])]


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
    xyz_str = ligand.xyz.get_xyz_file_format_string()
    with open(xyz_path, "w+") as f:
        f.write(xyz_str)
    os.system(f'obabel .xyz {xyz_path} .mol -O  {target_path}')


def tmp_clean_up(*args):
    for path in args:
        os.system(f"rm -f {path}")



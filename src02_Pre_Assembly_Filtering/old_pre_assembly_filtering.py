from src01.Molecule import RCA_Molecule
from src03_Assembly_Cian.building_block_utility import *
from src02_Pre_Assembly_Filtering.constants_pre_ass_filtering import get_boxes
from ase import io
from stk import *
import time

import logging
import stk
import numpy as np
from src03_Assembly_Cian.building_block_utility import rotate_tridentate_bb, rotate_tetradentate_bb, penta_as_tetra, \
    get_optimal_rotation_angle_tridentate, Bidentate_Rotator
from src03_Assembly_Cian.stk_utils import create_placeholder_Hg_bb
import src03_Assembly_Cian.stk_extension as stk_e
from src01.Molecule import RCA_Ligand



def visualize_complex():
    mol_ = io.read('../tmp/complex.xyz')
    ase_mol = RCA_Molecule(mol=mol_)
    ase_mol.view_3d()
    os.system("rm -f ../tmp/complex.xyz")

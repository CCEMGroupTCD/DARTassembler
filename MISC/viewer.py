import math
import numpy as np
from copy import deepcopy

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return [[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]]


with open("/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/output_test/all_max.xyz") as file:
    master_list= []
    i = 0
    num_atoms = 0

    tmp_molecule_dic = {}
    for line in file:
        word_list = line.split(" ")
        print(len(word_list))
        print(word_list)
        if len(word_list) == 1:
            master_list.append(deepcopy(tmp_molecule_dic))
            tmp_molecule_dic = {}
            tmp_molecule_dic.update({"atom_num": word_list[0].split("\n")[0]})
        elif len(word_list) == 8:
            tmp_molecule_dic.update({f"atom_{i}": [word_list[0], word_list[2], word_list[4], word_list[6]]})
        elif len(word_list) == 2:
            pass
        else:
            pass


        i +=1
    print("done")




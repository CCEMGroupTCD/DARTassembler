## former known as 'write_to_file'
import pickle
import os
from tqdm import tqdm
from ASE_Molecule import xyz_file


def get_all_xyzs():
    '''
    reads the two Cambridge database files
    :return: a list of strings, which essentially splits the database files into seperate molecules
    '''
    xyzs = list()

    for i in [1, 2, 3, 4]:
        with open(f"../database/tmQM_X{i}.xyz", "r") as f:
            lines = f.readlines()

        enum_counter = 0
        num_atoms = []
        for idx, line in enumerate(lines):
            if lines[idx + 1][:5] == "CSD_c":
                num_atoms.append(int(line))
            if idx + 2 == len(lines):
                break

        for n_at in num_atoms:
            xyzs.append(lines[enum_counter:enum_counter + n_at + 3][:-1])
            # splits up the chunks in the text file as per the number of atoms
            enum_counter += n_at + 3

    return xyzs


def safe_and_reset(k, relevant_xyzs):

    with open(f"../tmp/relevant_xyzs_{k}.pickle", "wb") as handle:
        pickle.dump(relevant_xyzs, handle)

    return k+1, dict()


def merge_dicts(k):

    full_dict = dict()

    for i in range(k):
        with open(f"../tmp/relevant_xyzs_{i}.pickle", "rb") as handle:
            new = pickle.load(handle)
            full_dict.update(new)
            os.remove(f"../tmp/relevant_xyzs_{i}.pickle")

    return full_dict


def read_tmqm_db():
    """
    :return: dict_ {csd_code:xyz_class}
    """
    xyzs = get_all_xyzs()

    relevant_xyzs = dict()
    k = 0

    for j, xyz in tqdm(enumerate(xyzs)):

        # fill relevant_xyzs
        if "CSD_code" in xyz[1]:
            csd_code = xyz[1].split("CSD_code")[1][3:9]
            atom_num = int(xyz[0])

            coordinates = dict()

            for i, str_ in enumerate(xyz[2:]):
                a = str_.split()
                atom_name = a[0]
                coords = [float(el) for el in a[1:]]
                coordinates[i] = [atom_name, coords]

            relevant_xyzs[csd_code] = xyz_file(csd_code=csd_code, atom_number=atom_num, coordinates=coordinates)

        if j % 100 == 98:
            k, relevant_xyzs = safe_and_reset(k, relevant_xyzs)

    full_dict = merge_dicts(k)

    with open("../data/new_xyzs.pickle", "wb") as handle:
        pickle.dump(full_dict, handle)



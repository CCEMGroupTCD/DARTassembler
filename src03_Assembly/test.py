import pickle
import collections
from tqdm import tqdm
from sympy import Point3D, Plane
import yaml
from yaml import SafeLoader


def get_sum_formula(self):
    sum_formula_1 = [a[0] for a in self.coordinates.values()]
    return collections.Counter(sum_formula_1)


def planar_check(self):
    fc = list()

    for index, information in self.coordinates.items():
        if self.ligand_to_metal[index] == 1:
            fc.append(information[1])

    c1, c2, c3, c4 = Point3D(fc[0]), Point3D(fc[1]), Point3D(fc[2]), Point3D(fc[3])
    E = Plane(c1, c2, c3)
    if round(E.distance(c4)) < 2:
        return True

    return False


if __name__ == "__main__":

    database_path = "../data/LigandDatabases/ligand_db_w_filters.pickle"

    with open(database_path, "rb") as handle:
        ligDB = pickle.load(handle)

    print(ligDB.filtered_database.keys())
    print(ligDB.full_ligand_dict.keys())



    with open("../data/LigandDatabases/ligand_db_w_filters.pickle", "rb") as handle:
        ligand_db = pickle.load(handle)

    tetra = ligand_db.filtered_database[4]

    for j, te in tqdm(enumerate(tetra)):
        if planar_check(te) is False:
            te.view_3d()
            input("press")

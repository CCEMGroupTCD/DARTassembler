import unittest

import os
import pickle

from src.read_database import get_all_xyzs, read_local_tmqm_db
from src.LigandDatabase import LigandDatabase


class Extraction(unittest.TestCase):

    def test_get_all_xyzs(self):
        """
        the given files also display how the expected outlook from the get_all_xyzs() looks like
        """

        data = list()
        data.append(open(f"../test/test_xyz_files/WIXKOE.xyz", "r").readlines())
        data.append(open(f"../test/test_xyz_files/DUCVIG.xyz", "r").readlines())
        data.append(open(f"../test/test_xyz_files/KINJOG.xyz", "r").readlines())

        extraction = get_all_xyzs()

        for i, el in enumerate(data):
            self.assertEqual(extraction[i], el)

    def test_read_local_tmqm(self):
        with open(f"../test/other_test_files/coordinate_dict_example.pickle", "rb") as file:
            data = pickle.load(file)

        db = read_local_tmqm_db()

        for test_csd, data_coords in data.items():
            for csd_code, coords in db.items():
                if csd_code == test_csd:
                    self.assertEqual(data_coords, coords)

    def test_planar_check(self):
        """
        Important: The ground truth corresponds to this exact test_set
        """

        if not os.path.exists("other_test_files/ligand_db_test.pickle"):
            ligand_db = LigandDatabase(TestSize=True)
            ligand_db.extract_ligands(denticity_numbers_of_interest=[2, 3, 4, 5],
                                      metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"]
                                      )
            ligand_db.add_monodentate_ligands()
            pickle.dump(ligand_db, open("other_test_files/ligand_db_test.pickle", "wb"))

        with open("other_test_files/ligand_db_test.pickle", "rb") as handle:
            database = pickle.load(handle)

        operating_dict = database.full_ligand_dict
        tridentate_list = operating_dict[3]

        planar_check_list = []
        for j, tri in enumerate(tridentate_list):
            planar_check_list.append(tri.planar_check())

            # for bugfixing:
            # print(f"Is planar? {tri.planar_check()}")
            # fc = list()

            # for index, information in tri.coordinates.items():
                # if tri.ligand_to_metal[index] == 1:
                    # fc.append(information[1])

            # tri.add_atom(symbol="He", coordinates=[0.0, 0.0, 0.0])
            # for el in fc:
                # print(el)
            # tri.view_3d()
            # print("\n")
            # input("press enter to cont")

        ground_truth = [True, True, True, True, False, True, True, False, False, True, True,
                        False, False, False, True, True, True, False, True, True,
                        False, True, False, True, False, True, False, True, False,
                        False, True, True, True, True, True, True, False, True,
                        True, True, True, True, True, False, True, True, False]

        self.assertEqual(planar_check_list[:len(ground_truth)], ground_truth)


if __name__ == '__main__':
    unittest.main()

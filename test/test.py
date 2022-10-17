import unittest

import pickle

from src.read_database import get_all_xyzs, read_local_tmqm_db


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



    def test_planar_type(self):
        pass


if __name__ == '__main__':
    unittest.main()

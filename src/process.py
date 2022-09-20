import read_database as preprocessing
import os
import Extracted_Molecule as csd
import RCA_Molecule
from mendeleev import element
from ligand_filter01_duplicates import *


def extract_from_db():
    if not os.path.exists("../data/new_xyzs.pickle"):
        preprocessing.read_tmqm_db()

    with open("../data/new_xyzs.pickle", "rb") as handle_:
        all_relevant_xyzs = pickle.load(handle_)

    return all_relevant_xyzs


def create_csd_ASEMol_dict(all_relevant_xyzs):

    dict_ = {}
    for csd_key, xyz in tqdm(all_relevant_xyzs.items()):
        with open("../tmp/tmp.xyz", "w+") as text_file:
            text_file.write(xyz.get_xyz_file_format_string())
        dict_[csd_key] = RCA_Molecule.RCA_Molecule(io.read("../tmp/tmp.xyz"))

    with open("../data/csd_ase_dict.pickle", "wb") as h:
        pickle.dump(dict_, h)


def extract_ligands(xyz_dict, denticity_numbers, metals_of_interest=None):
    """
    Runs the ligand extraction
    :param xyz_dict: dictionary with the .xyz information for the molecules
    :param denticity_numbers: denticity numbers of ligands, we wish to extract
    :param metals_of_interest: center metals of interest for the extraction
    :return: dict: {dent_number : list_of_ligands_w_that_denticity}
    """
    ligand_dict_ = {dent: [] for dent in denticity_numbers}

    for xyz_file_ in tqdm(xyz_dict.values()):
        try:
            test_mol = csd.Extracted_Molecule(xyz=xyz_file_)

            if metals_of_interest is None or element(int(test_mol.original_metal)).symbol in metals_of_interest:
                test_mol.extract_ligands(denticity_numbers=denticity_numbers)

                for lig in test_mol.ligands:
                    if lig.denticity in ligand_dict_.keys():
                        ligand_dict_[lig.denticity].append(lig)

        except Exception as ex:
            print(f"An Error occured: {ex}")
            pass

    return ligand_dict_


def add_monodentate_ligand(ligand_dict_):
    # Problem: No monodentate ligands
    #
    ligand_dict_[1] = []

    #
    # add monodentate ligands
    coordinates_ = {0: ["O", [0, 0, 1.4361]], 1: ["H", [0.2096, -0.5615, 2.1227]]}

    Hydroxi = RCA_Molecule.RCA_Ligand(xyz=RCA_Molecule.xyz_file(atom_number=2,
                                                                csd_code="Hydroxi",
                                                                coordinates=coordinates_),
                                      denticity=1,
                                      ligand_to_metal=[1, 0],
                                      name="Hydroxi"
                                      )

    ligand_dict_[1].append(Hydroxi)

    return ligand_dict_


def run(metals_of_interest: list, denticity_numbers: list, TestSize=False, **kwargs):
    """
    Runs the full script for the ligand extraction
    Parameters to set:
        metals_of_interest
        denticity_numbers
    """
    #
    # 1. DB extraction
    all_relevant_xyzs = extract_from_db()

    if TestSize is True:
        small_relevant_xyzs = dict()
        for i, (key, value) in enumerate(all_relevant_xyzs.items()):
            if i % 100 == 1:
                small_relevant_xyzs[key] = value
        all_relevant_xyzs = small_relevant_xyzs

    #
    # 2. (optional) csd_ASE_dict
    if "get_csd_Ase_dict" in kwargs.keys():
        if kwargs["get_csd_Ase_dict"]:
            create_csd_ASEMol_dict(all_relevant_xyzs)

    #
    # 3. extract ligands from DB molecules
    _ligand_dict = extract_ligands(xyz_dict=all_relevant_xyzs,
                                   denticity_numbers=denticity_numbers,
                                   metals_of_interest=metals_of_interest
                                   )

    # 3 (b) Add monodontate ligand, manually
    _ligand_dict = add_monodentate_ligand(_ligand_dict)

    with open("../data/ligand_dict_test.pickle", "wb") as h:
        pickle.dump(_ligand_dict, h)

    #
    #
    #
    # Filtering (all of these are pre-assembly filters)
    if "Duplicate_filter" in kwargs:
        if "Duplicate_filter":
            _ligand_dict = duplicant_filter(_ligand_dict)
        #

    if "Box_filter" in kwargs:
        if "Box_filter":
            pass
            # todo: ligand_dict = box_filter(ligand_dict)

    with open("../data/ligand_dict_filtered.pickle", "wb") as h:
        pickle.dump(_ligand_dict, h)



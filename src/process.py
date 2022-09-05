import get_csd_xyz_dict as preprocessing
import os
import pickle
import CSD_Molecule as csd
from tqdm import tqdm
from constants import PSE
import ASE_Molecule
from ase import io


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
        dict_[csd_key] = ASE_Molecule.ASE_Molecule(io.read("../tmp/tmp.xyz"))

    with open("../data/csd_ase_dict.pickle", "wb") as handle:
        pickle.dump(dict_, handle)


def extract_ligands(xyz_dict, denticity_numbers, metals_of_interest=None):
    """
    :param xyz_dict: dictionary with the .xyz information for the molecules
    :param denticity_numbers: denticity numbers of ligands, we wish to extract
    :param metals_of_interest: center metals of interest for the extraction
    :return: dict: {dent_number : list_of_ligands_w_that_denticity}
    """
    ligand_dict = {dent: [] for dent in denticity_numbers}

    for xyz_file in tqdm(xyz_dict.values()):
        try:
            test_mol = csd.CSD_Molecule(xyz=xyz_file)

            if metals_of_interest is None or PSE[test_mol.original_metal].symbol in metals_of_interest:
                test_mol.extract_ligands(denticity_numbers=denticity_numbers)

                for lig in test_mol.ligands:
                    if lig.denticity in ligand_dict.keys():
                        ligand_dict[lig.denticity].append(lig)

        except Exception as ex:
            print(f"An Error occured: {ex}")
            pass

    return ligand_dict


def add_monodentate_ligand(ligand_dict):
    # Problem: No monodentate ligands
    #
    # add monodentate ligand:
    property_dict_ = {
        'denticity':1,
        'ligand_to_metal': [1, 0],
        'original_metal': 80,
        'name': "Hydroxy"
    }
    coordinates_ = {0: ["O", [0, 0, 1.4361]], 1: ["H", [0.2096, -0.5615, 2.1227]]}

    xyz_file_ = ASE_Molecule.xyz_file(atom_number=2, csd_code="Hydroxi", coordinates=coordinates_)

    ligand_dict[1] = [ASE_Molecule.ASE_Ligand(xyz=xyz_file_, property_dict=property_dict_)]

    return ligand_dict


def run(metals_of_interest: list, denticity_numbers: list, Testing=False, get_csd_Ase_dict=False):
    """
    Runs the full script for the ligand extraction
    Parameters to set:
        metals_of_interest
        denticity_numbers
    """
    #
    # 1. DB extraction
    all_relevant_xyzs = extract_from_db()

    if Testing is True:
        small_relevant_xyzs = dict()
        for i, (key, value) in enumerate(all_relevant_xyzs.items()):
            if i % 100 == 1:
                small_relevant_xyzs[key] = value
        all_relevant_xyzs = small_relevant_xyzs

    #
    # 2. (optional) csd_ASE_dict
    if get_csd_Ase_dict is True:
        create_csd_ASEMol_dict(all_relevant_xyzs)

    #
    # 3. extract ligands from DB molecules
    ligand_dict = extract_ligands(xyz_dict=all_relevant_xyzs,
                                  denticity_numbers=denticity_numbers,
                                  metals_of_interest=metals_of_interest
                                  )

    # 3 (b) Add monodontate ligand, manually
    ligand_dict = add_monodentate_ligand(ligand_dict)

    with open("../data/ligand_dict.pickle", "wb") as handle:
        pickle.dump(ligand_dict, handle)

    #
    #
    # todo: From here on we need to work
    #
    # 4. remove duplicates
    # todo: muss glaub auch nochmal überarbeitet werden
    # clus.cluster_ligands(dent_number=4, depth=5, target_path="data")

    #
    # 5. assemble ligands


if __name__ == '__main__':
    run(Testing=True,
        get_csd_Ase_dict=False,
        metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"],
        denticity_numbers=[2, 3, 4, 5]
        )




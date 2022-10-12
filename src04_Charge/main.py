from src.read_database import read_tmqm_db_with_partial_charge
from tqdm import tqdm
from src04_Charge.Charge_Extraction import Extracted_molecule_w_charge
import pickle


def main():
    xyz_dict = read_tmqm_db_with_partial_charge()

    example_dict = {}

    for csd_code, coordinates in tqdm(xyz_dict.items()):
        try:
            molecule = Extracted_molecule_w_charge(coordinates=coordinates, csd_code=csd_code)
            molecule.extract_ligands()
            charge_dict = molecule.charge_dict

            example_dict[csd_code] = charge_dict

        except Exception as e:
            pass

        if len(example_dict) > 25:
            break

    with open("../data/Charge/charge_dict.pickle", "wb") as handle:
        pickle.dump(example_dict, handle)


if __name__ == "__main__":
    # main()

    path_to_pickle = "../data/Charge/charge_dict.pickle"

    with open(path_to_pickle, "rb") as handle:
        dict_ = pickle.load(handle)


    print("done")

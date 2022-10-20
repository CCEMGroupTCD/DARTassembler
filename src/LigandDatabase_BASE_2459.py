from src.read_database import read_local_tmqm_db
from src.Extracted_Molecule import Extracted_Molecule
from src.Molecule import RCA_Molecule, RCA_Ligand
from src.utilities import coordinates_to_xyz_str
from src.constants import get_monodentate_list
from mendeleev import element
import pickle
from tqdm import tqdm
from ase import io


class LigandDatabase:

    def __init__(self, TestSize=False):

        self.extracted_information = read_local_tmqm_db()

        if TestSize is True:
            small_relevant_xyzs = dict()
            for i, (key, value) in enumerate(self.extracted_information.items()):
                if i % 100 == 1:
                    small_relevant_xyzs[key] = value
            self.extracted_information = small_relevant_xyzs

        self.full_ligand_dict = {}
        self.ligand_dict = {}

        self.extraction_error_count = 0

    def __len__(self):
        """
        returns total length of the ligand_dict
        """
        if self.ligand_dict == {}:
            return 0
        else:
            len_ = 0
            for key, _list in self.ligand_dict.items():
                len_ += len(_list)
            return len_

    def __str__(self):
        if self.ligand_dict == {}:
            return ""
        else:
            str_ = ""
            for key, _list in self.ligand_dict.items():
                str_ += f"For denticiy {key} we have {len(_list)} ligands\n"
            return str_

    def get_information(self):
        print(str(self))

    @staticmethod
    def read_db():
        """
        # todo: this should be generalized to a more customizable input system
        return: dict of the type {csd_code: coordinates}
        """
        return read_local_tmqm_db()

    def get_molecule_lookup_dict(self, path: str = "../data"):
        """
        for the remainder we only have the extracted ligands from the molecules in the database
        if we would like to have a lookup table, call this method
        :path: where to store the lookup table
        """
        dict_ = {}
        for csd_key, coordinates in tqdm(self.extracted_information.items()):
            with open("../tmp/tmp.xyz", "w+") as text_file:
                text_file.write(coordinates_to_xyz_str(coordinates=coordinates))
            dict_[csd_key] = RCA_Molecule(io.read("../tmp/tmp.xyz"))

        with open(f"{path}/csd_molecule_lookup_file.pickle", "wb") as h:
            pickle.dump(dict_, h)

    def extract_ligands(self, denticity_numbers_of_interest, metals_of_interest=None):

        if isinstance(metals_of_interest, str):
            metals_of_interest = [metals_of_interest]

        if isinstance(denticity_numbers_of_interest, int):
            denticity_numbers_of_interest = [denticity_numbers_of_interest]

        #
        # init empty dict
        self.full_ligand_dict = {dent: [] for dent in denticity_numbers_of_interest}

        for csd_code, coordinates in tqdm(self.extracted_information.items()):
            try:
                molecule = Extracted_Molecule(coordinates=coordinates, csd_code=csd_code)

                if metals_of_interest is None or element(int(molecule.original_metal)).symbol in metals_of_interest:
                    # now we extract the ligands of the desired molecule
                    molecule.extract_ligands(denticity_numbers=denticity_numbers_of_interest)

                    # and add the ligands to our dict, if there are any
                    for lig in molecule.ligands:
                        if lig.denticity in self.full_ligand_dict.keys():
                            self.full_ligand_dict[lig.denticity].append(lig)

            except Exception as ex:
                print(f"An Error occured: {ex}")
                self.extraction_error_count += 1
                pass

        print(f"Extraction complete -- number of errors {self.extraction_error_count}")

    def add_monodentate_ligands(self):
        """
        problem is that we don't extract monodentates; so we just insert them customly from pre-defined constants
        """
        self.ligand_dict[1] = get_monodentate_list()

    def filter_duplicates(self):
        for denticity, ligand_list in self.full_ligand_dict.items():
            self.ligand_dict[denticity] = list(set(ligand_list))






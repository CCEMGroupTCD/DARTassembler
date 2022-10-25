from src.read_database import read_local_tmqm_db
from src.Extracted_Molecule import Extracted_Molecule
from src.Molecule import RCA_Molecule, RCA_Ligand
from src.utilities import coordinates_to_xyz_str
from src.constants import get_monodentate_list
from mendeleev import element
import pickle
from tqdm import tqdm
from ase import io
import numpy as np
import pandas as pd
from src.MoleculeDatabase import MoleculeDatabase
import itertools
from copy import deepcopy

def iter_over_ligand_dict(ligand_dict: dict, return_denticity: bool=False):
    for denticity, ligand_list in ligand_dict.items():
        for ligand in ligand_list:
            if return_denticity:
                yield ligand, denticity
            else:
                yield ligand
def total_number_of_ligands(ligand_dict: dict):
    num = 0
    for denticity, ligand_list in ligand_dict.items():
        num += len(ligand_list)
    return num
class LigandDatabase(MoleculeDatabase):

    def __init__(self, data_path, id_col, TestSize=False):
        MoleculeDatabase.__init__(self, data_path=data_path, id_col=id_col, TestSize=TestSize)
        # self.extracted_information = read_local_tmqm_db()
        #
        # if TestSize is True:
        #     small_relevant_xyzs = dict()
        #     for i, (key, value) in enumerate(self.extracted_information.items()):
        #         if i % 100 == 1:
        #             small_relevant_xyzs[key] = value
        #     self.extracted_information = small_relevant_xyzs

        self.full_ligand_dict = {}
        self.ligand_dict = {}

        self.extraction_error_count = 0

        # exclusively for later filter stages
        self.filtered_database_dict = None
        self.filtered_database = None

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

    # def get_molecule_lookup_dict(self, path: str = "../data"):
    #     """
    #     for the remainder we only have the extracted ligands from the molecules in the database
    #     if we would like to have a lookup table, call this method
    #     :path: where to store the lookup table
    #     """
    #     dict_ = {}
    #     for csd_key, coordinates in tqdm(self.extracted_information.items()):
    #         with open("../tmp/tmp.xyz", "w+") as text_file:
    #             text_file.write(coordinates_to_xyz_str(coordinates=coordinates))
    #         dict_[csd_key] = RCA_Molecule(io.read("../tmp/tmp.xyz"))
    #
    #     with open(f"{path}/csd_molecule_lookup_file.pickle", "wb") as h:
    #         pickle.dump(dict_, h)

    def extract_ligands(self, denticity_numbers_of_interest, metals_of_interest=None):
        
        if isinstance(metals_of_interest, str):
            metals_of_interest = [metals_of_interest]

        if isinstance(denticity_numbers_of_interest, int):
            denticity_numbers_of_interest = [denticity_numbers_of_interest]

        #
        # init empty dict
        self.full_ligand_dict = {dent: [] for dent in denticity_numbers_of_interest}
        self.get_all_Extracted_Molecules()

        for csd_code, molecule in tqdm(self.all_Extracted_Molecules.items(), desc='Extracting ligands'):
            try:
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
                self.extraction_errors[csd_code] = str(ex)
                pass

        print(f"Extraction complete -- number of errors {self.extraction_error_count}")

    def add_monodentate_ligands(self):
        """
        problem is that we don't extract monodentates; so we just insert them customly from pre-defined constants
        """
        self.ligand_dict[1] = get_monodentate_list()
    
    def get_all_denticities(self, ligand_dict):
        denticities = [denticity for denticity in ligand_dict.keys()]
        return denticities
    def filter_duplicates(self):
        print('Start filtering duplicates.')
        all_ligands = list(iter_over_ligand_dict(self.full_ligand_dict))
        
        unique_hashed_ligands = list(set(all_ligands))
        print(f'Number of unique hashed ligands: {len(unique_hashed_ligands)}.')
        
        all_hashes = [ligand.graph_hash for ligand in unique_hashed_ligands]
        all_ligands_by_hashes = {h: [] for h in all_hashes}
        for ligand in all_ligands:
            all_ligands_by_hashes[ligand.graph_hash].append(ligand)
        
        unique_ligands = deepcopy(unique_hashed_ligands)
        counter = 1
        for i, ligand1 in enumerate(unique_hashed_ligands):
            same_hashed_ligands = all_ligands_by_hashes[ligand1.graph_hash]
            
            for ligand2 in same_hashed_ligands:
                exactly_same_ligand = ligand1.name == ligand2.name
                if exactly_same_ligand:
                    continue
    
                counter += 1
                if ligand1 != ligand2:
                    print('Found different ligands with same hashes.')
                    unique_ligands.append(ligand2)
            print(f'{i}: {counter}')
        
        self.unique_ligand_dict = {denticity: [] for denticity in self.get_all_denticities(self.full_ligand_dict)}
        for unique_ligand in unique_ligands:
            unique_ligand.unique_name = 'unq_' + unique_ligand.name
            self.unique_ligand_dict[unique_ligand.denticity].append(unique_ligand)
            
        print(f'Number of unique ligands: {len(unique_ligands)}.')
        return self.unique_ligand_dict
    
    def get_df_of_all_ligands(self):
        """
        Returns a dataframe with name, denticity, CSD code and type of every ligand in the database.
        """
        ligand_props = []
        for ligand_list in self.full_ligand_dict.values():
            for l in ligand_list:
                ligand_props.append({
                    'name': l.name,
                    'csd_code': l.csd_code if hasattr(l, 'csd_code') else np.nan,
                    'original_metal_symbol': l.original_metal_symbol if hasattr(l, 'original_metal_symbol') else np.nan,
                    'denticity': l.denticity
                    })
        df = pd.DataFrame(ligand_props)
        return df

if __name__ == '__main__':

    ligand_db_file = "../data/LigandDatabases/ligand_db_test.pickle"
    
    with open(ligand_db_file, 'rb') as file:
        ligand_db = pickle.load(file)
        print('Loaded ligand db from pickle.')
    
    ligand_db.filter_duplicates()
    
    n_unique_ligands = total_number_of_ligands(ligand_db.unique_ligand_dict)
    n_total_ligands = total_number_of_ligands(ligand_db.full_ligand_dict)
    print(f'There are {n_unique_ligands} unique ligands out of a total of {n_total_ligands} ligands.')
    
    




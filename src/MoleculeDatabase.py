"""
Takes as input the tmqm downloaded data and makes them easier accessible.
"""
import pandas as pd
from pathlib import Path
import glob
from copy import deepcopy
import io
import numpy as np
from tqdm import tqdm
from Molecule import RCA_Molecule, RCA_Ligand
from Extracted_Molecule import Extracted_Molecule
import json

ATOMIC_PROPERTIES_COLUMN_SEPARATOR = '  ===  '

class MoleculeDatabase:
    
    def __init__(self, data_path: str, id_col: str, TestSize=False):
        self.data_path = data_path
        
        self.global_props_path = Path(self.data_path, 'global_mol_properties.csv')
        self.atomic_props_dir = Path(self.data_path, 'atomic_properties')
        
        self.id_col = id_col
        self.test_size = TestSize
        
        self.has_df = False
        self.has_atomic_props = False
    def get_global_properties(self) -> pd.DataFrame:
        """
        Returns a deepcopy of the internal global molecular properties dataframe.
        :return:
        """
        self.df = pd.read_csv(self.global_props_path)
        self.has_df = True
        return deepcopy(self.df)
    
    def get_all_molecular_ids(self) -> np.array:
        """
        Returns all molecular ids. The id is a unique string for each molecule, for example the CSD code.
        :return: molecular ids (np.array)
        """
        if not self.has_df:
            self.get_global_properties()
            
        ids = self.df[self.id_col]
        assert ids.is_unique, 'Molecular IDs are not unique.'
        ids = ids.values
        
        return ids
    
    def read_atomic_properties_file_header(self, filestr: str):
        lines = filestr.split('\n')
        num_atoms = int(lines[0])
    
        columns = lines[1].split(ATOMIC_PROPERTIES_COLUMN_SEPARATOR)
        mol_id = columns[0]
        comment = columns[-1]
        col_names = columns[1:-1]
        
        return num_atoms, mol_id, col_names, comment
        
    def read_single_atomic_properties_filestring(self, filestr: str,):
        num_atoms, mol_id, col_names, comment = self.read_atomic_properties_file_header(filestr)
        
        file = io.StringIO(filestr)
        txt_array = np.genfromtxt(file, skip_header=2, dtype=str)
        atoms, value_array = txt_array[:,0], txt_array[:,1:].astype(float)
        
        col_names.pop(0)
        values = {name: column.tolist() for name, column in zip(col_names, value_array.T)}
        atoms = [str(a) for a in atoms]     # numpy to regular python list of strings for json compatibility
        
        assert num_atoms == len(value_array) and num_atoms == len(atoms), 'Number of atoms specified in atomic properties file is not equal to the real number of atoms included.'
        return mol_id, atoms, values, comment
        
    def read_full_atomic_properties_file(self, path, sep='\n\n'):
        with open(path, 'r') as file:
            full_file = file.read()
        files = full_file.split(sep)
        
        atomic_props = {}
        for filestr in tqdm(files):
            mol_id, atoms, values, comment = self.read_single_atomic_properties_filestring(filestr)
            assert not mol_id in atomic_props, f'Molecular id {mol_id} not unique in file{path}.'
            values['atoms'] = atoms
            values['comment'] = comment
            atomic_props[mol_id] = values
        
        return atomic_props
        
    def get_atomic_properties(self) -> dict:
        """
        Returns a dictionary of all atomic properties of all molecules. The form is {mol_id: {'prop1': list}}.
        :return: atomic properties (dict)
        """
        all_atomic_props = {}
        
        pattern = Path(self.atomic_props_dir, '*.xyz')
        self.all_atomic_props_paths = sorted(glob.glob(str(pattern)))
        print(f'Found {len(self.all_atomic_props_paths)} atomic property files. Start reading in.')
        
        for atm_path in self.all_atomic_props_paths:
            atomic_props = self.read_full_atomic_properties_file(path=atm_path)
            
            assert not any([mol_id in all_atomic_props for mol_id in atomic_props.keys()]), 'Molecular ids of atomic property files not unique.'
            all_atomic_props.update(atomic_props)
        
        self.atomic_props = all_atomic_props
        self.has_atomic_props = True
        return all_atomic_props
    
    def save_atomic_properties(self, outpath=None):
        if not self.has_atomic_props:
            self.get_atomic_properties()
        
        self.atomic_properties_json = outpath or Path(self.atomic_props_dir, 'atomic_properties.json')
        with open(self.atomic_properties_json, 'w') as file:
            json.dump(deepcopy(self.atomic_props), file)
            
        print(f'Saved {self.atomic_properties_json}.')
        return
    
    def load_atomic_properties(self, json_path):
        with open(json_path, 'r') as file:
            self.atomic_props = json.load(file)
        
        return self.atomic_props
    
    # def get_Extracted_Molecule(self, mol_id) -> Extracted_Molecule:
    #     if not self.has_atomic_props:
    #         self.get_atomic_properties()
    #     atm = self.atomic_props[mol_id]
    #
    #
    #     return Extracted_Molecule()
    
    
        
            
        

if __name__ == '__main__':
    
    id_col = 'CSD_code'
    
    tmqm = MoleculeDatabase(data_path='../database/tmQM/data', id_col=id_col)
    # df = tmqm.get_global_properties()
    # atm = tmqm.get_atomic_properties()
    # ids = tmqm.get_all_molecular_ids()
    # mol = tmqm.get_Extracted_Molecule('WIXKOE')
    #tmqm.save_atomic_properties()
    print('Done!')

        
        
        

        
        
        


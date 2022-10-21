"""
Takes as input the tmqm downloaded data and makes them easier accessible.
"""
import pandas as pd
from pathlib import Path
from pymatgen.core import Composition

class MoleculeDatabase():
    
    def __init__(self, tmqm_dir_path):
        self.tmqm_dir_path = tmqm_dir_path
        self.tmqm_raw_data_path = Path(self.tmqm_dir_path, 'raw_data')
        
        self.data_path = Path(self.tmqm_dir_path, 'data')
        self.xyz_files_dir = Path(self.data_path, 'xyz_files')
        self.save_full_df_of_global_properties = Path(self.data_path, 'molecules.csv')

        self.all_xyz_path = Path(self.tmqm_raw_data_path, 'tmQM_X.xyz')
        self.all_properties_path = Path(self.tmqm_raw_data_path, 'tmQM_y.csv')

    def extract_and_save_each_xyz_file(self):
        with open(self.all_xyz_path, 'r') as file:
            full_file = file.read()
            all_xyz_files = full_file.split('\n\n')
            all_xyz_files = all_xyz_files[:-1]
            
        self.xyz_files_dir.mkdir(parents=True, exist_ok=True)
    
        all_molecule_infos = []
        for xyz_file_str in all_xyz_files:
            results = get_molecule_info_from_xyz_file(xyz_file_str)
            all_molecule_infos.append(results)
        
            CSD_code = results['CSD_code']
            outpath = Path(self.xyz_files_dir, f'{CSD_code}.xyz')
            with open(outpath, 'w') as outfile:
                outfile.write(xyz_file_str)
        
        df_all_molecules_from_xyz = pd.DataFrame(all_molecule_infos)
        
        return df_all_molecules_from_xyz
    
    def get_metal_from_stoichiometry(self, stoichiometry: str) -> str:
        comp = Composition(stoichiometry)
        metal = [str(el) for el in comp if el.is_metal]
        assert len(metal) == 1
        metal = metal[0]
        return metal
    
    def make_data_better_available(self):
        df_all_molecules_from_xyz = self.extract_and_save_each_xyz_file()
        
        df_all_properties = pd.read_csv(self.all_properties_path, delimiter=';')
        df = pd.merge(df_all_molecules_from_xyz, df_all_properties, on='CSD_code')
        
        df['stoichiometry'] = df['stoichiometry'].str.replace('\(.*\)$', '', regex=True) # remove charge string
        
        df['metal'] = df['stoichiometry'].apply(self.get_metal_from_stoichiometry)
        
        df.to_csv(self.save_full_df_of_global_properties, index=False)
        
        return df
        

if __name__ == '__main__':
    
    tmqm_dir_path = '/home/timo/PhD/projects/RCA/projects/LigandCharges/data/tmqm'
    
    tmqm = tmQM(tmqm_dir_path)
    df = tmqm.make_data_better_available()

    

        
        
        

        
        
        


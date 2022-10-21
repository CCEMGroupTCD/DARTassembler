"""
Takes as input the tmqm downloaded data and makes them easier accessible.
"""
import os
import pandas as pd
from pathlib import Path
from pymatgen.core import Composition
import numpy as np
import io
from tqdm import tqdm

def get_molecule_info_from_tmQM_xyz_file(xyz_file_str: str):
    """Parses the first two lines of the xyz file to return multiple global molecular properties.
    """
    lines = xyz_file_str.split('\n')
    num_atoms = int(lines[0])
    
    title_infos = lines[1].split('|')
    assert len(title_infos) == 5
    
    CSD_code = title_infos[0]
    assert CSD_code.startswith('CSD_code = ') and CSD_code.endswith(' ')
    CSD_code = CSD_code[11:-1]
    
    total_charge = title_infos[1]
    assert total_charge.startswith(' q = ') and total_charge.endswith(' ')
    total_charge = int(total_charge[5:-1])
    
    spin_mult = title_infos[2]
    assert spin_mult.startswith(' S = ') and spin_mult.endswith(' ')
    spin_mult = int(spin_mult[5:-1])
    
    stoichiometry = title_infos[3]
    assert stoichiometry.startswith(' Stoichiometry = ') and stoichiometry.endswith(' ')
    stoichiometry = stoichiometry[17:-1]
    
    metal_node_degree = title_infos[4]
    assert metal_node_degree.startswith(' MND = ')
    metal_node_degree = int(metal_node_degree[7:])
    
    results = {'CSD_code': CSD_code,
               'stoichiometry': stoichiometry,
               'total_charge': total_charge,
               'num_atoms': num_atoms,
               'metal_node_degree': metal_node_degree,
               'spin_multiplicity': spin_mult
               }
    
    return results

def break_file_into_chunks(single_file_path: str, save_dir_path: str, n_chunks: int, out_filename: str=None, break_on: str='\n\n') -> None:
    """
    Breaks a single text file with many subfiles (e.g. many xyz files) into multiple smaller files for github compatibility.
    :param single_file_path: Path to
    :param save_dir_path:
    :param n_chunks:
    :param out_filename:
    :param break_on:
    :return: None
    """
    single_file_path = Path(single_file_path)
    filename = single_file_path.name
    save_dir_path = Path(save_dir_path)
    save_dir_path.mkdir(parents=True, exist_ok=True)
    
    with open(single_file_path, 'r') as file:
        full_file = file.read()
        all_files = full_file.split(break_on)
        if all_files[-1] == '':
            all_files = all_files[:-1]
        
    n_saved_files = 0     # assertion test
    chunks = np.array_split(all_files, n_chunks)
    for i, chunk_files in enumerate(chunks):
        chunk_files_str = break_on.join(chunk_files)
        
        chunk_filename = f'{i}_' + filename if out_filename is None else f'{i}_' + out_filename
        save_file_path = Path(save_dir_path, chunk_filename)
        with open(save_file_path, 'w') as file:
            file.write(chunk_files_str)
        
        n_saved_files += len(chunk_files)
    assert n_saved_files == len(all_files)
    
    return

def array2tabularstr(array):
    return '\n'.join(['\t'.join(col) for col in array])

def get_atomic_properties_filestr(num_atoms, mol_id, column_names, values, sep='  ===  '):
    comment = ''
    column_names = [mol_id] + column_names + [comment]
    column_names_str = sep.join(column_names)
    header = f'{num_atoms}' + '\n' + column_names_str
    
    assert num_atoms == len(values)
    new_file_str = array2tabularstr(values)
    new_file_str = header + '\n' + new_file_str
    return new_file_str

def combine_atomic_property_files(in_file1: str, in_file2: str, out_file_path: str, column_names: list, num_atoms_and_mol_id_func, sep1='\n\n', sep2='\n\n', skip_header1=0, skip_header2=0, skip_footer1=0, skip_footer2=0, in_column_sep=None):
    """
    Combines several files with atomic properties into one file. Assumes that they all atoms and all molecules are sorted in exactly the same way. Currently still tmQM hardcoded (see comment).
    """
    atomlabel_col = 0
    with open(in_file1, 'r') as file1:
        with open(in_file2, 'r') as file2:
            files1 = file1.read().split(sep1)
            files2 = file2.read().split(sep2)
            
            last_file_empty = files1[-1] == ''
            if last_file_empty:
                assert files2[-1] == '', f'files1[-1]: {files1[-1]}, files2[-1]: {files2[-1]}'
                files1 = files1[:-1]
                files2 = files2[:-1]
            
            n_files = len(files1)
            assert n_files == len(files2)
            
            all_combined_files = []
            for f1, f2 in tqdm(zip(files1, files2)):
                header1 = '\n'.join(f1.split('\n')[0:skip_header1])
                header2 = '\n'.join(f2.split('\n')[0:skip_header2])
                assert header2 in header1, 'In the tmpQM, all headers of the charge file should be in the xyz, this is not given here.'
                
                columns1 = np.genfromtxt(io.StringIO(f1), skip_header=skip_header1, skip_footer=skip_footer1, delimiter=in_column_sep, dtype=str)
                columns2 = np.genfromtxt(io.StringIO(f2), skip_header=skip_header2, skip_footer=skip_footer2, delimiter=in_column_sep, dtype=str)
                
                assert len(columns1) == len(columns2), 'Numbers of atoms don\'t match. Maybe the order of molecules is not the same in both files?'
                assert all(columns1[:, atomlabel_col] == columns2[:, atomlabel_col]), 'Atom labels don\'t match. Maybe the order of atoms is not the same in both molecules?'
                
                columns2 = columns2[:,atomlabel_col+1:]
                all_columns = np.concatenate([columns1, columns2], axis=1)
                
                num_atoms, mol_id = num_atoms_and_mol_id_func(f1, f2)
                new_filestr = get_atomic_properties_filestr(
                                                            num_atoms=num_atoms,
                                                            mol_id=mol_id,
                                                            column_names=column_names,
                                                            values=all_columns
                                                            )
                all_combined_files.append(new_filestr)
                
    total_file = sep1.join(all_combined_files)
    with open(out_file_path, 'w') as outfile:
        outfile.write(total_file)
        
    return
    
class tmQM():
    
    def __init__(self, tmqm_dir_path):
        self.tmqm_dir_path = tmqm_dir_path
        self.tmqm_raw_data_path = Path(self.tmqm_dir_path, 'raw_data')
        
        self.data_path = Path(self.tmqm_dir_path, 'data')
        self.xyz_files_dir = Path(self.data_path, 'xyz_files')
        self.save_full_df_of_global_properties = Path(self.data_path, 'global_mol_properties.csv')

        self.all_xyz_path = Path(self.tmqm_raw_data_path, 'tmQM_X.xyz')
        self.all_properties_path = Path(self.tmqm_raw_data_path, 'tmQM_y.csv')
        self.all_bonds_path = Path(self.tmqm_raw_data_path, 'tmQM_X.BO')
        self.all_charges_path = Path(self.tmqm_raw_data_path, 'tmQM_X.q')
        
        self.atomic_properties_dir_path = Path(self.data_path, 'atomic_properties')
        self.atomic_properties_path = Path(self.atomic_properties_dir_path, 'atomic_properties.xyz')
    
    def extract_molecule_infos_from_xyz_files(self):
        with open(self.all_xyz_path, 'r') as file:
            full_file = file.read()
            all_xyz_files = full_file.split('\n\n')
            all_xyz_files = all_xyz_files[:-1]
    
        all_molecule_infos = []
        for xyz_file_str in all_xyz_files:
            results = get_molecule_info_from_tmQM_xyz_file(xyz_file_str)
            all_molecule_infos.append(results)
        
        df_all_molecules_from_xyz = pd.DataFrame(all_molecule_infos)
        
        return df_all_molecules_from_xyz
    def extract_and_save_each_xyz_file(self, batch_size=21666):
        self.xyz_files_dir.mkdir(parents=True, exist_ok=True)

        with open(self.all_xyz_path, 'r') as file:
            full_file = file.read()
            all_xyz_files = full_file.split('\n\n')
            all_xyz_files = all_xyz_files[:-1]
        
        all_molecule_infos = []
        for xyz_file_str in all_xyz_files:
            results = get_molecule_info_from_tmQM_xyz_file(xyz_file_str)
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

    def tmQM_num_atoms_and_mol_id_func(self, xyz_file_str: str, q_file_str: str):
        infos = get_molecule_info_from_tmQM_xyz_file(xyz_file_str)
        return infos['num_atoms'], infos['CSD_code']
    def make_data_better_available(self):
        print('Extract global molecular information from xyz files.')
        df_all_molecules_from_xyz = self.extract_molecule_infos_from_xyz_files()
        
        print('Merge with calculated global electronic information in a big dataframe.')
        df_all_properties = pd.read_csv(self.all_properties_path, delimiter=';')
        df = pd.merge(df_all_molecules_from_xyz, df_all_properties, on='CSD_code')
        
        # Sanitizing
        df['stoichiometry'] = df['stoichiometry'].str.replace('\(.*\)$', '', regex=True) # remove charge string
        df['metal'] = df['stoichiometry'].apply(self.get_metal_from_stoichiometry)
        
        print(f'Save global molecular information dataframe as {self.save_full_df_of_global_properties.name}.')
        self.data_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(self.save_full_df_of_global_properties, index=False)
        
        print('Merge xyz file and charge file into atomic properties file.')
        xyz_n_comment_rows = 2
        xyz_skip_footer = 0
        q_n_comment_rows = 1
        q_skip_footer = 1
        column_names = ['atom', 'x', 'y', 'z', 'partial_charge']
        self.atomic_properties_dir_path.mkdir(parents=True, exist_ok=True)
        combine_atomic_property_files(
                                        in_file1=self.all_xyz_path,
                                        in_file2=self.all_charges_path,
                                        out_file_path=self.atomic_properties_path,
                                        skip_header1=xyz_n_comment_rows,
                                        skip_header2=q_n_comment_rows,
                                        skip_footer1=xyz_skip_footer,
                                        skip_footer2=q_skip_footer,
                                        column_names=column_names,
                                        num_atoms_and_mol_id_func=self.tmQM_num_atoms_and_mol_id_func
                                        )
        
        n_chunks = 5
        print(f'Break large atomic properties file into {n_chunks} smaller chunks for github compatibility.')
        break_file_into_chunks(
                                single_file_path=self.atomic_properties_path,
                                save_dir_path=self.atomic_properties_dir_path,
                                n_chunks=n_chunks
                                )
        os.remove(self.atomic_properties_path)
        return df
        

if __name__ == '__main__':
    
    tmqm_dir_path = '/home/timo/PhD/projects/RCA/projects/CreateTMC/database/tmQM'
    
    tmqm = tmQM(tmqm_dir_path)
    df = tmqm.make_data_better_available()
    
    print('Done!')

    

        
        
        

        
        
        


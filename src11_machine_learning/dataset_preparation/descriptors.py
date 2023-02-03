"""
This file contains classes to calculate different ML descriptors of molecules.
"""
import warnings
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
import itertools
import ase
from tqdm import tqdm
from molSimplify.Informatics.autocorrelation import generate_full_complex_autocorrelations
from molSimplify.Classes.mol3D import mol3D
import RACS
from src01.utilities import flatten_list, identify_metal_in_atoms_list


def warn_if_nan_values(df):
    nan_columns = df.columns[df.isna().any()].tolist()
    n_nan = len(nan_columns)
    if n_nan > 0:
        warnings.warn(f'{n_nan} of {len(df.columns)} calculated features have NaN values for some molecules!')
# class RAC:
#     def __init__(self, molecules, depth=4, verbose=False, ff_optimization=False, steps=100):
#         self.molecules = molecules
#         self.depth = depth
#         self.verbose = verbose
#         self.ff_optimimization = ff_optimization
#         self.steps = steps
#     def input_to_mol3D(self, input):
#         mol = mol3D()
#         if input.endswith('.xyz'):
#             mol.readfromxyz(input)
#         else:
#             ff = 'uff' if self.ff_optimimization else False
#             mol.read_smiles(input, ff=ff, steps=self.steps)
#         return mol
#     def compute_molsimplify_graph_descriptors(self) -> list:
#         all_descriptors = []
#         for mol in tqdm(self.molecules):
#
#             mol = self.input_to_mol3D(mol)
#             desc = generate_full_complex_autocorrelations(mol,depth=self.depth,loud=self.verbose, oct=False)
#
#             # Doublecheck that columns are always in same order
#             try:
#                 assert all(np.array(colnames) == np.array(list(itertools.chain.from_iterable(desc['colnames']))))
#             except NameError:
#                 pass
#
#             colnames = list(itertools.chain.from_iterable(desc['colnames']))
#             results = list(itertools.chain.from_iterable(desc['results']))
#             all_descriptors.append(results)
#
#         headers = [f'RAC_{colname}' for colname in colnames]
#         df = pd.DataFrame(all_descriptors, columns=headers)
#         warn_if_nan_values(df)
#
#         return df
#
#     def compute_own_graph_descriptors(self):
#         all_descriptors = []
#         for mol in tqdm(self.molecules):
#             #A, _ = RACS.smiles_to_graph(mol)
#             RAC_dict = RACS.compute_full_RACs(smiles_str_=mol, max_distance=self.depth)
#             all_descriptors.append(RAC_dict)
#
#
#     def compute_graph_descriptors(self, method):
#         if method == 'molsimplify':
#             df = self.compute_molsimplify_graph_descriptors()
#         elif method == 'own':
#             df = self.compute_own_graph_descriptors()
#
#         return df

        
    
class RDKit_2D:
    def __init__(self, smiles):
        # self.smiles = smiles
        # self.errors = []
        # self.mols = []
        # for sm in smiles:
        #     sio = sys.stderr = StringIO()
        #     mol = Chem.MolFromSmiles(sm)
        #     error = sio.getvalue()
        #     error = error if not error == '' else np.nan
        #     self.errors.append(error)
        #
        # self.df_err = pd.DataFrame({'isosmiles': self.smiles, 'errors': self.errors})
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles
    
    def compute_2Drdkit(self):
        rdkit_2d_desc = []
        calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        for mol in tqdm(self.mols):
            ds = calc.CalcDescriptors(mol)
            rdkit_2d_desc.append(ds)
            
        header = [f'rdkit_{desc}' for desc in calc.GetDescriptorNames()]
        df = pd.DataFrame(rdkit_2d_desc, columns=header)
        warn_if_nan_values(df)
        
        return df

class SOAP_3D():

    def __init__(self, ase_molecules: list, r_cut: float=10.0, n_max: int=3, l_max: int=2, crossover: bool=False, average: str='off'):
        from dscribe.descriptors import SOAP

        self.molecules = ase_molecules
        self.species = self.get_all_species_from_molecules()

        self.soap = SOAP(
            species=self.species,
            periodic=False,
            r_cut=r_cut,
            n_max=n_max,
            l_max=l_max,
            crossover=crossover,
            average=average
            )

    def get_all_species_from_molecules(self):
        self.all_atoms_list = []
        for mol in self.molecules:
            self.all_atoms_list.append(list(mol.symbols))

        species = np.unique(flatten_list(self.all_atoms_list)).tolist()

        return species

    def calculate_descriptors(self, only_from_metal=False):
        if only_from_metal:
            positions = [[atoms.index(identify_metal_in_atoms_list(atoms))] for atoms in self.all_atoms_list]
        else:
            positions = None

        desc = self.soap.create(system=self.molecules, positions=positions)

        return desc




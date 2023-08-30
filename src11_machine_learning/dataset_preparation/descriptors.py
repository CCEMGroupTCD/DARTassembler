"""
This file contains classes to calculate different ML descriptors of molecules.
"""
import warnings
warnings.filterwarnings("ignore", message='"import openbabel" is deprecated, instead use "from openbabel import openbabel"')
from typing import Union, Tuple
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
import itertools
import ase
from pymatgen.core.periodic_table import Element as Pymatgen_Element
import sys
from io import StringIO
from tqdm import tqdm
# import RACS
from src01.utilities import flatten_list, identify_metal_in_atoms_list
from src01.utilities_graph import get_reindexed_graph, get_sorted_atoms_and_indices_from_graph
from src01.bond_orders import graph_to_smiles

def warn_if_nan_values(df):
    nan_columns = df.columns[df.isna().any()].tolist()
    n_nan = len(nan_columns)
    if n_nan > 0:
        warnings.warn(f'{n_nan} of {len(df.columns)} calculated features have NaN values for some molecules!')

class RAC:
    def __init__(self, depth=4, molecular_stats: list[str]=None, atom_stats: list[str]=None, element_label: str='node_label'):
        self.depth = depth
        self.element_label = element_label
        self.molecular_stats = molecular_stats or ['sum', 'std', 'min', 'max']
        self.atom_stats = atom_stats or ['sum']



    def input_molecule_to_graph(self, mol) -> nx.Graph:
        """
        Convert input molecule to graph. Input molecule can be a graph or a RCA_Molecule object.
        """
        if isinstance(mol, nx.Graph):
            # mol is already a graph
            graph = mol
        else:
            try:
                # mol can be made to graph
                graph = nx.Graph(mol)
            except TypeError:
                try:
                    graph = nx.Graph(mol.graph)
                except (AttributeError, TypeError):
                    raise ValueError(f'Could not convert input to graph: {mol}')

        graph = get_reindexed_graph(graph)

        # Check if graph has node labels
        if self.element_label not in graph.nodes[0]:
            raise ValueError(f'Could not find node labels in graph. Expected label: {self.element_label}')

        return graph

    def compute_own_graph_descriptors(self, molecules: list) -> pd.DataFrame:
        all_descriptors = []
        for mol in tqdm(molecules, desc='Computing RACs'):
            features, labels = self.molecule_autocorrelation(mol=mol, return_labels=True)
            RAC_dict = {label: feature for label, feature in zip(labels, features)}
            all_descriptors.append(RAC_dict)
        df = pd.DataFrame(all_descriptors)
        warn_if_nan_values(df)

        return df

    def atom_autocorrelation(self, graph, atom_index, prop_vector, stats=None):
        """
        Calculates the atom autocorrelation for a single atom in a molecule.

        Parameters
        ----------
        graph : networkx.Graph
            The molecule graph.
        atom_index : int
            The index of the atom.
        atom_props : str
            The atom properties to use for the autocorrelation.

        Returns
        -------
        autocorr : np.ndarray
            The autocorrelation matrix for the atom. The matrix has shape (depth + 1, n_stats), where n_stats is the number of statistics calculated.
        """
        if stats is None:
            stats = self.atom_stats
        autocorr = np.zeros((self.depth + 1, len(stats)))  # Initialize autocorrelation matrix

        for d in range(0, self.depth + 1):
            neighbors = nx.descendants_at_distance(graph, source=atom_index, distance=d)
            convolution_features = [prop_vector[atom_index] * prop_vector[neighbor] for neighbor in neighbors]
            for i, stat in enumerate(stats):
                if len(convolution_features) > 0:
                    autocorr[d,i] += getattr(np, stat)(convolution_features)
                else:
                    # If there are no neighbors at this distance just initialize to 0
                    autocorr[d,i] = 0

        return autocorr

    def get_prop_vector(self, graph, prop):
        """
        Returns the property matrix for a molecule. The property matrix is a matrix where each row corresponds to an atom and each column to a property.
        """
        elements, indices = get_sorted_atoms_and_indices_from_graph(graph, atom_label=self.element_label)
        elements = [Pymatgen_Element(element) for element in elements]
        valid_props = ['electronegativity', 'row', 'group', 'atomic_mass', 'electron_affinity', 'min_oxidation_state', 'max_oxidation_state', 'ionization_energy', 'nuclear_charge', 'ident', 'topology', 'size']

        if prop == 'electronegativity':
            prop_vector = [el.X for el in elements]
            label = 'chi'
        elif prop == 'row':
            prop_vector = [el.row for el in elements]
            label = 'R'
        elif prop == 'group':
            prop_vector = [el.group for el in elements]
            label = 'G'
        elif prop == 'atomic_mass':
            prop_vector = [el.atomic_mass for el in elements]
            label = 'M'
        elif prop == 'electron_affinity':
            prop_vector = [el.electron_affinity for el in elements]
            label = 'EA'
        elif prop == 'min_oxidation_state':
            prop_vector = [el.min_oxidation_state for el in elements]
            label = 'minOS'
        elif prop == 'max_oxidation_state':
            prop_vector = [el.max_oxidation_state for el in elements]
            label = 'maxOS'
        elif prop == 'ionization_energy':
            prop_vector = [el.ionization_energy for el in elements]
            label = 'IE'
        elif prop == 'nuclear_charge':
            prop_vector = [el.Z for el in elements]
            label = 'Z'
        elif prop == 'ident':
            prop_vector = [1 for el in elements]
            label = 'I'
        elif prop == 'topology':
            prop_vector = [graph.degree[atom] for atom in indices]
            label = 'T'
        elif prop == 'size':
            prop_vector = [el.atomic_radius for el in elements]
            label = 'S'
        else:
            raise ValueError(f'Invalid property: {prop}. Valid properties are: {valid_props}')

        return prop_vector, label

    def molecule_autocorrelation(self, mol, properties: list[str]=None, return_labels: bool=False) -> Union[np.array, Tuple[np.array, list]]:
        """
        Calculates the molecule features.

        Parameters
        ----------
        graph : networkx.Graph
            The molecule graph.

        Returns
        -------
        features : np.array
            The molecule features.
        labels : list
            The labels for the features.
        """
        graph = self.input_molecule_to_graph(mol)
        default_props = ['electronegativity', 'row', 'group', 'atomic_mass', 'electron_affinity', 'min_oxidation_state', 'max_oxidation_state', 'ionization_energy', 'nuclear_charge', 'ident', 'topology', 'size']
        if properties is None:
            properties = default_props
        if any(prop not in default_props for prop in properties):
            raise ValueError(f'Invalid property. Valid properties are: {default_props}')

        features = []
        labels = []
        for prop in properties:
            prop_vector, label = self.get_prop_vector(graph, prop)

            # Atom features is a 3d array with shape (n_atom, depth+1, n_stats)
            atom_features = np.array([self.atom_autocorrelation(graph, atom_index=i, prop_vector=prop_vector, stats=self.atom_stats) for i in graph.nodes])

            # Here we simply use multiple statistical measures of the atom features as the molecule features
            for stat in self.molecular_stats:
                # Take statistics over all atoms to make a fixed size feature vector
                stat_features = getattr(np, stat)(atom_features, axis=0)

                # Flatten the array and add to the features and labels in the correct order
                stat_labels = []
                for depth in range(stat_features.shape[0]):
                    for i, sub_stat in enumerate(self.atom_stats):
                        stat_labels.append(f'{label}-{depth}-{stat}-{sub_stat}')
                        features.append(stat_features[depth,i])
                labels.extend(stat_labels)

        features = np.array(features)

        if return_labels:
            return features, labels
        else:
            return features

    def compute_graph_descriptors(self, method):
        if method == 'molsimplify':
            df = self.compute_molsimplify_graph_descriptors()
        elif method == 'own':
            df = self.compute_own_graph_descriptors()

        return df

#     ########################### MolSimplify RACs ############################
#     This outcommented code is an old version that uses MolSimplify to compute the RACs. It is kept here for reference. However, there is one bug: When reading in the smiles strings, read_smiles() function adds more hydrogens, which is not what we want. If you really want to use it, be aware of this issue.
#     def input_to_mol3D(self, input):
#         mol = mol3D()
#         if input.endswith('.xyz'):
#             mol.readfromxyz(input)
#         else:
#             mol.read_smiles(input, ff=False)
#         return mol
#
#     def compute_molsimplify_graph_descriptors(self, molecules) -> list:
#         all_descriptors = []
#         for mol in tqdm(molecules):
#
#             mol = self.input_to_mol3D(mol)
#             desc = generate_full_complex_autocorrelations(mol, depth=self.depth, oct=False)
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
# def compute_RAC_from_graph(smiles, depth=4, return_colnames=False):
#     mol = mol3D()
#     mol.read_smiles(smiles, ff=False)
#     desc = generate_full_complex_autocorrelations(mol,depth=depth, loud=False, oct=False)
#     colnames = list(itertools.chain.from_iterable(desc['colnames']))
#     results = list(itertools.chain.from_iterable(desc['results']))
#
#     if return_colnames:
#         return results, colnames
#     else:
#         return results
    
class RDKit_2D:
    def __init__(self, smiles):
        self.original_smiles_provided = smiles
        self.smiles = []
        self.errors = []
        self.mols = []
        for sm in smiles:
            # sio = sys.stderr = StringIO()
            mol = Chem.MolFromSmiles(sm, sanitize=False)
            # error = sio.getvalue()
            error = mol is None
            if error:
                self.errors.append(error)
            else:
                self.errors.append(np.nan)
                self.mols.append(mol)
                self.smiles.append(sm)
        sys.stderr = sys.__stderr__

        self.df_err = pd.DataFrame({'isosmiles': self.original_smiles_provided, 'errors': self.errors})

    
    def compute_2Drdkit(self):
        rdkit_2d_desc = []
        calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        for mol in tqdm(self.mols):
            ds = calc.CalcDescriptors(mol)
            rdkit_2d_desc.append(ds)
            
        header = [f'rdkit_{desc}' for desc in calc.GetDescriptorNames()]
        df = pd.DataFrame(rdkit_2d_desc, columns=header, index=self.smiles)
        warn_if_nan_values(df)
        
        return df

class SOAP_3D():

    def __init__(self, ase_molecules: list, r_cut: float=10.0, n_max: int=3, l_max: int=2, crossover: bool=False, average: str='off', weighting=None):
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
            average=average,
            weighting=weighting
            )

    def get_all_species_from_molecules(self):
        self.all_atoms_list = []
        for mol in self.molecules:
            self.all_atoms_list.append(list(mol.symbols))

        species = np.unique(flatten_list(self.all_atoms_list)).tolist()

        return species

    def calculate_descriptors(self, only_from_metal: bool=False, positions=None):
        """
        Calculate the SOAP descriptors for the molecules in the dataset.
        @param only_from_metal: If True, positions will be set to the position of the metal atom in the molecule. Only used if positions is None.
        @param positions: If not None, SOAP descriptors will be calculated only from the given positions. List of (x,y,z) coordinates. If None, all atoms will be used as centers for SOAP descriptors.
        """
        if positions is None and only_from_metal:
            positions = [[atoms.index(identify_metal_in_atoms_list(atoms))] for atoms in self.all_atoms_list]

        try:
            desc = self.create_SOAP_descriptors(system=self.molecules, positions=positions)
        except ValueError:
            warnings.warn('Could not calculate SOAP descriptors with specified positions because some ligands do not have these atoms. Using all molecule atoms as SOAP centers instead.')
            desc = self.create_SOAP_descriptors(system=self.molecules, positions=None)

        return desc

    def create_SOAP_descriptors(self, system, positions):
        try:
            desc = self.soap.create(system=system, centers=positions)
        except TypeError:
            # Older version of dscribe
            desc = self.soap.create(system=system, positions=positions)

        return desc

    def get_descriptor_names(self):
        n_features = self.soap.get_number_of_features()
        n_per_element = n_features // len(self.species)
        self.names = np.array(['x'*50]*n_features)

        for el in self.species:
            loc = self.soap.get_location([el, el])
            el_names = [f'soap_{el}_{i}' for i in range(n_per_element)]
            self.names[loc] = el_names

        return self.names

if __name__ == '__main__':

    # Example graph:
    # 4 -- 1
    #      |
    # 2 -- 0
    #      |
    # 5 -- 3
    G = nx.Graph()
    G.add_nodes_from([(0, {'node_label': 'C'}),
                      (1, {'node_label': 'C'}),
                      (2, {'node_label': 'H'}),
                      (3, {'node_label': 'C'}),
                      (4, {'node_label': 'C'}),
                      (5, {'node_label': 'C'})])
    G.add_edges_from([(0, 1), (0, 2), (0, 3), (1, 4), (3, 5)])
    own_result, labels = RAC(depth=4).molecule_autocorrelation(mol=G, return_labels=True)
    own_result = list(own_result)

    molsimplify_result = list(compute_RAC_from_graph('[H][C]([C][C])([C][C])', depth=4))

    df = pd.DataFrame({'own': own_result, 'molsimplify': molsimplify_result}, index=labels)
    df['diff'] = df['own'] - df['molsimplify']

    print('Done!')
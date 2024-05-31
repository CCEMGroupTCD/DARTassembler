"""
This file contains an own implementation of RAC features (after the paper from Kulik who introduced them).
"""
import warnings
warnings.filterwarnings("ignore", message='"import openbabel" is deprecated, instead use "from openbabel import openbabel"')
from typing import Union, Tuple
import networkx as nx
import numpy as np
import pandas as pd
from pymatgen.core.periodic_table import Element as Pymatgen_Element
from tqdm import tqdm
from DARTassembler.src.ligand_extraction.utilities_graph import get_reindexed_graph, get_sorted_atoms_and_indices_from_graph

def warn_if_nan_values(df):
    nan_columns = df.columns[df.isna().any()].tolist()
    n_nan = len(nan_columns)
    if n_nan > 0:
        warnings.warn(f'{n_nan} of {len(df.columns)} calculated features have NaN values for some molecules!')

class RAC:
    """
    Class for computing the RAC descriptors of molecules.

    """
    def __init__(self, depth: int=4, molecular_stats: list[str]=None, atom_stats: list[str]=None, element_label: str='node_label'):
        """
        :param depth: The depth of the autocorrelation, i.e. how often the convolution is applied.
        :param molecular_stats: List of strings such as ['sum', 'std', 'min', 'max'] specifying which statistics to calculate over the atom autocorrelation features to get the molecular features.
        :param atom_stats: List of strings such as ['sum', 'std', 'min', 'max'] specifying which statistics to calculate over the convolution features to get the atom autocorrelation features.
        :param element_label: The label of the node attribute of the networkx graph that specifies the element of the atom.

        Simple usage example:
        ``` python
            import networkx as nx
            G = nx.Graph()
            G.add_nodes_from([(0, {'node_label': 'C'}),
                              (1, {'node_label': 'C'}),
                              (2, {'node_label': 'H'}),
                              (3, {'node_label': 'C'}),
                              (4, {'node_label': 'C'}),
                              (5, {'node_label': 'C'})])
            G.add_edges_from([(0, 1), (0, 2), (0, 3), (1, 4), (3, 5)])
            features, labels = RAC(depth=4).molecule_autocorrelation(mol=G, return_labels=True)
        ```
        """
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

    def molecule_autocorrelation(self, mol: nx.Graph, properties: list[str]=None, return_labels: bool=False) -> Union[np.array, Tuple[np.array, list]]:
        """
        Calculates the molecule features.

        :param mol: The molecule as a networkx graph.
        :param properties: List of strings. Atomic properties of the graph to use for the autocorrelation. If None, the following are used and automatically calculated from pymatgen: ['electronegativity', 'row', 'group', 'atomic_mass', 'electron_affinity', 'min_oxidation_state', 'max_oxidation_state', 'ionization_energy', 'nuclear_charge', 'ident', 'topology', 'size']. This should also support to use properties which are in the graph pre-calculated for each atom, such as atomic charges or atomic energies. In this case, provide the label name of the property as a string here.
        :param return_labels: If True, the labels of the features are also returned as second element of the tuple.

        Example usage with own calculated properties, here xtb atomic charges in the graph attribute 'charge':
        ``` python
            import networkx as nx
            G = nx.Graph()
            G.add_nodes_from([(0, {'node_label': 'C', 'charge': 0.5}),
                              (1, {'node_label': 'C', 'charge': 0.5}),
                              (2, {'node_label': 'H', 'charge': 0.1}),
                              (3, {'node_label': 'C', 'charge': -0.5}),
                              (4, {'node_label': 'C', 'charge': -0.5}),
                              (5, {'node_label': 'C', 'charge': -0.5})])
            G.add_edges_from([(0, 1), (0, 2), (0, 3), (1, 4), (3, 5)])
            features, labels = RAC(depth=4).molecule_autocorrelation(mol=G, properties=['charge'], return_labels=True)
        ```
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
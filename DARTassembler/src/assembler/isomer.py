"""
This file contains the classes and methods that are used to process the input data and generate the assembled transition metal complex isomers.
"""
# Dart Assembler imports
from DARTassembler.src.metalig.geometry import align_donor_atoms
from DARTassembler.src.metalig.mol import BaseMolecule, Ligand
from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.misc.io import load_json
from DARTassembler.src.constants import chem

# Scientific package imports
from scipy.optimize import differential_evolution, linear_sum_assignment, brute
from sklearn.cluster import MeanShift, estimate_bandwidth
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import cdist
import numpy as np

# Standard library imports
from typing import Dict, Any, List, Optional, Tuple, Union
from collections import defaultdict
import itertools

# Matplotlib imports
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

# Plotly/Dash imports
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import dash

# Data manipulation imports
from networkx.readwrite import json_graph
from ase.visualize import view
from ase import Atoms
import networkx as nx
import pathlib as pl
import pandas as pd
import ase

pio.renderers.default = 'browser'

# Todo: we need to cite the original source of the covalent radii data (I am not sure where it was obtained)
elem_cov_radii = {'H': 0.32, 'He': 0.46, 'Li': 1.33, 'Be': 1.02, 'B': 0.85, 'C': 0.75, 'N': 0.71, 'O': 0.63, 'F': 0.64, 'Ne': 0.67, 'Na': 1.55, 'Mg': 1.39, 'Al': 1.26, 'Si': 1.16, 'P': 1.11,
                  'S': 1.03, 'Cl': 0.99, 'Ar': 0.96, 'K': 1.96, 'Ca': 1.71, 'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22, 'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11, 'Ni': 1.1, 'Cu': 1.12, 'Zn': 1.18,
                  'Ga': 1.24, 'Ge': 1.21, 'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17, 'Rb': 2.1, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25,
                  'Pd': 1.2, 'Ag': 1.28, 'Cd': 1.36, 'In': 1.42, 'Sn': 1.4, 'Sb': 1.4, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31, 'Cs': 2.32, 'Ba': 1.96, 'La': 1.8, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74,
                  'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68, 'Gd': 1.69, 'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.7, 'Lu': 1.62, 'Hf': 1.52, 'Ta': 1.46, 'W': 1.37, 'Re': 1.31,
                  'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33, 'Tl': 1.44, 'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01, 'Ac': 1.86, 'Th': 1.75,
                  'Pa': 1.69, 'U': 1.7, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66, 'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65, 'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57, 'Db': 1.49,
                  'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34, 'Mt': 1.29, 'Ds': 1.28, 'Rg': 1.21, 'Cn': 1.22, 'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75, 'Ts': 1.65, 'Og': 1.57}

class Isomer(BaseMolecule):
    """
    generates a BaseMolecule object from an ASE Atoms object and
    a list of ligands
    """

    def __init__(self,
                 atoms: Atoms,
                 ligands: List[Ligand],
                 default_graph: bool = True,
                 ligand_target_vectors: List[List[float]] = None,
                 ligand_origins: List[List[float]] = None,
                 metal_centers: Union[List[List[Union[str, List[float]]]], str] = None):

        self._ani = None
        orig_atoms = atoms.copy()
        self.DART_atoms = orig_atoms
        self.ligands = ligands
        self.default_graph = default_graph
        self.ligand_target_vectors = ligand_target_vectors
        self.ligand_origins = ligand_origins
        self.metal_centers = metal_centers
        self.ligand_info = self._get_ligand_info()
        self.warning = ''

        if not len(self.ligands) == len(self.ligand_target_vectors) == len(self.ligand_origins):
            raise ValueError(f"Fatal Error: The input of ligands, target vectors, ligand origins and metal centers must have the same length."
                             f"Respective lenghts given: [ligands: {len(self.ligands)}, target_vectors: {len(self.ligand_target_vectors)}, "
                             f"ligand_origins: {len(self.ligand_origins)}.")

        assert [len(vector_set) for vector_set in self.ligand_target_vectors] == [ligand.n_eff_denticities for ligand in self.ligands], \
            ValueError("Fatal Error: The number of target vectors must match the number of ligand donor atoms.")

        super().__init__(
            atomic_props=atoms.copy(),
            graph=self._get_graph(),
            global_props=None)
        print("done")



    @classmethod
    def from_dict(cls, data: Dict[str, Any]):
        """
        Creates an Isomer object from a dictionary in the correct format.
        :param data: Dictionary containing the Isomer data
        :return: Isomer object
        """
        return cls(**data)

    @classmethod
    def from_json(cls, filepath):
        """
        Loads an Isomer object from a .json file.
        :param filepath: Path to the .json file
        :return: Isomer object
        """
        data = load_json(filepath)
        return cls.from_dict(data)

    def _get_ligand_info(self) -> Dict[str, Any]:
        """
        Extract information about the ligands in this isomer.
        - unique_names: A list of unique ligand names that appears in the same
        order as they appear in the ASE Atoms object.
        :return: Dictionary containing:
                - unique_names: List of unique ligand names.
                - geometries: List of ligand geometries.
                - donors: List of donor atoms for each ligand, represented as strings.
                - charges: List of predicted charges for each ligand.
                - stoichiometries: List of ligand stoichiometries.
        """

        return {'unique_names': [ligand.unique_name for ligand in self.ligands],
                'geometries': [ligand.geometry for ligand in self.ligands],
                'donors': ['-'.join(sorted(ligand.donor_elements)) for ligand in self.ligands],
                'charges': [ligand.charge for ligand in self.ligands],
                'stoichiometries': [lig.stoichiometry for lig in self.ligands]}

    def get_metal_symbols(self) -> List[str]:
        """
        Get the symbols of the metal centers in this isomer.
        :return: List of metal symbols
        """
        return self.atoms[self._get_metal_idc()].get_chemical_symbols()

    def _get_graph(self):
        """
        Generate a graph where (default):
            - Every metal is connected to every other metal.
            - Every metal is connected to every ligand donor atom.
            - Ligand donor atoms are not connected to each other.
        :return: networkx.Graph
        """

        graph = nx.Graph()
        ligand_tags = self.DART_atoms.get_array('ligand_idx_tags')
        vector_keys = self.DART_atoms.get_array("vector_keys")
        metal_idc = self._get_metal_idc()
        donor_idc = self._get_donor_idc()


        # --- Step 1: Add all atoms as graph nodes ---
        for i, atom in enumerate(self.DART_atoms):
            graph.add_node(i, node_label=atom.symbol)

        # --- Step 3: Connect all metal atoms to each other ---
        for i in metal_idc:
            for j in metal_idc:
                if i < j:
                    graph.add_edge(i, j)

        # --- Step 4: Connect each metal to all donor atoms ---
        for metal_idx in metal_idc:
            for donor_idx in donor_idc:
                graph.add_edge(metal_idx, donor_idx)

        # --- Step 5: Add intra-ligand bonds (all atoms with same ligand tag > 0) ---
        for ligand_idx, ligand in enumerate(self.ligands):
            # Get atom indices in self.DART_atoms that belong to this ligand
            atom_indices = np.where(ligand_tags == ligand_idx + 1)[0]  # ligand indices start at 1

            if len(atom_indices) != len(ligand.graph.nodes):
                raise ValueError(f"Mismatch: ligand {ligand_idx} has {len(ligand.graph.nodes)} nodes but {len(atom_indices)} atoms in Atoms object.")

            # Map ligand-local node indices to global atom indices
            index_map = dict(zip(ligand.graph.nodes, atom_indices))

            for u, v in ligand.graph.edges():
                graph.add_edge(index_map[u], index_map[v])

        graph = nx.relabel_nodes(graph, {n: int(n) for n in graph.nodes})

        # Add node labels for BaseMolecule
        for i, atom in enumerate(self.DART_atoms):
            graph.add_node(i, node_label=atom.symbol)

        return graph


    def _get_metal_idc(self) -> List[int]:
        """
        Get the indices of metal nodes from the ase object
        :return: List of indices of metal nodes in the graph
        """
        ligand_tags = self.DART_atoms.get_array('ligand_idx_tags')
        metal_idc = np.where(ligand_tags == 0)[0]
        return metal_idc.tolist()

    def _get_ligand_idc(self):
        # Todo this functionality has not been implemented yet.
        pass

    def _get_donor_idc(self) -> List[int]:
        """
        Get the indices of donor atoms (where donor_flag == 1).
        :return: List of indices of donor atoms
        """
        donor_flags = self.DART_atoms.get_array('donor_flag')
        donor_idc = np.where(donor_flags == 1)[0]
        return donor_idc.tolist()

    def to_dict(self):
        """
        Converts the Isomer object to a fully JSON-serializable dictionary.
        """
        # Note: I seem to have encountered an issue with non serializable types in the graph data.
        # Note: This function sanitizes the graph data to ensure it can be serialized.
        def sanitize(obj):
            if isinstance(obj, dict):
                return {k: sanitize(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [sanitize(v) for v in obj]
            elif isinstance(obj, (np.integer, np.int64)):
                return int(obj)
            elif isinstance(obj, (np.floating, np.float64)):
                return float(obj)
            else:
                return obj

        raw_graph_data = json_graph.node_link_data(self.graph)
        clean_graph_data = sanitize(raw_graph_data)

        return {
            "atomic_props": {
                "x": [float(x) for x in self.DART_atoms.positions[:, 0]],
                "y": [float(y) for y in self.DART_atoms.positions[:, 1]],
                "z": [float(z) for z in self.DART_atoms.positions[:, 2]],
                "atoms": [str(s) for s in self.DART_atoms.get_chemical_symbols()],
            },
            "graph": clean_graph_data,
            "global_props": self.global_props,
            "metal_idc": [int(i) for i in self._get_metal_idc()],
            "donor_idc": [int(i) for i in self._get_donor_idc()],
            "ligand_idc": self._get_ligand_idc(),  # make sure this returns a serializable type too
            "ligand_info": sanitize(self.ligand_info),
        }

    def view_3d_graph_interactive(self) -> None:
        """
        Interactive 3D visualization of the ASE Atoms object and bonding graph using Plotly.
        """
        atoms = self.DART_atoms
        graph = self._get_graph()
        pos = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()

        # Individual arrays
        ligand_idx_tags = atoms.get_array('ligand_idx_tags')
        donor_flags = atoms.get_array('donor_flag')

        # Color and size scheme using CPK colors
        colors = []
        sizes = []
        for symbol, ligand_tag, donor_tag in zip(symbols, ligand_idx_tags, donor_flags):
            color = self.get_cpk_colors().get(symbol, 'grey')  # fallback to grey
            colors.append(color)

            # Size logic based on tags
            if ligand_tag == 0:  # metal
                sizes.append(40)
            elif donor_tag == 1:  # donor
                sizes.append(30)
            else:
                sizes.append(20)
        # Atom coordinates
        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]

        # --- Main atom spheres (no text) ---
        node_trace = go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(size=sizes, color=colors, opacity=1.0),
            name='Atoms',
            hoverinfo='text',
            text=[f"{i}: {s}" for i, s in enumerate(symbols)]
        )

        # --- Optional labels (separate trace) ---
        label_trace = go.Scatter3d(
            x=x, y=y, z=z,
            mode='text+markers',  # Add markers to trigger legend
            text=[f"{i}: {s}" for i, s in enumerate(symbols)],
            textposition="top center",
            textfont=dict(size=12, color='black'),
            marker=dict(size=1, color='rgba(0,0,0,0)'),  # Fully transparent marker
            name='Atom Labels',
            showlegend=True,
            hoverinfo='none',
            visible='legendonly'  # Hidden by default, toggleable
        )

        edge_x_solid, edge_y_solid, edge_z_solid = [], [], []
        edge_x_dashed, edge_y_dashed, edge_z_dashed = [], [], []

        metal_idc = set(self._get_metal_idc())
        donor_idc = set(self._get_donor_idc())

        for u, v in graph.edges():
            if u in graph.nodes and v in graph.nodes:
                x0, y0, z0 = pos[u]
                x1, y1, z1 = pos[v]
                is_hashed = (
                        (u in metal_idc and v in donor_idc) or
                        (v in metal_idc and u in donor_idc) or
                        (u in metal_idc and v in metal_idc and u != v)
                )

                if is_hashed:
                    edge_x_dashed += [x0, x1, None]
                    edge_y_dashed += [y0, y1, None]
                    edge_z_dashed += [z0, z1, None]
                else:
                    edge_x_solid += [x0, x1, None]
                    edge_y_solid += [y0, y1, None]
                    edge_z_solid += [z0, z1, None]

        # Solid bonds (non donor-metal)
        edge_trace_solid = go.Scatter3d(
            x=edge_x_solid, y=edge_y_solid, z=edge_z_solid,
            mode='lines',
            line=dict(width=5, color='gray'),
            hoverinfo='none',
            name='Normal Bonds'
        )

        # Dashed donor-metal bonds
        edge_trace_dashed = go.Scatter3d(
            x=edge_x_dashed, y=edge_y_dashed, z=edge_z_dashed,
            mode='lines',
            line=dict(width=10, color='cornflowerblue', dash='solid'),
            hoverinfo='none',
            name='Donor–Metal Bonds'
        )

        # Assemble and render
        fig = go.Figure(data=[edge_trace_dashed, edge_trace_solid, node_trace, label_trace])
        fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
            ),
            margin=dict(l=0, r=0, t=0, b=0),
            showlegend=True,  # <-- Enable legend
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.2,
                xanchor="center",
                x=0.5
            ),
        )
        fig.show()

    @staticmethod
    def get_cpk_colors():
        """
        Returns a dictionary mapping element symbols to their CPK colors.
        :return:
        """
        return {
            'H': '#D3D3D3', 'D': '#FFFFC0', 'T': '#FFFFA0', 'He': '#D9FFFF', 'Li': '#CC80FF',
            'Be': '#C2FF00', 'B': '#FFB5B5', 'C': '#909090', 'C-13': '#505050', 'C-14': '#404040',
            'N': '#3050F8', 'N-15': '#105050', 'O': '#FF0D0D', 'F': '#90E050', 'Ne': '#B3E3F5',
            'Na': '#AB5CF2', 'Mg': '#8AFF00', 'Al': '#BFA6A6', 'Si': '#F0C8A0', 'P': '#FF8000',
            'S': '#FFFF30', 'Cl': '#1FF01F', 'Ar': '#80D1E3', 'K': '#8F40D4', 'Ca': '#3DFF00',
            'Sc': '#E6E6E6', 'Ti': '#BFC2C7', 'V': '#A6A6AB', 'Cr': '#8A99C7', 'Mn': '#9C7AC7',
            'Fe': '#E06633', 'Co': '#F090A0', 'Ni': '#50D050', 'Cu': '#C88033', 'Zn': '#7D80B0',
            'Ga': '#C28F8F', 'Ge': '#668F8F', 'As': '#BD80E3', 'Se': '#FFA100', 'Br': '#A62929',
            'Kr': '#5CB8D1', 'Rb': '#702EB0', 'Sr': '#00FF00', 'Y': '#94FFFF', 'Zr': '#94E0E0',
            'Nb': '#73C2C9', 'Mo': '#54B5B5', 'Tc': '#3B9E9E', 'Ru': '#248F8F', 'Rh': '#0A7D8C',
            'Pd': '#006985', 'Ag': '#C0C0C0', 'Cd': '#FFD98F', 'In': '#A67573', 'Sn': '#668080',
            'Sb': '#9E63B5', 'Te': '#D47A00', 'I': '#940094', 'Xe': '#429EB0', 'Cs': '#57178F',
            'Ba': '#00C900', 'La': '#70D4FF', 'Ce': '#FFFFC7', 'Pr': '#D9FFC7', 'Nd': '#C7FFC7',
            'Pm': '#A3FFC7', 'Sm': '#8FFFC7', 'Eu': '#61FFC7', 'Gd': '#45FFC7', 'Tb': '#30FFC7',
            'Dy': '#1FFFC7', 'Ho': '#00FF9C', 'Er': '#00E675', 'Tm': '#00D452', 'Yb': '#00BF38',
            'Lu': '#00AB24', 'Hf': '#4DC2FF', 'Ta': '#4DA6FF', 'W': '#2194D6', 'Re': '#267DAB',
            'Os': '#266696', 'Ir': '#175487', 'Pt': '#D0D0E0', 'Au': '#FFD123', 'Hg': '#B8B8D0',
            'Tl': '#A6544D', 'Pb': '#575961', 'Bi': '#9E4FB5', 'Po': '#AB5C00', 'At': '#754F45',
            'Rn': '#428296', 'Fr': '#420066', 'Ra': '#007D00', 'Ac': '#70ABFA', 'Th': '#00BAFF',
            'Pa': '#00A1FF', 'U': '#008FFF', 'Np': '#0080FF', 'Pu': '#006BFF', 'Am': '#545CF2',
            'Cm': '#785CE3', 'Bk': '#8A4FE3', 'Cf': '#A136D4', 'Es': '#B31FD4', 'Fm': '#B31FBA',
            'Md': '#B30DA6', 'No': '#BD0D87', 'Lr': '#C70066', 'Rf': '#CC0059', 'Db': '#D1004F',
            'Sg': '#D90045', 'Bh': '#E00038', 'Hs': '#E6002E', 'Mt': '#EB0026'
        }


class IsomerFactory:
    def __init__(self, ligands: Dict[int, Ligand], target_vectors: List[Dict[Any, List[float]]], ligand_origins: List[List[float]], metal_origins: List[List[float]],
                 metal_types: List[str], filter_duplicate_isomers: bool = True, filter_clashing_structures: bool = True,
                 filter_clashing_structures_cov_radii_buffer: float = -0.3,
                 check_metal_clashes: bool = True,
                 filter_duplicate_isomers_method: str = "fingerprint",
                 filter_duplicate_isomers_grid_size: int = 9,
                 isomer_comparison_mode: str = "max_diff",
                 isomer_comparison_grouping_mode: str = "cluster",
                 isomer_comparison_grouping_cutoff: float = 1.0,
                 swap_groups: Optional[List[int]] = None):
        """
        Generates novel transition metal complexes from ligands and metals
        :param ligands:         Dictionary of ligand objects
        :param target_vectors:  List of dictionaries containing target vectors for each ligand
        :param ligand_origins:  List of ligand origins
        :param metal_origins:   List of metal origins
        :param metal_types:     List of metal types

        """
        # Define the class variables
        self.ligands = ligands  # for each ligand, get the ligand object
        self.target_vectors = target_vectors  # List of target vectors i.e [[0, 0, 0], [0, 0, 1], ...]
        self.ligand_origins = ligand_origins  # List of ligand_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_origins = metal_origins  # List of metal_origins i.e [[0, 0, 0], [0, 0, 1], ...]
        self.metal_types = metal_types  # List of metal types i.e ['Ru', 'Mn', ...]
        self.filter_duplicate_isomers = filter_duplicate_isomers  # Boolean to determine if duplicate isomers should be filtered
        self.filter_clashing_structures = filter_clashing_structures  # Boolean to determine if clashing structures should be filtered
        self.filter_clashing_structures_cov_radii_buffer = filter_clashing_structures_cov_radii_buffer  # Buffer for clashing structures filtering
        self.check_metal_clashes = check_metal_clashes  # Boolean to determine if metal clashes should be checked
        self.filter_duplicate_isomers_method = filter_duplicate_isomers_method  # Method for filtering duplicate isomers
        self.filter_duplicate_isomers_grid_size = filter_duplicate_isomers_grid_size
        self.isomer_comparison_mode = isomer_comparison_mode
        self.isomer_comparison_grouping_mode = isomer_comparison_grouping_mode
        self.isomer_comparison_grouping_cutoff = isomer_comparison_grouping_cutoff
        self.swap_groups = swap_groups

        # Validate the input
        self._validate_input()
        self._validate_swap_groups()



    def _validate_input(self) -> None:
        """
        Validates the input to the IsomerFactory class
        :raises: ValueError: If the input is invalid
        :return: None
        """
        if len(self.ligands) != len(self.target_vectors) or len(self.ligands) != len(self.ligand_origins):
            raise ValueError(
                f"Fatal Error: Ligand objects [{len(self.ligands)}], target vectors [{len(self.target_vectors)}], and ligand origins [{len(self.ligand_origins)}] must have the same length")

        if len(self.metal_origins) != len(self.metal_types):
            raise ValueError(
                f"Fatal Error: Metal origins [length: {self.metal_origins}] and metal types [length: {self.metal_types}] must have the same length for multi-metallic systems.", )

    def _validate_swap_groups(self) -> None:
        """
        Validates the swap groups for ligands. Ligands of the same 'swap_group' must have the same effective ligand
        coordination number (n_eff_denticities) to be considered for swapping. Will raise an error if this is not the case.
        :raises: ValueError: If the swap groups are not valid
        :return: None
        """
        if any(swap_tag is None for swap_tag in self.swap_groups) and not all(swap_tag is None for swap_tag in self.swap_groups):
            raise ValueError("Fatal Error: If a swap_group is specified for any ligand, it MUST be specified for all ligands.")

        elif all(swap_tag is None for swap_tag in self.swap_groups):
            print("No swap groups specified. Using effective ligand coordination numbers (n_eff_denticities) as swap groups.")
            self.swap_groups = [ligand.n_eff_denticities for ligand in self.ligands]

        else:
            # All swap_groups are present; now check if ligands that are specified by the user to be swapped have the
            # same effective ligand coordination number (n_eff_denticities) and make sense to be swapped
            list1 = self.swap_groups
            list2 = [ligand.n_eff_denticities for ligand in self.ligands]
            groups = defaultdict(list)

            for i, v in enumerate(list1):
                groups[v].append(i)

            for group_id, indices in groups.items():
                elcns = {list2[i] for i in indices}
                if len(elcns) > 1:
                    details = ", ".join(f"{self.ligands[i].unique_name} (n_eff_denticities={list2[i]})" for i in indices)
                    raise ValueError(f"Fatal Error: Ligands in swap_group {group_id} must all have the same elcn, but found mismatches: {details}")
        print("Swap groups validated successfully.")
        return None

    def generate(self) -> List[Atoms]:
        # Generate all possible isomers by exchanging compatible ligands
        ligand_lists = self._assign_ligands_to_vectors(
            ligands=list(self.ligands),
            vectors=self.target_vectors,
            swap_groups=self.swap_groups
        )

        # Loop through each ligand set permutation and generate each isomer
        all_isomers = []
        all_ligands = []  # List to store the DARTLigand objects used in each isomer
        for ligand_swapped_combo in ligand_lists:
            aligned_ligands = []  # list of a complexes ligands i.e. [l1, l2, [l3-i1, l3-i2], l4, ...] where i1, i2, etc. correspond to a unique coordination of a ligand within a binding site
            aligned_DARTLigands = []  # list of DARTLigand objects used in the complex
            for idx, (ligand, ligand_target_vectors, origin) in enumerate(zip(ligand_swapped_combo, self.target_vectors, self.ligand_origins)):

                # Retrieve the ligand's geometry (ASE atoms object) and donor atom indices (List[List[index]).
                geometry, donor_atoms = ligand.get_isomers_effective_ligand_atoms_with_effective_donor_indices()

                n_atoms = len(geometry)
                donor_idc = list(self._flatten(ligand._get_denticities_and_hapticities_idc()))
                donor_flags = np.zeros(n_atoms, dtype=int)
                donor_flags[donor_idc] = 1

                geometry.new_array('n_denticity_tags', np.full(n_atoms, ligand.n_eff_denticities, dtype=int))
                geometry.new_array('ligand_idx_tags', np.full(n_atoms, idx + 1, dtype=int))
                geometry.new_array('donor_flag', donor_flags)

                # Convert target vector dictionary values to numpy arrays.
                vector_labels = list(ligand_target_vectors.keys())
                target_vectors = [np.array(v) for v in ligand_target_vectors.values()]

                # Align the donor atoms of the ligand to the target vectors
                ligand_isomers, donor_atoms_ordered, rssd = self.try_all_geometrical_isomer_possibilities(atoms=geometry,
                                                                                                          donor_idc=donor_atoms[0],
                                                                                                          vector_labels=vector_labels,
                                                                                                          target_vectors=target_vectors)

                # Translate the ligand to its correct location in the complex
                for ligand_isomer in ligand_isomers:
                    # Note: This method assumes that the ligand has been pre-translated to 0,0,0
                    ligand_isomer.set_positions(ligand_isomer.get_positions() + np.array(origin))

                # Remove dummy atoms (e.g., "Cu") from haptic ligands.
                cleaned_isomers = self._remove_haptic_dummy_atom(
                    atoms_list=ligand_isomers,
                    dummy_atom="Cu",
                    donor_atoms_idc=ligand.hapdent_idc
                )
                # Append the rotated ligands to the list
                aligned_ligands.append(cleaned_isomers)
                aligned_DARTLigands.append([ligand for _ in range(len(cleaned_isomers))])

            all_isomers.append(aligned_ligands)
            all_ligands.append(aligned_DARTLigands)

        # Assemble the complexes using the ligand permutations
        all_isomers, ligands_used = self._gen_all_isomers(ligands=all_isomers, DARTLigands=all_ligands)

        dart_isomers = [Isomer(atoms=isomer,
                               ligands=ligands,
                               default_graph=True,
                               ligand_target_vectors=self.target_vectors,
                               ligand_origins=self.ligand_origins,
                               metal_centers=self.metal_origins
                               ) for isomer, ligands in zip(all_isomers, ligands_used)]

        return dart_isomers

    @staticmethod
    def try_all_geometrical_isomer_possibilities(
            atoms: ase.Atoms,
            donor_idc: list[int],
            target_vectors: list[np.ndarray],
            vector_labels: list[str] = None
    ) -> tuple[list[ase.Atoms], list[list[int]], float]:
        """
        Tries out all combinations of mapping the provided donor atoms to the target vectors, for each target vector in the list of target vectors. All combinations are tried out and all isomers with the lowest rssd are returned. Usually, there will be several with the same rssd.
        :param: atoms: ase.Atoms() object of the ligand
        :param: vector_labels: List of labels for the target vectors. If provided, these labels will be assigned to the donor atoms in the isomers.
        :param: donor_idc: Indices of the donor atoms in `atoms`
        :param: target_vectors: List of 3D vectors of shape (n,3)
        :return: Tuple of:
            - List of ASE Atoms objects of the best isomers
            - List of lists of indices of the best isomers
            - The root sum of squared differences (RSSD) of the best isomers
        """
        target_vectors = np.array(target_vectors)
        if target_vectors.ndim == 1:
            target_vectors = target_vectors[None, :]
        assert target_vectors.shape[-1] == 3, f'Wrong dimension of target_vectors. It\'s {target_vectors.shape[-1]} but must be 3.'

        data = []
        n = len(target_vectors)
        donor_idc_permutations = list(itertools.permutations(donor_idc, n))
        for idc in donor_idc_permutations:
            idc = list(idc)  # Convert tuple to list to allow indexing later on with the indices saved in `data`
            isomer, rssd = align_donor_atoms(atoms, donor_idc=idc, target_vectors=target_vectors, return_rssd=True)

            # Assign vector_key to each donor atom in the permutation
            # This allows us to later attribute each donor atom with a specific target vector
            if vector_labels:
                vector_keys_array = np.full(len(isomer), "", dtype=object)
                for atom_idx, label in zip(idc, vector_labels):
                    vector_keys_array[atom_idx] = label
                isomer.new_array("vector_keys", vector_keys_array)

            data.append((rssd, idc, isomer))
        best_rssd, _, _ = min(data, key=lambda x: x[0])

        all_best_rssd = []
        all_best_idc = []
        all_best_isomers = []
        for rssd, vectors, isomer in data:
            if np.isclose(rssd, best_rssd):
                all_best_rssd.append(rssd)
                all_best_idc.append(vectors)
                all_best_isomers.append(isomer)
        best_rssd = np.mean(all_best_rssd)

        return all_best_isomers, all_best_idc, best_rssd

    def _flatten(self, nested):
        """
        Recursively flattens a nested list or tuple into a flat generator of elements.
        :param: nested (list or tuple): Arbitrarily nested list/tuple of elements (e.g., indices)
        :return: Individual flattened elements
        """
        for item in nested:
            if isinstance(item, (list, tuple)):
                yield from self._flatten(item)
            else:
                yield item

    @staticmethod
    def _assign_ligands_to_vectors(ligands: List[Ligand],
                                   vectors: List[Dict[str, List[float]]],
                                   swap_groups: List[int]) -> List[List[Ligand]]:
        """
        Assigns ligands to vector entries based on swap group IDs.
        Ligands in the same swap group can be swapped among vector entries assigned the same swap group ID.

        :param ligands:       A list of ligand objects.
        :param vectors:       A list of dictionaries, each representing a coordination site with arbitrary keys.
        :param swap_groups:   A list of integers where the i-th value defines the swap group for the i-th ligand.
                              Ligands with the same group ID can be swapped; those with different IDs cannot.
                              This list must match the length and order of the `vectors` list.
        :return:              A list of ligand groupings for each vector entry.
                              Each entry contains the swappable ligands for that vector,
                              used later for permutation generation.
        :raise:               LoggedValueError if swap_groups are inconsistent with ligand or vector count.
        """

        # Validate swap_groups length
        if len(vectors) != len(swap_groups):
            raise ValueError("Fatal Error: swap_groups must match the number of target vectors.")
        if len(ligands) != len(swap_groups):
            raise ValueError("Fatal Error: swap_groups must match the number of ligands.")

        # Step 1: Build a mapping of swap_group ID → list of ligands that belong to that group.
        group_to_ligands = defaultdict(list)
        for ligand, group_id in zip(ligands, swap_groups):
            group_to_ligands[group_id].append(ligand)

        # Step 2: Build the vector-to-ligand group mapping.
        # For each vector entry, assign the appropriate ligand group based on the vector’s corresponding swap group ID.
        _list = [group_to_ligands[group_id] for group_id in swap_groups]

        # Step 3: Generate all valid permutations of ligand assignments.
        # For each vector position, choose one ligand from the allowed group — ensuring no ligand is reused.
        results = []
        current_permutation = []
        used = set()

        def backtrack(index: int):
            # Base case: if all vector positions have been filled, record the permutation.
            if index == len(_list):
                results.append(current_permutation.copy())
                return

            # Iterate over allowed ligands for the current position.
            for ligand in _list[index]:
                if ligand in used:
                    continue  # Skip ligands already assigned in this permutation.
                used.add(ligand)
                current_permutation.append(ligand)
                backtrack(index + 1)
                # Backtrack: remove the ligand from current permutation and reuse it elsewhere.
                current_permutation.pop()
                used.remove(ligand)

        backtrack(0)
        return results

    def _add_metals(self, ligand_structure: Atoms) -> Atoms:
        """
        Adds the metals to the complex, ensuring they are first in the Atoms object.
        :param ligand_structure: Atoms object representing the ligand structure
        :return: Atoms object with metals prepended
        """
        # Start with an empty Atoms object to accumulate metals
        metal_atoms = Atoms()

        for idx, (metal_type, metal_origin) in enumerate(zip(self.metal_types, self.metal_origins)):
            # Create a single-atom ASE structure for the metal
            metal = Atoms(symbols=metal_type, positions=[metal_origin])
            for tag, value in {'n_denticity_tags': 0, 'ligand_idx_tags': 0, 'donor_flag': 0}.items():
                metal.new_array(tag, np.array([value], dtype=int))
            metal_atoms = self.merge_atoms_with_arrays(metal_atoms, metal)
            print(f"Metal [{idx + 1}/{len(self.metal_types)}]: [{metal_type}] added at origin [{metal_origin}]")

        # Combine metals first, then ligand structure
        combined = metal_atoms + ligand_structure
        return combined

    def _gen_all_isomers(
            self,
            ligands: List[List[List[Atoms]]],
            DARTLigands: List[List[List[Ligand]]]
    ) -> Tuple[List[Atoms], List[List[Ligand]]]:
        """
        Generate all possible isomers from a list of ligands which have multiple isomers.

        :param ligands: Nested list of ASE ligand structures.
        :param DARTLigands: Nested list of corresponding Ligand objects.
        :return: Tuple of (assembled ASE Atoms objects, per-isomer list of Ligand objects used).
        """
        isomers = []
        ligand_compositions = []

        print(f"Generating isomers from {len(ligands)} ligand slots.")

        for ligand_atoms_list, ligand_obj_list in zip(ligands, DARTLigands):
            combinations_atoms = list(itertools.product(*ligand_atoms_list))
            combinations_ligands = list(itertools.product(*ligand_obj_list))

            for idx, (combo_atoms, combo_ligands) in enumerate(zip(combinations_atoms, combinations_ligands)):
                combined = combo_atoms[0].copy()
                for ligand_atoms in combo_atoms[1:]:
                    combined = self.merge_atoms_with_arrays(combined, ligand_atoms)

                combined = self._add_metals(ligand_structure=combined)
                isomers.append(combined)
                ligand_compositions.append(list(combo_ligands))

                print(f"Isomer {idx + 1}/{len(combinations_atoms)} assembled.")

        print(f"Total number of assembled isomers: {len(isomers)}.")
        return isomers, ligand_compositions

    @staticmethod
    def merge_atoms_with_arrays(base: Atoms, to_add: Atoms) -> Atoms:
        """
        Merges two ASE Atoms objects, combining their positions, numbers, and custom arrays.
        :param base: The base Atoms object to merge into.
        :param to_add: The Atoms object to add to the base.
        :return: A new Atoms object containing the merged data.
        """
        combined = base.copy()
        combined += to_add

        # Merge custom arrays (ignoring positions/numbers etc. already handled by +=)
        for name in set(base.arrays.keys()).union(to_add.arrays.keys()):
            if name in {'positions', 'numbers'}:
                continue  # Already handled
            a1 = base.arrays.get(name)
            a2 = to_add.arrays.get(name)

            if a1 is not None and a2 is not None:
                merged = np.concatenate([a1, a2], axis=0)
            elif a1 is not None:  # a2 missing
                pad_shape = (len(to_add),) + a1.shape[1:]
                merged = np.concatenate([a1, np.zeros(pad_shape, dtype=a1.dtype)], axis=0)
            elif a2 is not None:  # a1 missing
                pad_shape = (len(base),) + a2.shape[1:]
                merged = np.concatenate([np.zeros(pad_shape, dtype=a2.dtype), a2], axis=0)
            else:
                continue  # neither exists, shouldn't happen
            combined.set_array(name, merged)

        return combined

    @staticmethod
    def _remove_haptic_dummy_atom(atoms_list: List[Atoms], dummy_atom: str, donor_atoms_idc: Tuple[Tuple[int]]):
        """
        Removes the dummy atom from the generated isomers
        :return: List[Atoms]
        """
        # Check to see if there is haptic coordination
        haptic_coordination = False
        for donor_atoms in donor_atoms_idc:
            if type(donor_atoms) == tuple:
                haptic_coordination = True
                break
            else:
                pass

        # If there is no haptic coordination, return the atoms list as is
        if not haptic_coordination:
            return atoms_list

        # If there is haptic coordination, remove the dummy atom from the donor atoms
        else:
            for atoms in atoms_list:
                dummy_idc = [i for i, atom in enumerate(atoms) if atom.symbol == dummy_atom]
                dummy_idc.sort(reverse=True)  # This is important so that the larger index is removed first so as not to change the index of the other atoms
                for dummy_idx in dummy_idc:
                    atoms.pop(dummy_idx)
            return atoms_list


class AxialOptModifier:
    def __init__(self, isomers: List[Isomer], opt: bool = True):
        """
        This class will take a list of ASE Atoms objects and optimize mono-coordinating ligands around their coordination axis
        """
        self.input_isomers = isomers
        self.opt_command = opt
        self.output_isomers = []
        print(f"AxialOpt initialized with {len(self.input_isomers)} Isomer objects.")



    def modify(self):
        """
        Optimize the rotation of each mono-coordinating ligand around their respective coordination axis simultaneously
        """
        # Set a random seed for reproducibility
        np.random.seed(42)

        # If opt_command is False, return the input isomers without optimization
        if not self.opt_command:
            return self.input_isomers

        # Loop through each of the inputted complexes
        for isomer in self.input_isomers:

            atoms = isomer.DART_atoms.copy()
            target_vectors = isomer.ligand_target_vectors
            ligand_origins = isomer.ligand_origins

            # Run the global optimizer.
            bounds = [[0, 360] for _ in target_vectors]
            result = differential_evolution(
                self.objective_function, bounds=bounds,
                args=(target_vectors, ligand_origins, atoms)
            )

            # Retrieve the relevant array
            ligand_idx_tags = atoms.get_array('ligand_idx_tags')
            n_denticity_tags = atoms.get_array('n_denticity_tags')

            unique_tags_indices_set = np.unique(ligand_idx_tags)
            unique_tags_indices_set = unique_tags_indices_set[unique_tags_indices_set != 0]

            # Loop through the unique second tags
            for angle, axis, origin, tag in zip(list(result.x), target_vectors, ligand_origins, unique_tags_indices_set):
                # Get indices where the second tag is the current tag
                indices = np.where(ligand_idx_tags == tag)[0]

                # Check if any of these indices have a first tag not equal to 1, if so, skip
                if np.any(n_denticity_tags[indices] != 1):
                    continue

                atoms = self.rotate(atoms=atoms, vector=np.array(list(axis.values())[0]), origin=np.array(origin), idc=indices, angle=angle).copy()
                print(f"Rotated ligand with tag {tag} by {angle:.2f} degrees around vector {axis} at origin {origin}.")

            # Create new Isomer from rotated structure
            new_dart = Isomer(
                atoms=atoms.copy(),
                ligands=isomer.ligands,
                ligand_target_vectors=isomer.ligand_target_vectors,
                ligand_origins=isomer.ligand_origins,
                metal_centers=isomer.metal_centers,
                default_graph=isomer.default_graph
            )

            # Append the new Isomer to the output complexes
            self.output_isomers.append(new_dart)

        print(f"Optimized {len(self.output_isomers)} complexes with mono-coordinating ligand rotations.")
        return self.output_isomers

    def visualize_structures(self):
        """
        Visualize input and output complexes interleaved using ASE's GUI.
        Each input is followed by its corresponding optimized output.
        """
        if not self.output_isomers:
            print("No output complexes found. Run opt_mono_rotation() first.")
            return

        structures_to_view = []

        for i, (in_isomer, out_isomer) in enumerate(zip(self.input_isomers, self.output_isomers)):
            input_copy = in_isomer.DART_atoms.copy()
            input_copy.info["label"] = f"Input {i}"

            output_copy = out_isomer.DART_atoms.copy()
            output_copy.info["label"] = f"Optimized {i}"

            structures_to_view.extend([input_copy, output_copy])  # interleave input/output

        print(f"Launching viewer for {len(structures_to_view)} structures...")
        view(structures_to_view)

    def objective_function(self, x: np.ndarray, vectors_in: List[np.array], origins_in: List[np.array], TMC_in: Atoms, ):
        """
        Objective function to optimize the position of the ligands in the TMC complex.
        :param: x:  Array of angles to rotate each ligand around its respective vector.
        :param: vectors_in: List of vectors for each ligand, where each vector is a dictionary with a single key-value pair.
        :param: origins_in: List of origins for each ligand, where each origin is a list of three floats.
        :param: TMC_in: ASE Atoms object representing the transition metal complex (TMC) to be optimized.
        :return: float: The penalty score based on interatomic distances.
        """
        # Generate a copy of the input complex
        TMC_worker = TMC_in.copy()

        # Retrieve the relevant array
        ligand_tags = TMC_worker.get_array("ligand_idx_tags")
        denticity_tags = TMC_worker.get_array("n_denticity_tags")

        # Get the unique set of all nonzero second tags (each unique tag represents a unique ligand, a zero tag represents a metal).
        unique_tags = [tag for tag in np.unique(ligand_tags) if tag != 0]

        # Loop through each unique tag and apply the rotation
        for tag, angle, axis, origin in zip(unique_tags, list(x), vectors_in, origins_in):

            # Get indices where the second tag in the list is equal to the current "tag" (essentially the indices of the atoms in this particular ligand)
            indices = np.where(ligand_tags == tag)[0]

            # Check if any of these indices have a first tag not equal to 1, if so, skip
            # This ensures only ligands that have a effective coordination number of 1 are rotated
            if np.any(denticity_tags[indices] != 1):
                continue  # Skip this ligand group

            TMC_worker = self.rotate(atoms=TMC_worker, vector=np.array(list(axis.values())[0]), origin=np.array(origin), idc=indices, angle=angle).copy()

        # Get the interatomic distance matrix
        distance_matrix = TMC_worker.get_all_distances()

        # Set the diagonal to a large number to avoid self-interaction (or np.inf)
        np.fill_diagonal(distance_matrix, np.inf)

        # Calculate the penalty: for each pair with d <= 4, add 1/d^2 to the penalty
        penalty = np.sum(1.0 / (distance_matrix ** 2))

        return penalty

    @staticmethod
    def rotate(atoms: Atoms, vector: np.array, origin: np.array, idc: List[int], angle: int):
        """
        Rotate the atoms in the Atoms object (only atoms with indices=idc) around the vector by the specified angle.
        :param atoms: Atoms object to rotate.
        :param vector: vector to rotate around.
        :param origin: origin of the rotation.
        :param idc: indices of the atoms to rotate.
        :param angle: the angle to rotate the atoms by in degrees.
        :return: an ase.Atoms object with the rotated atoms.
        """

        # Normalize rotation vector
        vector = np.asarray(vector, dtype=float)
        vector /= np.linalg.norm(vector)

        # Create rotation object
        rotation = R.from_rotvec(np.radians(angle) * vector)

        # Copy the atoms object to avoid modifying the original
        rotated_atoms = atoms.copy()

        # Apply rotation to selected atoms
        for i in idc:
            pos = atoms.positions[i] - origin  # Translate to origin
            rotated_pos = rotation.apply(pos) + origin  # Rotate and translate back
            rotated_atoms.positions[i] = rotated_pos

        return rotated_atoms


class DuplicateIsomerFilter:
    """
    Class to reduce the number of isomers based on alignment or fingerprint similarity.
    """

    def __init__(self, isomers: List[Isomer],
                 method: str = "fingerprint",
                 grid_size=9,
                 isomer_comparison_mode: str = "max_diff",
                 isomer_comparison_grouping_mode: str = "cluster",  # 'cluster' or 'cutoff'
                 fingerprint_grouping_cutoff: float = 1.0,
                 metal_centres: List[List[float]] = None,
                 energy_heuristic_mode: str = "max",
                 ):
        """
        Initialize the isomer reduction class.
        :param: isomers: The list of ASE Atoms objects representing isomers.
        :param: method: The method to use for reduction, either 'alignment' or 'fingerprint'.
        :param: grid_size: The number of grid points when scanning from 0 to 360 for exact alignment.
        :param: isomer_comparison_mode: The mode for comparing fingerprints, e.g., 'max_diff', etc.
        :param: isomer_comparison_grouping_mode: The mode for grouping fingerprints, either 'cluster' or 'cutoff'.
        :param: fingerprint_grouping_cutoff: The cutoff value for grouping fingerprints when using 'cutoff' mode.
        """
        self.isomers = isomers
        self.method = method.lower()
        self.grid_size = grid_size  # The number of grid points when scanning from 0 to 360 for exact alignment
        self.isomer_comparison_mode = isomer_comparison_mode
        self.isomer_comparison_grouping_mode = isomer_comparison_grouping_mode
        self.fingerprint_grouping_cutoff = fingerprint_grouping_cutoff
        self.metal_centres = metal_centres
        self.diff_matrix = None  # Placeholder for the fingerprint difference matrix
        self.energy_heuristic_mode = energy_heuristic_mode
        self.output_isomers = []
        self.similarity_cutoff_used = None



    def filter(self) -> List[Atoms]:
        """
        Reduce the number of isomers based on the specified method.
        :return: unique isomers as a list of ASE Atoms objects
        """

        if self.method == "alignment":
            self.output_isomers = self._reduce_by_alignment()
        elif self.method == "fingerprint":
            self.output_isomers = self._reduce_by_fingerprint()
        else:
            raise ValueError(f"Fatal Error: Unsupported reduction method '{self.method}'. "
                             f"Supported methods are 'alignment' and 'fingerprint'.")
        print(f"Reduced isomers from {len(self.isomers)} to {len(self.output_isomers)} using method '{self.method}'.")
        return self.output_isomers

    def _reduce_by_alignment(self) -> List[Atoms]:
        """
        Reduce isomers by aligning them and calculating RMSD or another distance metric.
        :return:
        """

        n = len(self.isomers)
        diff_matrix = np.zeros((n, n))
        counter = 0
        total = n * (n - 1) // 2

        for i in range(n):
            isomer1 = self.isomers[i]
            for j in range(i + 1, n):
                isomer2 = self.isomers[j]
                print(f"Comparing isomers [{i}] and [{j}] ... [{counter}/{total}]")
                score = self.align_isomers(isomer1.DART_atoms, isomer2.DART_atoms)
                diff_matrix[i, j] = score
                diff_matrix[j, i] = score
                counter += 1

        self.diff_matrix = diff_matrix

        group_labels_matrix = self._analyze_similarity(
            self.diff_matrix,
            quantile=0.2,
            method=self.isomer_comparison_grouping_mode,
            cutoff=self.fingerprint_grouping_cutoff
        )

        # Assign 'duplicate' warning to all but the first occurrence in each group
        visited = set()
        for i in range(n):
            if i in visited:
                continue
            for j in range(n):
                if i != j and group_labels_matrix[i, j] == "Close":
                    self.isomers[j].warning += "duplicate-"
                    visited.add(j)

        self.output_isomers = self.isomers  # return all isomers, tagged if duplicate
        return self.output_isomers

    def energy_heuristic(self, stat_atoms: ase.Atoms, rot_atoms: ase.Atoms):
        """
        Heuristic to determine the energy of the isomer to be minimized.

        For each element type, pair each atom in stat_atoms with an atom in rot_atoms
        (both having the same number of atoms for that element) such that:
          - If mode="sum": the sum of distances is minimized.
          - If mode="max": the maximum paired distance is returned.

        :param: stat_atoms: ase.Atoms object (stationary isomer)
        :param: rot_atoms: ase.Atoms object (rotated isomer)
        :return: The energy metric (either sum or max of paired distances).
        """
        # Initialize the total energy and overall max distance.
        total_energy = 0.0
        overall_max_distance = 0.0

        # Get chemical symbols from both isomers
        stat_symbols = stat_atoms.get_chemical_symbols()
        rot_symbols = rot_atoms.get_chemical_symbols()

        # Precompute indices for each element.
        stat_index_dict = {element: [i for i, sym in enumerate(stat_symbols) if sym == element]
                           for element in set(stat_symbols)}
        rot_index_dict = {element: [i for i, sym in enumerate(rot_symbols) if sym == element]
                          for element in set(rot_symbols)}

        # Loop over each unique element in the stationary isomer.
        for element, stat_indices in stat_index_dict.items():
            if element not in rot_index_dict:
                raise ValueError(f"Element {element} missing in rotated isomer")
            rot_indices = rot_index_dict[element]
            if len(stat_indices) != len(rot_indices):
                raise ValueError(f"Mismatch in number of '{element}' atoms between isomers: "
                                 f"{len(stat_indices)} vs. {len(rot_indices)}")
            # Get positions for the current element.
            stat_positions = stat_atoms.positions[stat_indices]
            rot_positions = rot_atoms.positions[rot_indices]

            # Compute the distance matrix using C-optimized function.
            distance_matrix = cdist(stat_positions, rot_positions)

            # Solve the assignment problem to pair atoms optimally
            row_ind, col_ind = linear_sum_assignment(distance_matrix)

            if self.energy_heuristic_mode == "sum":
                total_energy += distance_matrix[row_ind, col_ind].sum()
            elif self.energy_heuristic_mode == "max":
                element_max = distance_matrix[row_ind, col_ind].max()
                overall_max_distance = max(overall_max_distance, element_max)
            else:
                raise ValueError("Invalid mode. Use 'sum' or 'max'.")

        return total_energy if self.energy_heuristic_mode == "sum" else overall_max_distance

    def align_isomers(self, stationary_atoms: ase.Atoms, rotated_atoms: ase.Atoms) -> float:
        """
        Using Scipy's global optimizer to find the optimal rotation that minimizes the distance
        between like atoms in the two isomers. When the two isomers are optimally aligned,
        the RMSD or another measure of inter-atomic distance can be calculated.

        This method, with a large enough grid density, is considered to be the ground truth
        for the optimal rotation of the two isomers. Other Scipy global optimizers can be
        easily integrated here for testing purposes.

        :param stationary_atoms: ASE Atoms object to remain fixed.
        :param rotated_atoms: ASE Atoms object to be rotated.
        :return: float — alignment score (as defined by energy_heuristic or objective_function)
        """
        assert hasattr(self, "metal_centres"), ValueError("Fatal Error: metal_centres must be defined before calling align_isomers. ")
        print(f"Aligning isomers based on {len(self.metal_centres)} metal centre(s).")

        # Here the number of metal centres and the fact that metal centres of different isomers
        # must always be aligned is taken into account to determine if the isomers are similar or not
        if len(self.metal_centres) == 1:
            # There are 3 axes which an isomer can be rotated around to align it with another isomer
            bounds = [[0, 360] for _ in range(3)]
            # The 3 cardinal axes, properly centered on the metal centre
            axes = np.eye(3)  # Standard Cartesian axes (x, y, z)
            print("Performing 3D brute-force alignment over x, y, z axes.")

        elif len(self.metal_centres) == 2:
            # The isomer can only be rotated around the metal-metal axis to determine if the isomers are similar or not
            bounds = [[0, 360] for _ in range(1)]
            axis_vector = np.array(self.metal_centres[1]) - np.array(self.metal_centres[0])
            axis_vector /= np.linalg.norm(axis_vector)  # Normalize
            axes = [axis_vector]
            print(f"Performing 1D brute-force alignment around axis: {axis_vector.tolist()}")

        elif len(self.metal_centres) >= 3:
            # 3 metal centres means each isomer is fixed in space and their geometries can be directly compared
            # We simply return the energy heuristic
            print("Three or more metal centres detected — skipping alignment and using direct heuristic.")
            return self.energy_heuristic(stat_atoms=stationary_atoms, rot_atoms=rotated_atoms)

        else:
            raise ValueError(f"Fatal Error: Unsupported number of metal centres ({len(self.metal_centres)}). ")

        # Perform brute-force global optimization using the configured grid density
        result_angles = brute(
            func=self.objective_function,
            ranges=bounds,
            args=(stationary_atoms, rotated_atoms, axes),
            Ns=self.grid_size
        )

        min_score = self.objective_function(np.array(result_angles), stationary_atoms, rotated_atoms, axes)
        print(f"Best alignment score after brute-force search: {min_score:.4f}")
        return min_score

    def objective_function(self, x: np.ndarray, atoms1: ase.Atoms, atoms2: ase.Atoms, axes: np.array):
        """
        Objective function to optimize the position of the ligands in the TMC complex.
        :param x:
        :param atoms1:
        :param atoms2:
        :param axes:
        :return:
        """
        # Copy the input atoms to avoid modifying the original
        stationary_isomer = atoms1.copy()
        rotated_isomer = atoms2.copy()

        # Precompute the combined rotation matrix.
        R_total = self.combined_rotation_matrix(x, axes)

        # Apply the rotation to the rotated isomer
        self.apply_combined_rotation(atoms=rotated_isomer,
                                     R_total=R_total,
                                     center=np.array(self.metal_centres[0]))

        # calculate the energy heuristic that will be minimized
        val = self.energy_heuristic(stat_atoms=stationary_isomer,
                                    rot_atoms=rotated_isomer)
        return val

    @staticmethod
    def combined_rotation_matrix(angles, axes):
        """
        Compute the combined rotation matrix from a list of angles and corresponding axes.
        :param: angles: List of angles in degrees.
        :param: axes: List of rotation axes.
        """
        # Start with identity matrix.
        R_total = np.eye(3)
        for angle, axis in zip(angles, axes):
            # Create a rotation for this angle and axis.
            r = R.from_rotvec(np.deg2rad(angle) * np.array(axis))
            # Combine rotations. Note that the order matters!
            R_total = r.as_matrix() @ R_total
        return R_total

    @staticmethod
    def apply_combined_rotation(atoms, R_total, center):
        """
        Apply the combined rotation to the positions of an ASE Atoms object.
        :param: atoms: ASE Atoms object to rotate.
        :param: R_total: Combined rotation matrix.
        :param: center: Center of rotation
        """
        # Shift positions relative to center, apply rotation, then shift back.
        shifted = atoms.positions - center
        atoms.positions = center + (shifted @ R_total.T)

    def _reduce_by_fingerprint(self) -> List[Atoms]:
        """
        Reduce isomers using fingerprint-based similarity clustering.
        Select one representative per cluster labeled 'Close'.
        :return: List of unique isomers after fingerprint-based reduction.
        """
        self.diff_matrix = self._compute_fingerprint_matrix(self.isomers)

        # Cluster labels: "Close" or "Far" for each pair
        group_labels_matrix = self._analyze_similarity(self.diff_matrix, quantile=0.2,
                                                       method=self.isomer_comparison_grouping_mode,
                                                       cutoff=self.fingerprint_grouping_cutoff)

        n = len(self.isomers)
        visited = set()

        for i in range(n):
            if i in visited:
                continue
            for j in range(n):
                if i != j and group_labels_matrix[i, j] == "Close":
                    self.isomers[j].warning += "duplicate-"
                    visited.add(j)

        self.output_isomers = self.isomers
        return self.output_isomers

    def _analyze_similarity(self, matrix: np.ndarray, quantile: float = 0.2, method: str = "cluster", cutoff: Optional[float] = None) -> np.ndarray:
        """
        Analyze pairwise similarity between isomers using either MeanShift clustering or a hard cutoff.

        :param matrix: 2D numpy array representing the fingerprint difference matrix.
        :param quantile: Quantile for estimating bandwidth in MeanShift mode.
        :param method: 'meanshift' or 'cutoff' to define grouping method.
        :param cutoff: Optional float threshold to use when method='cutoff'.
        :return: 2D numpy array of group labels ('Close' or 'Far').
        """
        if method == "cluster":
            # Flatten the matrix values for clustering.
            # Upper triangle indices are used to avoid redundancy.
            triu_indices = np.triu_indices_from(matrix, k=1)
            values = matrix[triu_indices].reshape(-1, 1)
            bandwidth = estimate_bandwidth(values, quantile=quantile, n_samples=len(values))
            ms = MeanShift(bandwidth=0.5)
            ms.fit(values)
            cluster_labels = ms.labels_
            cluster_centers = ms.cluster_centers_

            # Identify the cluster whose center is closest to zero.
            close_cluster_idx = np.argmin(np.abs(cluster_centers))

            # Extract the true cutoff as the largest value still labeled "Close"
            close_values = values[np.where(cluster_labels == close_cluster_idx)]
            self.similarity_cutoff_used = float(np.max(close_values))  # ensure it's a float, not a 0-d array

            # Reconstruct group label matrix
            full_labels = np.full(matrix.shape, "Far", dtype=object)
            full_labels[triu_indices] = np.where(cluster_labels == close_cluster_idx, "Close", "Far")
            full_labels[(triu_indices[1], triu_indices[0])] = full_labels[triu_indices]  # mirror
            np.fill_diagonal(full_labels, "Far")  # enforce diagonal is "Far"

            return full_labels

        elif method == "cutoff":
            if cutoff is None:
                raise ValueError("Fatal Error: Cutoff value must be provided when using 'cutoff' method.")
            group_labels = np.where(matrix <= cutoff, "Close", "Far")
            self.similarity_cutoff_used = cutoff
            return group_labels

        else:
            raise ValueError(
                f"Fatal Error: Unsupported similarity analysis method '{method}'. "
                f"Supported options: 'cluster', 'cutoff'.")

    def _compute_fingerprint_matrix(self, isomers: list) -> np.ndarray:
        """
        Efficiently compute the upper-triangle of the fingerprint difference matrix.
        This result is then reflected to the lower triangle to create a full symmetric matrix.
        """
        # Generate fingerprints for each isomer
        fingerprints = [self._compute_sorted_distance_fingerprint(isomer) for isomer in isomers]
        # Ensure all fingerprints have the same length
        assert all(len(fp) == len(fingerprints[0]) for fp in fingerprints)
        # Initialize a square matrix to hold the differences
        n = len(fingerprints)
        diff_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):  # Only compute upper triangle
                diff = self._fingerprint_comparison(fingerprints[i], fingerprints[j], mode=self.isomer_comparison_mode)
                diff_matrix[i, j] = diff
                diff_matrix[j, i] = diff  # Reflect to lower triangle

        return diff_matrix

    @staticmethod
    def _compute_sorted_distance_fingerprint(isomer) -> np.ndarray:
        """
        Compute an inter-atomic distance matrix for an isomer and sort the distances under to two conditions:
        1. entries in the matrix (atom-atom distances) are sorted in ascending order.
        2. the rows of the matrix are sorted in lexicographical order of the element symbols.
            i.e. C-C, C-H, H-H, etc. This ensures like inter-atomic distances are compared between isomers.
        The resulting 1D array acts as a fingerprint for the isomer which can be compared to other isomers.

        :param isomer: ASE Atoms object.
        :return: sorted inter-atomic distance fingerprint as a 1D numpy array.
        """

        # Get the positions and symbols of the atoms in the isomer.
        positions = isomer.DART_atoms.get_positions()  # shape (N, 3)
        symbols = isomer.DART_atoms.get_chemical_symbols()  # list of element symbols

        # Dictionary to collect distances keyed by element pair, e.g. ('C', 'H')
        pair_distances = {}

        # Compute distances only for the upper triangle (i < j)
        for i in range(len(symbols)):
            for j in range(i + 1, len(symbols)):
                # Create a key that is independent of the order of the atoms i.e ('C', 'H') == ('H', 'C') becomes ('C', 'H')
                key = tuple(sorted((symbols[i], symbols[j])))
                d = np.linalg.norm(positions[i] - positions[j])
                # Checks if the key exists if not the key is created and the distance is appended to the list
                pair_distances.setdefault(key, []).append(d)

        # For a consistent fingerprint, the distances for each element pair is sorted
        # and then concatenate them in a fixed order (e.g., lexicographical order of the keys)
        fingerprint_parts = [
            d for key in sorted(pair_distances)
            for d in sorted(pair_distances[key])
        ]
        return np.array(fingerprint_parts)

    @staticmethod
    def _fingerprint_comparison(fp1: np.ndarray, fp2: np.ndarray, mode: str = "max_diff"):
        """
        Compute the maximum absolute difference between two sorted fingerprint vectors.
        :param fp1: 1D numpy array fingerprint for isomer 1.
        :param fp2: 1D numpy array fingerprint for isomer 2.
        :param mode: Comparison mode: 'max_diff', 'sum_diff', 'mean_diff', or 'rmsd'.
        :return: Maximum absolute difference.
        """
        print(f"Comparing fingerprints with mode '{mode}'.")
        if mode == "max_diff":
            return np.max(np.abs(fp1 - fp2))
        elif mode == "sum_diff":
            return np.sum(np.abs(fp1 - fp2))
        elif mode == "mean_diff":
            return np.mean(np.abs(fp1 - fp2))
        elif mode == "rmsd":
            return np.sqrt(np.mean((fp1 - fp2) ** 2))
        else:
            raise ValueError(
                f"Fatal Error: Unsupported fingerprint comparison mode '{mode}'. "
                f"Supported modes are 'max_diff', 'sum_diff', 'mean_diff', and 'rmsd'.")

    def plot_fingerprint_difference_matrix(self,
                                           write_svg: bool = True,
                                           plot_plotly: bool = True,
                                           colorscale_min: float = 0.0,
                                           colorscale_mid: float = 0.5,
                                           color_scale_max: float = 1.0,
                                           min_color: Optional[dict] = None,
                                           mid_color: Optional[dict] = None,
                                           max_color: Optional[dict] = None,
                                           cell_label_mode: str = "value",  # "value", "group", or "none"
                                           ) -> None:
        """
        Visualizes the fingerprint difference matrix as a heatmap using Matplotlib and/or Plotly.

        :param: write_svg: If True, saves the Matplotlib figure as an SVG file.
        :param: plot_plotly: If True, creates an interactive Plotly heatmap.
        :param: colorscale_min: Minimum value for the color scale.
        :param: colorscale_mid: Midpoint value for the color scale.
        :param: color_scale_max: Maximum value for the color scale.
        :param: min_color: Dictionary with RGB values for the minimum color.
        :param: mid_color: Dictionary with RGB values for the midpoint color.
        :param: max_color: Dictionary with RGB values for the maximum color.
        :param: cell_label_mode: 'value' (float), 'group' (Close/Far), or 'none'.
        """

        assert cell_label_mode in {"value", "group", "none"}, ValueError(
            f"Fatal Error: Unsupported cell label mode '{cell_label_mode}'. "
            "Supported modes are 'value', 'group', and 'none'.")

        # Fallback defaults for colors
        min_color = {"r": 238, "g": 100, "b": 97} if min_color is None else min_color
        mid_color = {"r": 255, "g": 255, "b": 255} if mid_color is None else mid_color
        max_color = {"r": 12, "g": 171, "b": 185} if max_color is None else max_color

        if self.diff_matrix is None:
            self.filter()

        labels_list = [f"{i}" for i in range(len(self.isomers))]
        df = pd.DataFrame(self.diff_matrix, index=labels_list, columns=labels_list)

        group_labels_matrix = self._analyze_similarity(
            matrix=self.diff_matrix,
            method=self.isomer_comparison_grouping_mode,
            cutoff=self.fingerprint_grouping_cutoff
        )

        # --- Matplotlib ---
        if write_svg:
            cmap = LinearSegmentedColormap.from_list("custom_rwb", [
                (colorscale_min, tuple(np.array(list(min_color.values())) / 255)),
                (colorscale_mid, tuple(np.array(list(mid_color.values())) / 255)),
                (color_scale_max, tuple(np.array(list(max_color.values())) / 255))
            ])
            fig, ax = plt.subplots(figsize=(8, 8))
            cax = ax.imshow(df.values, cmap=cmap, vmin=0, vmax=1.0)

            def get_contrasting_text_color(value, vmin=0.0, vmax=0.5):
                norm_val = (value - vmin) / (vmax - vmin)
                r, g, b = cmap(norm_val)[:3]
                luminance = 0.299 * r + 0.587 * g + 0.114 * b
                return 'black' if luminance > 0.5 else 'white'

            for i in range(len(df)):
                for j in range(len(df)):
                    label = ""
                    if cell_label_mode == "value":
                        label = f"{df.values[i, j]:.2f}"
                    elif cell_label_mode == "group":
                        label = group_labels_matrix[i, j]

                    if label:
                        color = get_contrasting_text_color(df.values[i, j])
                        ax.text(j, i, label, ha='center', va='center', color=color, fontsize=8)

            ax.set_xticks(np.arange(len(df)))
            ax.set_yticks(np.arange(len(df)))
            ax.set_xticklabels(labels_list, rotation=0, ha='center', fontsize=12)
            ax.set_yticklabels(labels_list, rotation=0, ha='right', fontsize=12)
            ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

            fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04, label="Difference")
            ax.set_title("Comparison Matrix (Matplotlib)")
            plt.tight_layout()
            plt.savefig("Isomer_comparison_matrix.svg", format="svg")
            plt.close()

        # --- Plotly ---
        if plot_plotly:
            color_scale = [
                [colorscale_min, f"rgb({min_color['r']},{min_color['g']},{min_color['b']})"],
                [colorscale_mid, f"rgb({mid_color['r']},{mid_color['g']},{mid_color['b']})"],
                [color_scale_max, f"rgb({max_color['r']},{max_color['g']},{max_color['b']})"]
            ]
            self.launch_interactive_heatmap(df, group_labels_matrix, color_scale, cell_label_mode)

    def launch_interactive_heatmap(self, df: pd.DataFrame, group_labels_matrix: np.ndarray,
                                   color_scale: list, cell_label_mode: str = "value"):
        """
        Launch an interactive Dash app to render the heatmap and enable isomer alignment viewing.
        """
        app = dash.Dash(__name__)
        fig = px.imshow(
            df,
            color_continuous_scale=color_scale,
            zmin=0, zmax=1.0,
            title="Comparison Matrix (Interactive)",
            labels={"x": "Isomer Index", "y": "Isomer Index", "color": "Difference"}
        )

        # Get the upper triangle indices (excluding the diagonal)
        triu_indices = np.triu_indices_from(df.values, k=1)

        # Extract the values from the upper triangle (excluding diagonal)
        upper_triangle_values = df.values[triu_indices]

        # Create the histogram
        fig_hist = px.histogram(
            x=upper_triangle_values,
            nbins=400,
            title="Histogram of Comparison Values (Upper Triangle Only)",
            labels={'x': 'Difference Value', 'y': 'Count'}
        )

        if self.similarity_cutoff_used is not None:
            fig_hist.add_vline(
                x=self.similarity_cutoff_used,
                line_width=3,
                line_dash="dash",
                line_color="red",
                annotation_text=f"Cutoff = {self.similarity_cutoff_used:.2f}",
                annotation_position="top left"
            )

        annotations = []
        for i in range(len(df)):
            for j in range(len(df)):
                value = df.values[i, j]
                label = ""
                if cell_label_mode == "value":
                    label = f"{value:.2f}"
                elif cell_label_mode == "group":
                    label = group_labels_matrix[i, j]

                if label:
                    # Compute contrast-aware color
                    r, g, b = [int(c) for c in color_scale[min(2, int(3 * value))][1][4:-1].split(',')]
                    luminance = 0.299 * r + 0.587 * g + 0.114 * b
                    font_color = "black" if luminance > 186 else "white"
                    annotations.append(dict(
                        x=j, y=i, text=label,
                        showarrow=False,
                        font=dict(color=font_color, size=10),
                        xanchor="center", yanchor="middle"
                    ))

        fig.update_layout(
            annotations=annotations,
            dragmode="zoom",
            hovermode="closest",
            height=800,
            width=800
        )

        app.layout = dash.html.Div([
            dash.dcc.Graph(id="heatmap", figure=fig),
            dash.dcc.Graph(id="histogram", figure=fig_hist),

        ])

        @app.callback(
            dash.Output("heatmap", "figure"),
            dash.Input("heatmap", "clickData")
        )
        def display_click(clickData):
            if clickData:
                point = clickData["points"][0]
                i_idx = int(point["y"])
                j_idx = int(point["x"])
                print(f"Clicked cell: ({i_idx}, {j_idx}) — triggering isomer viewer.")
                self.view_isomer_alignment(i_idx, j_idx)
            return fig  # Return same figure to avoid layout refresh

        app.run(debug=False, use_reloader=False)

    def view_isomer_alignment(self, index1: int, index2: int, grid_size: int = None) -> None:
        """
        Visualize two isomers (by index) and their optimal alignment using ASE's viewer.
        Frame 0: isomer1
        Frame 1: isomer2
        Frame 2: overlaid view after rotating isomer2 to align with isomer1.

        :param: grid_size: The number of grid points when scanning from 0 to 360 for exact alignment.
        :param: index1: Index of the reference isomer in self.isomers
        :param: index2: Index of the isomer to be aligned and overlaid
        """
        assert 0 <= index1 < len(self.isomers), f"Index1 out of range: {index1}"
        assert 0 <= index2 < len(self.isomers), f"Index2 out of range: {index2}"

        isomer1 = self.isomers[index1].DART_atoms.copy()
        isomer2 = self.isomers[index2].DART_atoms.copy()
        isomer2_aligned = self.isomers[index2].DART_atoms.copy()

        # Determine alignment axes
        if len(self.metal_centres) == 1:
            bounds = [[0, 360] for _ in range(3)]
            axes = np.eye(3)
        elif len(self.metal_centres) == 2:
            bounds = [[0, 360]]
            axis_vector = np.array(self.metal_centres[1]) - np.array(self.metal_centres[0])
            axis_vector /= np.linalg.norm(axis_vector)
            axes = [axis_vector]
        elif len(self.metal_centres) >= 3:
            bounds = None
            axes = None
            print("Three or more metal centres — skipping rotation; showing structures unaligned.")
        else:
            raise ValueError("Fatal Error: Invalid number of metal centres for alignment.")

        # Perform rotation if applicable
        if len(self.metal_centres) < 3:
            result_angles = brute(
                func=self.objective_function,
                ranges=bounds,
                args=(isomer1, isomer2_aligned, axes),
                Ns=grid_size if grid_size else self.grid_size
            )
            R_total = self.combined_rotation_matrix(result_angles, axes)
            self.apply_combined_rotation(atoms=isomer2_aligned, R_total=R_total, center=np.array(self.metal_centres[0]))

        # Create overlaid image: isomer1 + rotated isomer2
        overlaid = isomer1.copy() + isomer2_aligned.copy()

        # Assign tags so colors are distinguishable
        for atom in isomer1:
            atom.tag = 1
        for atom in isomer2:
            atom.tag = 2
        for atom in overlaid:
            atom.tag = 3

        view([isomer1, isomer2, overlaid], viewer="ase")


class IsomerClashFilter:
    """
    Filters out isomers that have atomic clashes based on covalent radii.
    Optionally considers ligand-metal clashes separately.
    """

    def __init__(
            self,
            isomers: List[Isomer],
            covalent_radii: dict,
            buffer: float = -0.3,
            check_metal_clashes: bool = True
    ):
        """
        :param isomers: List of ASE Atoms objects.
        :param covalent_radii: Dictionary mapping atomic symbols to covalent radii.
        :param buffer: Float to tune clash threshold (default -0.3 Å).
        :param check_metal_clashes: Whether to consider metal-ligand clashes (default True).
        """
        self.isomers = isomers
        self.cov_radii = covalent_radii
        self.buffer = buffer
        self.check_metal_clashes = check_metal_clashes
        self.non_clashing_isomers = []
        self.rejected_isomers = []
        print(f"Initialized ClashFilter with {len(self.isomers)} isomers, "
              f"buffer: {self.buffer}, "
              f"check_metal_clashes: {self.check_metal_clashes}")



    def filter(self) -> Tuple[List[Atoms], List[int]]:
        """
        Filters out isomers that have atomic clashes.
        :return: Tuple of (filtered isomers, indices of retained isomers in original list)
        """
        self.non_clashing_isomers = []
        self.rejected_isomers = []

        for idx, isomer in enumerate(self.isomers):
            atoms = isomer.DART_atoms
            if self.has_clash(atoms):
                isomer.warning += "clashing-"
                self.rejected_isomers.append(isomer)
                print(f"Clash detected in isomer [{idx}] of [{len(self.isomers)}], tagging as 'clash'.")
            else:
                self.non_clashing_isomers.append(isomer)

        print(f"Tagged {len(self.rejected_isomers)} isomers as 'clash'.")
        return self.isomers  # all isomers are returned

    def has_clash(self, atoms: Atoms) -> bool:
        """
        Determines whether any atomic clashes exist in the structure.
        Atom pairs are checked based on their arraus.

        :param atoms: ASE Atoms object with a information stored in array.
        :return: True if a clash is found, False otherwise.
        """
        if 'ligand_idx_tags' not in atoms.arrays:
            raise ValueError(
                "Fatal Error: 'ligand_idx_tags' array not found in atoms object. "
                "Ensure that the atoms object has been properly tagged.")

        ligand_tags = atoms.get_array('ligand_idx_tags')  # shape: (n_atoms,)
        n_atoms = len(atoms)
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                tag_i = ligand_tags[i]
                tag_j = ligand_tags[j]

                # Skip same-ligand checks
                if tag_i == tag_j:
                    continue

                # Skip ligand-metal checks if disabled
                if not self.check_metal_clashes and (tag_i == 0 or tag_j == 0):
                    continue

                dist = np.linalg.norm(positions[i] - positions[j])
                r1 = self.cov_radii.get(symbols[i], 1.0)
                r2 = self.cov_radii.get(symbols[j], 1.0)

                if dist < (r1 + r2 + self.buffer):
                    print(
                        f"Clash detected between {symbols[i]}{i} (tag {tag_i}) and "
                        f"{symbols[j]}{j} (tag {tag_j}) at distance {dist:.2f} Å."
                    )
                    return True

        return False

    def view_rejected_isomers(self):
        """
        View the rejected isomers in ASE GUI
        """
        if not self.rejected_isomers:
            print("No rejected isomers to view.")
            return

        rejected_atoms = [isomer.DART_atoms for isomer in self.rejected_isomers]
        view(rejected_atoms)


class AtomsCombiner:
    def __init__(self, base_atoms: Atoms, xyz_path: Union[str, pl.Path]):
        """
        Initialize with a base Atoms object and a path to an .xyz file.

        :param base_atoms: ASE Atoms object
        :param xyz_path: Path to .xyz file (str or pathlib.Path)
        """
        if not isinstance(base_atoms, Atoms):
            raise ValueError("Fatal Error: base_atoms must be an ASE Atoms object.")

        self.base_atoms = base_atoms
        self.xyz_path = pl.Path(xyz_path)
        self.xyz_atoms = self._load_xyz()



    def _load_xyz(self) -> Atoms:
        """Read the .xyz file and return an ASE Atoms object."""
        if not self.xyz_path.exists():
            return None
        atoms = ase.io.read(str(self.xyz_path))  # ASE does not accept Path directly
        if not isinstance(atoms, Atoms):
            raise ValueError(f"Fatal Error: The file {self.xyz_path} does not contain valid ASE Atoms data.")
        return atoms

    def combine(self) -> Atoms:
        """Return a new Atoms object with the base and xyz Atoms combined."""
        if self.xyz_atoms is None:
            print("No auxiliary structure found in .xyz file; returning base atoms only.")
            return self.base_atoms
        return self.base_atoms + self.xyz_atoms




class BiTransRotationModifier:
    def __init__(self, isomers: List[Isomer], opt: bool = True):
        """
        This class will take a list of ASE Atoms objects and optimize mono-coordinating ligands around their coordination axis
        """
        self.input_isomers = isomers
        self.opt_command = opt
        self.output_isomers = []
        print(f"AxialOpt initialized with {len(self.input_isomers)} Isomer objects.")



    def modify(self):
        """
        Optimize the rotation of each mono-coordinating ligand around their respective coordination axis simultaneously
        """
        # Set a random seed for reproducibility
        np.random.seed(42)

        # If opt_command is False, return the input isomers without optimization
        if not self.opt_command:
            return self.input_isomers

        # Loop through each of the inputted complexes
        for isomer in self.input_isomers:

            atoms = isomer.DART_atoms.copy()
            target_vectors = isomer.ligand_target_vectors
            ligand_origins = isomer.ligand_origins

            # Run the global optimizer.
            bounds = [[0, 360] for _ in target_vectors]
            result = differential_evolution(
                self.objective_function, bounds=bounds,
                args=(target_vectors, ligand_origins, atoms)
            )

            # Retrieve the relevant array
            ligand_idx_tags = atoms.get_array('ligand_idx_tags')
            n_denticity_tags = atoms.get_array('n_denticity_tags')

            unique_tags_indices_set = np.unique(ligand_idx_tags)
            unique_tags_indices_set = unique_tags_indices_set[unique_tags_indices_set != 0]

            # Loop through the unique second tags
            for angle, axis, origin, tag in zip(list(result.x), target_vectors, ligand_origins, unique_tags_indices_set):
                # Get indices where the second tag is the current tag
                indices = np.where(ligand_idx_tags == tag)[0]

                # Check if any of these indices have a first tag not equal to 1, if so, skip
                if np.any(n_denticity_tags[indices] != 1):
                    continue

                atoms = self.rotate(atoms=atoms, vector=np.array(list(axis.values())[0]), origin=np.array(origin), idc=indices, angle=angle).copy()
                print(f"Rotated ligand with tag {tag} by {angle:.2f} degrees around vector {axis} at origin {origin}.")

            # Create new Isomer from rotated structure
            new_dart = Isomer(
                atoms=atoms.copy(),
                ligands=isomer.ligands,
                ligand_target_vectors=isomer.ligand_target_vectors,
                ligand_origins=isomer.ligand_origins,
                metal_centers=isomer.metal_centers,
                default_graph=isomer.default_graph
            )

            # Append the new Isomer to the output complexes
            self.output_isomers.append(new_dart)

        print(f"Optimized {len(self.output_isomers)} complexes with mono-coordinating ligand rotations.")
        return self.output_isomers

    def visualize_structures(self):
        """
        Visualize input and output complexes interleaved using ASE's GUI.
        Each input is followed by its corresponding optimized output.
        """
        if not self.output_isomers:
            print("No output complexes found. Run opt_mono_rotation() first.")
            return

        structures_to_view = []

        for i, (in_isomer, out_isomer) in enumerate(zip(self.input_isomers, self.output_isomers)):
            input_copy = in_isomer.DART_atoms.copy()
            input_copy.info["label"] = f"Input {i}"

            output_copy = out_isomer.DART_atoms.copy()
            output_copy.info["label"] = f"Optimized {i}"

            structures_to_view.extend([input_copy, output_copy])  # interleave input/output

        print(f"Launching viewer for {len(structures_to_view)} structures...")
        view(structures_to_view)

    def objective_function(self, x: np.ndarray, vectors_in: List[np.array], origins_in: List[np.array], TMC_in: Atoms, ):
        """
        Objective function to optimize the position of the ligands in the TMC complex.
        :param: x:  Array of angles to rotate each ligand around its respective vector.
        :param: vectors_in: List of vectors for each ligand, where each vector is a dictionary with a single key-value pair.
        :param: origins_in: List of origins for each ligand, where each origin is a list of three floats.
        :param: TMC_in: ASE Atoms object representing the transition metal complex (TMC) to be optimized.
        :return: float: The penalty score based on interatomic distances.
        """
        # Generate a copy of the input complex
        TMC_worker = TMC_in.copy()

        # Retrieve the relevant array
        ligand_tags = TMC_worker.get_array("ligand_idx_tags")
        denticity_tags = TMC_worker.get_array("n_denticity_tags")

        # Get the unique set of all nonzero second tags (each unique tag represents a unique ligand, a zero tag represents a metal).
        unique_tags = [tag for tag in np.unique(ligand_tags) if tag != 0]

        # Loop through each unique tag and apply the rotation
        for tag, angle, axis, origin in zip(unique_tags, list(x), vectors_in, origins_in):

            # Get indices where the second tag in the list is equal to the current "tag" (essentially the indices of the atoms in this particular ligand)
            indices = np.where(ligand_tags == tag)[0]

            # Check if any of these indices have a first tag not equal to 1, if so, skip
            # This ensures only ligands that have a effective coordination number of 1 are rotated
            if np.any(denticity_tags[indices] != 1):
                continue  # Skip this ligand group

            TMC_worker = self.rotate(atoms=TMC_worker, vector=np.array(list(axis.values())[0]), origin=np.array(origin), idc=indices, angle=angle).copy()

        # Get the interatomic distance matrix
        distance_matrix = TMC_worker.get_all_distances()

        # Set the diagonal to a large number to avoid self-interaction (or np.inf)
        np.fill_diagonal(distance_matrix, np.inf)

        # Calculate the penalty: for each pair with d <= 4, add 1/d^2 to the penalty
        penalty = np.sum(1.0 / (distance_matrix ** 2))

        return penalty

    @staticmethod
    def rotate(atoms: Atoms, vector: np.array, origin: np.array, idc: List[int], angle: int):
        """
        Rotate the atoms in the Atoms object (only atoms with indices=idc) around the vector by the specified angle.
        :param atoms: Atoms object to rotate.
        :param vector: vector to rotate around.
        :param origin: origin of the rotation.
        :param idc: indices of the atoms to rotate.
        :param angle: the angle to rotate the atoms by in degrees.
        :return: an ase.Atoms object with the rotated atoms.
        """

        # Normalize rotation vector
        vector = np.asarray(vector, dtype=float)
        vector /= np.linalg.norm(vector)

        # Create rotation object
        rotation = R.from_rotvec(np.radians(angle) * vector)

        # Copy the atoms object to avoid modifying the original
        rotated_atoms = atoms.copy()

        # Apply rotation to selected atoms
        for i in idc:
            pos = atoms.positions[i] - origin  # Translate to origin
            rotated_pos = rotation.apply(pos) + origin  # Rotate and translate back
            rotated_atoms.positions[i] = rotated_pos

        return rotated_atoms

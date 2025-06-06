"""
Class for extracting ligands from a database of complexes.
"""
import functools
import re
import warnings

import jsonlines
import networkx as nx
import pandas as pd
from copy import deepcopy

import pysmiles
from ase import Atoms

from DARTassembler.src.constants.chem import Element
from DARTassembler.src.metalig.db import LigandDB, BaseDB
from pathlib import Path
from tqdm import tqdm
import json
from typing import Union
import numpy as np
from datetime import datetime
from collections import Counter
from DARTassembler.src._extraction.dataloader import DataLoader
from DARTassembler.src.metalig.mol import Ligand, BaseMolecule
from DARTassembler.src.metalig.utils_graph import get_only_complex_graph_connected_to_metal, get_graph_fragments, \
    get_sorted_atoms_and_indices_from_graph, find_node_in_graph_by_label, get_graph_hash, view_graph, \
    graph_to_dict_with_node_labels, graph_from_graph_dict
from DARTassembler.src.misc.io import NumpyEncoder, iterate_over_json, get_n_entries_of_json_db
from DARTassembler.src._extraction.utils import get_charges_of_unique_ligands, update_ligand_with_charge_inplace, \
    Composition
from DARTassembler.src.metalig.utils import get_duration_string, series2namedtuple, identify_metal_in_ase_mol, \
    make_None_to_NaN, is_between, update_dict_with_warning_inplace
from DARTassembler.src._extraction.test_constants import CHARGE_BENCHMARKED_COMPLEXES
from DARTassembler.src._extraction.constants import odd_n_electrons_warning, similar_molecule_with_diff_n_hydrogens_warning
from collections import defaultdict
# from memory_profiler import profile as mem_profile
# from test.profiling import profile as line_profile
# from test.profiling import print_stats

class LigandExtraction:

    def __init__(self, database_path: str,
                 data_store_path: str,
                 exclude_not_fully_connected_complexes: bool = True,
                 testing: Union[bool, int] = False,
                 graph_strat: str = "default",
                 exclude_charged_complexes: bool = False,
                 only_complexes_with_os: bool = False,
                 unique_ligand_id: str = 'graph_hash_with_metal',
                 store_database_in_memory: bool = False,
                 ):

        self.ligand_to_unique_ligand = None
        self.unique_ligand_info_props = None
        self.grouped_ligands = None
        self.database_path = ""
        self.data_store_path = ""
        self.exclude_not_fully_connected_complexes = None
        self.testing = None

        self.df_full_ligand_db = None
        self.df_unique_ligand_db = None
        self.df_complex_db = None

        self.test_complexes = None
        self.n_pred_charges = None
        self.store_database_in_memory = store_database_in_memory
        self.exclude_charged_complexes = exclude_charged_complexes
        self.only_complexes_with_os = only_complexes_with_os
        self.test_complexes = CHARGE_BENCHMARKED_COMPLEXES if not testing == False else []
        self.unique_ligand_id = unique_ligand_id

        self.excluded_complex_ids = defaultdict(list)
        self.graph_strat = graph_strat

        self.check_and_set_init_input(
            database_path=database_path,
            data_store_path=data_store_path,
            exclude_not_fully_connected_complexes=exclude_not_fully_connected_complexes,
            testing=testing,
        )

        # Define properties which will be passed to the Linear Charge Solver
        self.LCS_needed_complex_props = ['mol_id', 'metal_oxi_state', 'charge', 'ligand_names']
        self.LCS_needed_ligand_props = ['stoichiometry', 'name', 'n_protons', 'unique_name']

        self.input_complexes_json = Path(self.data_store_path, 'tmQMG.json')
        self.output_complexes_json = Path(self.data_store_path, 'complex_db.json')
        self.tmp_output_complexes_json = Path(self.data_store_path, 'TEMPORARY_complex_db.json')
        self.unique_ligands_json = Path(data_store_path, 'tmQM_Ligands_unique.json')
        self.full_ligands_json = Path(data_store_path, 'tmQM_Ligands_full.json')

    def check_and_set_init_input(self,
                                 database_path: str,
                                 data_store_path: str,
                                 exclude_not_fully_connected_complexes: bool,
                                 testing: Union[bool, int]
                                 ):
        database_path = Path(str(database_path))
        data_store_path = Path(str(data_store_path))

        if not database_path.exists():
            raise ValueError(f'Path given as `database_path` doesn\'t exist: {database_path}')

        if not data_store_path.exists():
            print(f'Path given as `data_store_path` ({data_store_path} doesn\'t exist yet. Make new directory.')
            data_store_path.mkdir(parents=True, exist_ok=True)

        if not (isinstance(testing, int) or isinstance(testing, bool) or isinstance(testing, list)):
            raise ValueError(f'Input variable `testing` must be int or bool but is {type(testing)}.')


        self.database_path = database_path
        self.data_store_path = data_store_path
        self.exclude_not_fully_connected_complexes = exclude_not_fully_connected_complexes
        self.testing = testing

        return

    def reorder_input_complexes(self, db_dict: dict, first_complexes: list):
        """
        Reorders the complexes in db_dict so that the complexes specified in first_complexes are the first ones in the dictionary and the others follow. This is useful to include specific complexes when running a test run with less complexes and doing a charge benchmark, so that the complexes which are benchmarked definitely appear.
        This was implemented using this reordering so that if self.testing is a number less than the number of specified complexes, it would still be respected and not all specified complexes would be included.
        :param db_dict: dictionary of complexes
        :param first_complexes: list of complexes to be shifted to the front of the dictionary
        """
        first_complexes = [c_id for c_id in first_complexes if c_id in set(db_dict.keys())]
        print(f'Include {len(first_complexes)} specified test complexes.')

        complexes_not_in_test_complexes = [c_id for c_id in db_dict.keys() if not c_id in set(first_complexes)]
        new_complex_order = first_complexes + complexes_not_in_test_complexes
        db_dict = {c_id:db_dict[c_id] for c_id in new_complex_order}

        return db_dict


    def load_input_data_to_json(self,
                                overwrite_atomic_properties: bool = False,
                                **kwargs):
        """
        Establish and safe the database as json for simple loading.
        """
        db_dict = DataLoader(database_path_=self.database_path, overwrite=overwrite_atomic_properties).data_for_molDB

        # Reorder complexes so that the benchmarked complexes are the first ones, to get a good statistic even when testing only with few complexes
        if self.test_complexes:
            db_dict = self.reorder_input_complexes(db_dict=db_dict, first_complexes=self.test_complexes)

        input_complex_db = ComplexDB.from_json(
                                            json_=db_dict,
                                            type_="Complex",
                                            max_number=self.testing,
                                            graph_strategy=self.graph_strat,
                                            **kwargs
                                            )

        input_complex_db._to_json(path=self.input_complexes_json, json_lines=True)

        return

    def print_excluded_complexes(self):
        print('Excluded complexes:')
        for reason, c_ids in self.excluded_complex_ids.items():
            print(f'    - {reason}: {len(c_ids)} ({len(c_ids) / self.n_input_complexes_before_filtering*100:.2g}%)')

        print(f'  New number of input complexes: {self.n_complexes}')

        return

    def prefilter_if_complex_valid(self, c_id, c) -> bool:
        """
        Filters input complexes in `self.complex_db` by multiple criteria without needing information about ligands.
        """
        if c.n_donors == 0:
            self.excluded_complex_ids['Metal ion'].append(c_id)
            return False

        if c.has_fragment(frag='O'):
            self.excluded_complex_ids['Has unconnected O'].append(c_id)
            return False

        if c.has_fragment(frag='H'):
            self.excluded_complex_ids['Has unconnected H'].append(c_id)
            return False

        if c.has_fragment(frag=['O', 'O']):
            self.excluded_complex_ids['Has unconnected O2'].append(c_id)
            return False

        if c.has_fragment(frag=['H', 'O']):
            self.excluded_complex_ids['Has unconnected OH'].append(c_id)
            return False

        if not c.global_props['is 3d']:
            self.excluded_complex_ids['Is not 3D'].append(c_id)
            return False

        if not c.has_consistent_stoichiometry_with_CSD():
            self.excluded_complex_ids['Inconsistent CSD stoichiometry'].append(c_id)
            return False

        # if not c.has_consistent_stoichiometry_with_smiles(smiles=c.global_props['smiles'], ignore_element_count=True, print_warnings=False):
        #     self.excluded_complex_ids['Inconsistent smiles elements'].append(c_id)
        #     return False

        if c.global_props['smiles'] is None and not c.has_unknown_bond_orders:
            self.excluded_complex_ids['No smiles without bad bonds'].append(c_id)
            return False

        # if not c.complex_is_biggest_fragment(allow_complexes_greater_than=10):
        #     self.excluded_complex_ids['Complex is counter ion'].append(c_id)
        #     return False

        if not 'H' in c.atomic_props['atoms']:
            self.excluded_complex_ids['Complex has no H'].append(c_id)
            return False

        if not 'C' in c.atomic_props['atoms']:
            self.excluded_complex_ids['Complex has no C'].append(c_id)
            return False

        if c._count_atoms_with_n_bonds(element='C', n_bonds=1) > 0:
            self.excluded_complex_ids['C atom with only 1 bond'].append(c_id)
            return False

        min_dist, _, _ = c.get_interatomic_distances()
        if min_dist < 0.5:
            self.excluded_complex_ids['Atoms closer than 0.5A'].append(c_id)
            return False

        min_dist, _, _ = c.get_interatomic_distances(skip_elements='H')
        if min_dist < 0.85:
            self.excluded_complex_ids['Heavy atoms closer than 0.85A'].append(c_id)
            return False

        if self.only_complexes_with_os and not c.has_metal_os():
            self.excluded_complex_ids['No metal OS'].append(c_id)
            return False

        if self.exclude_not_fully_connected_complexes and not c.fully_connected:
            self.excluded_complex_ids['Not fully connected'].append(c_id)
            return False

        if self.exclude_charged_complexes and c.charge != 0:
            self.excluded_complex_ids['Is charged'].append(c_id)
            return False

        return True

    def postfilter_if_ligands_valid(self, c_id: str, comp: RCA_Complex) -> bool:
        """
        Function for filtering complexes after the ligand extraction. Return False to exclude that complex.
        """
        if not hasattr(comp, 'ligands'):
            raise ValueError(f'The complex {comp.mol_id} has no attribute `ligands`.')
        elif comp.ligands == []:
            raise ValueError(f'The ligand list of the complex {comp.mol_id} is empty.')

        if comp.count_ligands_with_stoichiometry(atoms=['O'], only_connected=True) >= 3:
            self.excluded_complex_ids['More than 3 O ligands'].append(c_id)
            return False

        if comp.count_ligands_with_stoichiometry(atoms=['N'], only_connected=True) >= 2:
            self.excluded_complex_ids['More than 2 N ligands'].append(c_id)
            return False

        if comp.count_ligands_with_stoichiometry(atoms=['C']) > 0:
            self.excluded_complex_ids['Ligand which is just C'].append(c_id)
            return False

        if comp.count_n_unconnected_ligands(max_n_atoms=1) > 5:
            self.excluded_complex_ids['More than 5 unconnected ligands'].append(c_id)
            return False

        if comp._count_coordinating_atoms_with_distance_to_metal_greater_than(distance=1.9, element='O', max_n_atoms=1) > 0:
            self.excluded_complex_ids['Oxygen ligand more than 1.9A away from metal'].append(c_id)
            return False

        alkalis = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
        if comp.count_atoms_in_ligands(atoms=alkalis, only_if_connected_to_metal=True) > 0:
            self.excluded_complex_ids['Alkali metal in ligand'].append(c_id)
            return False

        noble_gases = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'Og']
        if comp.count_atoms_in_ligands(atoms=noble_gases, only_if_connected_to_metal=True) > 0:
            self.excluded_complex_ids['Noble gas in ligand'].append(c_id)
            return False

        heavy_metals = ['Tl', 'Pb', 'Bi', 'Po', 'Nh', 'Fl', 'Mc', 'Lv']
        if comp.count_atoms_in_ligands(atoms=heavy_metals, only_if_connected_to_metal=True) > 0:
            self.excluded_complex_ids['Heavy metal in ligand'].append(c_id)
            return False

        # Exclude likely cage structures of metals which are no real ligands.
        metals = ['B', 'Al', 'Ga', 'In', 'Tl', 'Nh', 'Si', 'Ge', 'Sn', 'Pb', 'Fl', 'As', 'Sb', 'Bi', 'Mc', 'Te', 'Po', 'Lv']
        if comp.count_ligands_containing_only(atoms=metals, denticity_range=[1, np.inf], n_atoms_range=[2, np.inf], except_elements=['H']):
            self.excluded_complex_ids['Metal cage structures'].append(c_id)
            return False

        return True

    def get_ligand_class_properties_of_complex(self, complex_, props: list) -> list:
        """
        Returns a list of dicts with the specified properties for all ligands of this complex.
        @param complex_: RCA_Complex with ligands for which to get the properties
        @param props: list of properties. Must be the exact name of a property of the Ligand class.
        @return: list of dicts of all specified properties
        """""
        ligand_infos = []
        for lig in complex_.ligands:
            infos = {prop: value for prop, value in vars(lig).items() if prop in props}
            ligand_infos.append(infos)

        return ligand_infos

    def extract_ligands(self):
        """
        Extracts ligands from complexes and writes the complexes with ligands to an intermediate json file.
        """
        # Properties of the ligand which should be recorded in the global ligand dataframe.
        important_ligand_props = ['name', 'stoichiometry', 'n_protons', 'graph_hash', 'n_donors', 'has_all_bond_orders_valid',
                                    self.unique_ligand_id, 'original_metal_symbol', 'original_metal_os',
                                    'was_connected_to_metal', 'parent_complex_id', 'n_hydrogens', 'parent_complex_id',
                                    'heavy_atoms_graph_hash_with_metal']

        df_full_ligand_db = []
        df_complex_db = []
        if self.store_database_in_memory:
            self.complex_db = {}
        self.n_input_complexes_before_filtering = 0

        with jsonlines.open(self.tmp_output_complexes_json, mode='w', dumps=functools.partial(json.dumps, cls=NumpyEncoder)) as complex_writer:
            for idx, (csd_code, comp_dict) in tqdm(enumerate(iterate_over_json(self.input_complexes_json, show_progress=False)), desc='Extracting ligands from complexes'):
                if not self.testing or idx < self.testing:
                    self.n_input_complexes_before_filtering += 1

                    # Make complex class from dict
                    comp = RCA_Complex.from_dict(dict_=comp_dict)

                    # Normalize the complex
                    comp._remove_node_features_from_molecular_graphs_inplace()
                    comp._normalize_multigraph_into_graph_inplace()

                    # Filter out bad complexes
                    # All complexes which are filtered out are added to `self.excluded_complex_ids` with the reason why they were excluded.
                    complex_is_valid = self.prefilter_if_complex_valid(c_id=csd_code, c=comp)

                    if complex_is_valid:
                        # Extract ligands from the complex and add as property `comp.ligands`
                        comp.de_assemble()

                        # Do another filtering step, this time using information from the ligands.
                        complex_ligands_are_valid = self.postfilter_if_ligands_valid(c_id=csd_code, comp=comp)

                        if complex_ligands_are_valid:

                            # Record important information for dataframes of complexes and ligands
                            ligand_infos = self.get_ligand_class_properties_of_complex(
                                                                                        complex_=comp,
                                                                                        props=important_ligand_props
                                                                                        )
                            df_full_ligand_db.extend(ligand_infos)
                            df_complex_db.append({
                                                        'mol_id': comp.mol_id,
                                                        'stoichiometry': comp.stoichiometry,
                                                        'metal_oxi_state': comp.metal_oxi_state,
                                                        'charge': comp.charge,
                                                        'ligand_names': [lig.name for lig in comp.ligands],
                                                        'metal': comp.metal,
                                                    })

                            if not self.store_database_in_memory:
                                comp._append_to_file(key=csd_code, writer=complex_writer)    # Write to jsonlines file
                            else:
                                self.complex_db[csd_code] = comp        # Store in memory


        if self.store_database_in_memory:
            self.complex_db = ComplexDB(self.complex_db)
            self.full_ligand_db = self.build_full_ligand_db(copy=False)

        # Important: make dataframes of important information for ligands, complexes, unique ligands.
        # That way, we don't need to keep the entire database in memory but can still access the information we need.
        # ligands
        self.df_full_ligand_db = pd.DataFrame(df_full_ligand_db).set_index('name', drop=False)
        self.n_full_ligands = len(self.df_full_ligand_db)
        # complexes
        self.df_complex_db = pd.DataFrame(df_complex_db).set_index('mol_id', drop=False)
        self.n_complexes = len(self.df_complex_db)
        # unique ligands
        self.df_unique_ligand_db = self.get_unique_ligand_df()
        self.n_unique_ligands = len(self.df_unique_ligand_db)

        # Get useful mapping from name to unique name for later use
        self.ligand_to_unique_ligand = self.get_ligand_to_unique_ligand_dict()

        # Add `unique name` column to full ligand df
        df_ligand_to_unique_ligand = pd.DataFrame.from_dict(self.ligand_to_unique_ligand, orient='index', columns=['unique_name'])
        self.df_full_ligand_db = pd.merge(self.df_full_ligand_db, df_ligand_to_unique_ligand, left_index=True, right_index=True, how='left')
        # Add `n_ligand_instances` column to full ligand df
        self.df_full_ligand_db = pd.merge(self.df_full_ligand_db, self.df_unique_ligand_db['n_ligand_instances'], left_on='unique_name', right_index=True, how='left')

        self.print_excluded_complexes()

        return

    # @profile
    def get_unique_ligand_df(self) -> pd.DataFrame:
        """
        Returns a dataframe with all unique ligands from the full ligand db.
        """
        self.grouped_ligands = self.group_same_ligands()
        df_same_graph_hash = self.df_full_ligand_db.loc[:, ['graph_hash', 'n_donors', 'graph_hash_with_metal']].groupby('graph_hash').agg(lambda x: x.unique().tolist())
        unique_ligands = {}
        for grouped_ligand in tqdm(self.grouped_ligands.to_dict(orient='index').values(), desc="Building unique ligand dataframe"):
            same_ligand_names = grouped_ligand['name']
            name = self.choose_unique_ligand_representative_from_all_same_ligands(same_ligands=same_ligand_names)
            uname = 'unq_' + name

            unique_ligands[uname] = self.df_full_ligand_db.loc[name].to_dict()
            unique_ligands[uname].update({
                                            'unique_name': uname,
                                            'same_ligand_names': same_ligand_names
                                            })


            # Add useful statistical information of all ligands for this unique ligand
            graph_hash = unique_ligands[uname]['graph_hash']
            denticities = df_same_graph_hash.loc[graph_hash, 'n_donors']
            metals = self.df_full_ligand_db.loc[same_ligand_names, 'original_metal_symbol']
            count_metals = metals.value_counts().sort_values(ascending=False).to_dict()
            n_graph_hashes = len(df_same_graph_hash.loc[graph_hash, 'graph_hash_with_metal'])
            assert not 0 in denticities, 'The n_donors for unconnected ligands is assumed to be -1 but here there appears a 0.'
            has_unconnected_ligands = -1 in denticities
            unique_ligand_infos = {
                                    'n_ligand_instances': len(same_ligand_names),
                                    'same_graph_denticities': denticities,
                                    'metal_counts': count_metals,
                                    'n_same_graph_denticities': len(denticities),
                                    'n_metals': len(count_metals),
                                    'n_same_graphs': n_graph_hashes,
                                    'has_unconnected_ligands': has_unconnected_ligands,
                                    'all_ligands_metals': metals.tolist(),
                                    }
            unique_ligands[uname].update(unique_ligand_infos)

        df = pd.DataFrame.from_dict(unique_ligands, orient='index')
        # for updating the ligands from complex and full ligands db later
        self.unique_ligand_info_props = list(unique_ligand_infos.keys())

        return df

    def get_ligand_to_unique_ligand_dict(self) -> dict:
        """
        Returns a dict with ligand names as keys and unique ligand names as values.
        """
        ligand_to_unique_ligand = {}
        for ulig in self.df_unique_ligand_db.itertuples():
            uname = ulig.unique_name
            for name in ulig.same_ligand_names:
                ligand_to_unique_ligand[name] = uname

        return ligand_to_unique_ligand

    def build_full_ligand_db(self, copy: bool=True):
        full_ligands = {}
        for c in tqdm(self.complex_db.db.values(), 'Build full ligand db'):
            for lig in c.ligands:
                name = lig.name
                if copy:
                    full_ligands[name] = deepcopy(lig)
                else:
                    full_ligands[name] = lig
        full_ligand_db = LigandDB(full_ligands)
        return full_ligand_db

    def save_full_ligand_db(self):
        self.full_ligand_db._to_json(self.full_ligands_json)

        return

    def group_same_ligands(self, groupby: Union[str, None] = None) -> list:
        """
        Groups all ligands into groups with the same unique ligand.
        @param groupby: list of properties in self.basic_ligand_infos. If None, the default is used.
        @return: dataframe of grouped ligands
        """
        if groupby is None:
            groupby = self.unique_ligand_id

        grouped_ligands = self.df_full_ligand_db.groupby(groupby, sort=False).agg(list)

        return grouped_ligands


    def choose_unique_ligand_representative_from_all_same_ligands(self,
                                                                  same_ligands,
                                                                  strategy='good_bond_orders',
                                                                  ) -> str:
        if isinstance(same_ligands, dict):
            same_ligands = list(same_ligands.keys())

        if strategy == 'first':
            # Just take the first entry.
            name = same_ligands[0]
        elif strategy == 'good_bond_orders':
            # Preferably take ligands which have good bond orders
            name = same_ligands[0]
            same_ligand_props = self.df_full_ligand_db['has_all_bond_orders_valid'].loc[same_ligands]
            for lig_name, has_good_bond_orders in same_ligand_props.items():
                if has_good_bond_orders:
                    name = lig_name
                    break
        else:
            raise ValueError(
                f'Unknown strategy `{strategy}` to choose the unique ligand representative from all same ligands.')

        return name

    def get_unique_ligand_from_ligand(self, ligand: Ligand) -> Ligand:
        """
        Returns the unique ligand, given a normal ligand. The unique ligand has some additional and some deleted properties in contrast to the normal ligand.
        @param ligand: normal ligand
        @return: unique ligand
        """
        ulig = deepcopy(ligand)

        uname = self.ligand_to_unique_ligand[ligand.name]
        df_same_ligands = self.df_full_ligand_db.loc[self.df_full_ligand_db['unique_name'] == uname]
        same_ligand_names = df_same_ligands['name'].tolist()

        ulig.unique_name = uname
        ulig.all_ligand_names = same_ligand_names

        identical_ligand_info = {}
        identical_ligand_info['name'] = same_ligand_names
        identical_ligand_info['original_metal_symbol'] = df_same_ligands['original_metal_symbol'].tolist()
        identical_ligand_info['original_metal_os'] = df_same_ligands['original_metal_os'].tolist()
        original_complex_ids = df_same_ligands['parent_complex_id'].tolist()
        identical_ligand_info['original_complex_charge'] = self.df_complex_db.loc[original_complex_ids, 'charge'].tolist()
        identical_ligand_info['parent_complex_id'] = df_same_ligands['parent_complex_id'].tolist()
        ulig.identical_ligand_info = identical_ligand_info

        # Add information which is costly to calculate and therefore calculated only for unique ligands instead of all ligands
        # ulig.has_planar_donor_atoms = ulig.planar_check()

        # Set unique ligand stats information
        props = self.df_unique_ligand_db.loc[uname, self.unique_ligand_info_props].to_dict()
        for prop, val in props.items():
            setattr(ulig, prop, val)

        # Delete attributes which make sense only for ligands but not for unique ligands
        del ulig.is_chosen_unique_ligand

        ulig.same_graph_charges = self.all_charges_with_same_graph_hash[ulig.graph_hash]
        ulig.n_pred_charges = len(ulig.same_graph_charges)

        most_common_n_H = self.all_most_common_n_H[ulig.heavy_atoms_graph_hash_with_metal]
        ulig.common_graph_with_diff_n_hydrogens = bool(most_common_n_H != ulig.n_hydrogens)
        if ulig.common_graph_with_diff_n_hydrogens:
            ulig.add_warning(similar_molecule_with_diff_n_hydrogens_warning)

        ulig.n_electrons = ulig.n_protons - ulig.charge
        ulig.odd_n_electron_count = bool(ulig.n_electrons % 2 == 1)
        if ulig.odd_n_electron_count:
            ulig.add_warning(odd_n_electrons_warning)

        ulig.has_warnings = bool(len(ulig.warnings) > 0)

        return ulig

    @staticmethod
    def update_ligand_with_unique_ligand_information_inplace(
                                                                lig,
                                                                ulig,
                                                                share_properties=None,
                                                                collect_properties=None
                                                                ):
        if collect_properties is None:
            collect_properties = {}
        if share_properties is None:
            share_properties = []

        for prop in share_properties:
            value = deepcopy(getattr(ulig, prop))
            setattr(lig, prop, value)

        # Collect properties from unique ligand in a dictionary in the full ligands.
        for new_prop, old_props in collect_properties.items():
            info_dict = {prop: deepcopy(getattr(ulig, prop)) for prop in old_props}
            setattr(lig, new_prop, info_dict)

        lig.is_chosen_unique_ligand = ulig.name == lig.name

        return

    def ensure_complex_db(self):
        try:
            self.complex_db
        except AttributeError:
            self.complex_db = ComplexDB.from_json(self.output_complexes_json)

        return

    def ensure_unique_ligand_db(self):
        try:
            self.unique_ligand_db
        except AttributeError:
            self.unique_ligand_db = LigandDB.from_json(self.unique_ligands_json)

        return

    def ensure_full_ligand_db(self):
        try:
            self.full_ligand_db
        except AttributeError:
            self.full_ligand_db = LigandDB.from_json(self.full_ligands_json)

        return

    def update_complex_db_with_information(self,
                                           share_properties: list = [],
                                           collect_properties: dict = {}
                                           ):
        self.ensure_complex_db()
        charges = self.df_unique_ligand_db.to_dict(orient='index')

        # Update ligands with unique ligand information
        for c in tqdm(self.complex_db.db.values(), 'Update complex db with unique ligand information'):
            for lig in c.ligands:
                uname = self.ligand_to_unique_ligand[lig.name]
                ulig = series2namedtuple(self.df_unique_ligand_db.loc[uname])
                self.update_ligand_with_unique_ligand_information_inplace(
                    lig=lig,
                    ulig=ulig,
                    share_properties=share_properties,
                    collect_properties=collect_properties
                )
                update_ligand_with_charge_inplace(lig, charges=charges)

            # Update global props with some useful information
            c.global_props['n_ligands'] = len(c.ligands)
            c.global_props['n_unique_ligands'] = len(set(lig.unique_name for lig in c.ligands))
            n_ligands_occurring_once = sum(
                [lig.unique_ligand_information['n_ligand_instances'] == 1 for lig in c.ligands])
            c.global_props['n_ligands_occurring_once'] = n_ligands_occurring_once
            c.global_props['frac_ligands_occurring_once'] = n_ligands_occurring_once / len(c.ligands)

        return

    def update_databases_with_charges(self, df_ligand_charges: pd.DataFrame):
        charges = df_ligand_charges.set_index('unique_name').to_dict(orient='index')

        self.ensure_unique_ligand_db()
        not_intersecting_ligands = set(self.unique_ligand_db.db.keys()).symmetric_difference(set(charges.keys()))
        print(f'Charges could not be assigned due to missing OS: {len(not_intersecting_ligands)}/{len(self.unique_ligand_db.db)}')

        for ulig in self.unique_ligand_db.db.values():
            update_ligand_with_charge_inplace(ulig, charges)

        self.ensure_complex_db()
        for c in self.complex_db.db.values():
            for lig in c.ligands:
                update_ligand_with_charge_inplace(lig, charges)

        return

    def ensure_input_complex_db_exists(self,
                                       overwrite_atomic_properties: bool,
                                       use_existing_input_json: bool,
                                       **kwargs
                                       ):
        if use_existing_input_json:
            if not self.input_complexes_json.exists():
                print(
                    f'WARNING: Cannot use existing input json of complexes because path not found: {self.input_complexes_json}. Reload xzy, global properties and graph data instead.')
                self.load_input_data_to_json(overwrite_atomic_properties=overwrite_atomic_properties, **kwargs)
            else:
                # Check if the existing input json contains enough complexes for the testing parameter
                n_complexes_in_json = get_n_entries_of_json_db(self.input_complexes_json)

                if not self.testing is None and (self.testing > n_complexes_in_json or self.testing == False):
                    print(
                        f'WARNING: Cannot use existing input json of complexes because it contains less complexes than the testing parameter. Reload xzy, global properties and graph data instead.')
                    self.load_input_data_to_json(overwrite_atomic_properties=overwrite_atomic_properties, **kwargs)
        else:
            self.load_input_data_to_json(overwrite_atomic_properties=overwrite_atomic_properties, **kwargs)

        return

    def get_complex_dict_for_LCS(self) -> dict:
        """
        Returns a dictionary of all complexes with only the needed properties for the Linear Charge Solver.
        """
        charge_complexes = self.df_complex_db[self.LCS_needed_complex_props].to_dict(orient='index')

        for c_id, c in charge_complexes.items():
            charge_complexes[c_id]['ligands'] = []
            for name in c['ligand_names']:
                lig = self.df_full_ligand_db.loc[name]
                lig_props = {prop: lig.loc[prop] for prop in self.LCS_needed_ligand_props}
                charge_complexes[c_id]['ligands'].append(lig_props)
            # Delete ligand names because they are not needed anymore
            del charge_complexes[c_id]['ligand_names']

        return charge_complexes

    def calculate_ligand_charges(self, max_iterations=None):
        charge_complexes = self.get_complex_dict_for_LCS()
        df_ligand_charges = get_charges_of_unique_ligands(all_complexes=charge_complexes, max_iterations=max_iterations)

        return df_ligand_charges

    def assign_charges_to_unique_ligands(self, max_charge_iterations: Union[int, None]):
        """
        Assigns charges to the unique ligands in the database. Currently, only the Linear Charge Solver method is used for this.
        """
        print('\nCHARGE CALCULATION:')
        df_ligand_charges = self.calculate_ligand_charges(max_iterations=max_charge_iterations)
        df_ligand_charges = df_ligand_charges.set_index('unique_name')
        df_ligand_charges = df_ligand_charges[
            [col for col in df_ligand_charges.columns if not col in self.LCS_needed_ligand_props]]
        self.df_unique_ligand_db = self.df_unique_ligand_db.join(df_ligand_charges)
        self.df_unique_ligand_db['is_confident'] = self.df_unique_ligand_db['is_confident'].fillna(False)
        self.n_pred_charges = self.df_unique_ligand_db['charge'].notna().sum()

        return

    def iterate_over_complexes(self):
        """
        Iterates over all complexes in the database and yields them one by one.
        """
        if self.store_database_in_memory:
            for c_id, c in self.complex_db.db.items():
                yield c_id, c
        else:
            with jsonlines.open(self.tmp_output_complexes_json, 'r') as reader:
                for line in reader:
                    c_id = line['key']
                    c = RCA_Complex.from_dict(dict_=line['value'])
                    yield c_id, c

        return

    def update_and_save_ligand(self, lig, writer=None):
        """
        Updates the ligand with unique ligand information and charges and saves it to disk.
        """
        share_properties = ['unique_name']
        collect_properties = {}

        # Update ligands with unique ligand information
        uname = self.ligand_to_unique_ligand[lig.name]
        ulig_props = series2namedtuple(self.df_unique_ligand_db.loc[uname])
        self.update_ligand_with_unique_ligand_information_inplace(
            lig=lig,
            ulig=ulig_props,
            share_properties=share_properties,
            collect_properties=collect_properties
        )
        update_ligand_with_charge_inplace(lig, charges=self.charge_dict)

        # Write ligand to disk as json
        if writer is not None:
            lig._append_to_file(key=lig.name, writer=writer)

        return

    def update_and_save_unique_ligand(self, lig, writer):
        ulig = self.get_unique_ligand_from_ligand(ligand=lig)
        uname = ulig.unique_name
        ulig._append_to_file(key=uname, writer=writer)
        if self.store_database_in_memory:
            self.unique_ligand_db[uname] = ulig

        return

    def update_and_save_complex(self, c, writer):
        # Update global props of complex with some useful information about unique ligands
        df_complex_ligands = self.df_full_ligand_db.loc[self.df_full_ligand_db['parent_complex_id'] == c.mol_id, ['unique_name', 'n_ligand_instances']]
        n_ligands = len(df_complex_ligands)
        c.global_props['n_ligands'] = n_ligands
        c.global_props['n_unique_ligands'] = len(set(df_complex_ligands['unique_name']))
        n_ligands_occurring_once = sum(df_complex_ligands['n_ligand_instances'] == 1)
        c.global_props['n_ligands_occurring_once'] = n_ligands_occurring_once
        c.global_props['frac_ligands_occurring_once'] = n_ligands_occurring_once / n_ligands

        # Write complex data to disk as json
        c._append_to_file(key=c.mol_id, writer=writer)

        return
    def save_databases_to_json(self):
        """
        Saves the databases to jsonlines files. At the same time, additional properties are calculated and stored. This function can read complexes from a temporary jsonlines file if the database is not stored in memory to reduce memory usage. This function does both the calculation of additional properties and the saving of the databases to jsonlines files so that each complex needs to be read only once.
        """

        # Prepare some data for the calculation of additional properties.
        self.all_most_common_n_H = self.df_unique_ligand_db.groupby(['heavy_atoms_graph_hash_with_metal'])[
            'n_hydrogens'].agg(lambda x: Counter(x).most_common(1)[0][0])
        self.all_charges_with_same_graph_hash = self.df_unique_ligand_db.groupby('graph_hash')['charge'].agg(
            lambda x: dict(Counter(x)))
        self.charge_dict = self.df_unique_ligand_db.to_dict(orient='index')

        if self.store_database_in_memory:
            self.unique_ligand_db = {}

        # Open jsonlines files for writing of unique ligands and complexes
        encoder = functools.partial(json.dumps, cls=NumpyEncoder)
        with jsonlines.open(self.output_complexes_json, mode='w', dumps=encoder) as complex_writer:
            with jsonlines.open(self.unique_ligands_json, mode='w', dumps=encoder) as ulig_writer:

                # Iterate once over all complexes
                desc = 'Writing databases to disk'
                for c_id, c in tqdm(self.iterate_over_complexes(), desc=desc, total=self.n_complexes):
                    for lig in c.ligands:
                        self.update_and_save_ligand(lig=lig, writer=None)
                        if lig.is_chosen_unique_ligand:
                            self.update_and_save_unique_ligand(lig=lig, writer=ulig_writer)
                    self.update_and_save_complex(c, writer=complex_writer)

                # Delete temporary json file
                self.tmp_output_complexes_json.unlink()

        if self.store_database_in_memory:
            self.unique_ligand_db = LigandDB(self.unique_ligand_db)

        return

    def run_ligand_extraction(self,
                              overwrite_atomic_properties: bool = True,
                              use_existing_input_json: bool = True,
                              max_charge_iterations: Union[int, None] = 10,
                              **kwargs
                              ):
        """
        Runs the entire ligand extraction process from reading in the .xzy files to optionally assigning charges.
        """
        start = datetime.now()

        # Read in complexes from xyz files and graphs
        self.ensure_input_complex_db_exists(overwrite_atomic_properties=overwrite_atomic_properties,
                                            use_existing_input_json=use_existing_input_json,
                                            **kwargs)

        # Extract ligands from complexes and save them in property 'ligands' of each complex.
        self.extract_ligands()

        # Assign charges to unique ligands
        self.assign_charges_to_unique_ligands(max_charge_iterations=max_charge_iterations)

        # Save complexes, ligands, unique ligands into a jsonlines file each
        self.save_databases_to_json()

        duration = get_duration_string(start)
        print(f'\nDuration of extraction: {duration}')
        print(f'Ligand database with charges (n={self.n_pred_charges}/{self.n_unique_ligands}) established successfully!')

        return


class RCA_Complex(BaseMolecule):

    def __init__(self,
                 mol: Atoms = None,
                 atomic_props: dict = None,
                 global_props: dict = None,
                 pymat_mol=None,
                 graph=None,
                 graph_creating_strategy="default",
                 has_ligands=True,
                 reindex_graph: bool = False,
                 other_props={},
                 **kwargs):

            super().__init__(
                            mol=mol,
                            atomic_props=atomic_props,
                            global_props=global_props,
                            pymat_mol=pymat_mol,
                            graph=graph,
                            graph_creating_strategy=graph_creating_strategy,
                            has_ligands=has_ligands,
                            reindex_graph=reindex_graph,
                            other_props=other_props,
                            **kwargs
                             )
            self._check_if_molecule_valid()
            self.shift_metal_to_origin()

            self.metal = identify_metal_in_ase_mol(self.atoms)
            self.metal_atomic_number = Element(self.metal).atomic_number
            self.metal_oxi_state = make_None_to_NaN(self.global_props['metal_oxi_state'])
            self.charge = make_None_to_NaN(self.global_props['charge'])
            self.metal_idx = self.atomic_props['atoms'].index(self.metal)
            self.metal_position = (self.atomic_props['x'][self.metal_idx], self.atomic_props['y'][self.metal_idx], self.atomic_props['z'][self.metal_idx])
            self.fully_connected = nx.is_connected(self.graph)

            self.donor_indices = sorted(self.graph.neighbors(self.metal_idx))
            # self.donor_elements = [self.atomic_props['atoms'][idx] for idx in self.donor_indices]
            self.n_donors = len(self.donor_indices)

            self.metal_position = (self.atomic_props['x'][self.metal_idx], self.atomic_props['y'][self.metal_idx], self.atomic_props['z'][self.metal_idx])
            self.fully_connected = nx.is_connected(self.graph)
            self.mol_id = self.global_props['CSD_code']

            # Set parameters for octahedral complexes
            self.is_octahedral = self.check_octahedral()
            mean_distance, sd_distance, max_angular_deviation = self.calculate_distortion_parameters()
            self.global_props['donors_mean_dist'] = mean_distance
            self.global_props['donors_sd_dist'] = sd_distance
            self.global_props['oct_max_ang_dev'] = max_angular_deviation

            self.add_additional_complex_information_to_global_props()

            return

    def has_metal_os(self):
        return ~np.isnan(self.metal_oxi_state)

    def calculate_distortion_parameters(self) -> tuple[float, float, float]:
        """
        Function to calculate distortion parameters of the complex
        :return: sd_distance (standard deviation of neighbour distances from the metal atom) and max_angular_deviation (maximum deviation from 90 or 180 degrees in bond angles)
        """
        # Calculate neighbour distances and positions
        metal_pos = np.array(self.metal_position)
        neighbours = self.donor_indices
        neighbour_positions = []
        neighbour_distances = []
        for neighbour_idx in neighbours:
            neighbour_pos = np.array([self.atomic_props['x'][neighbour_idx], self.atomic_props['y'][neighbour_idx], self.atomic_props['z'][neighbour_idx]])
            neighbour_positions.append(neighbour_pos)
            neighbour_distances.append(np.linalg.norm(neighbour_pos - metal_pos))

        # Calculate standard deviation of neighbour distances
        if len(neighbour_distances) > 0:
            mean_distance = float(np.mean(neighbour_distances))
            sd_distance = float(np.std(neighbour_distances))
        else:
            mean_distance = np.nan
            sd_distance = np.nan

        # Calculate maximum angular deviation for octahedral complexes
        if len(neighbours) == 6:
            max_angular_deviation = 0
            for i in range(len(neighbour_positions)):
                for j in range(i+1, len(neighbour_positions)):
                    # Calculate the vector from the metal to the neighbours
                    vector_i = neighbour_positions[i] - metal_pos
                    vector_j = neighbour_positions[j] - metal_pos

                    # Calculate the angle between these vectors
                    cos_angle = np.dot(vector_i, vector_j) / (np.linalg.norm(vector_i) * np.linalg.norm(vector_j))
                    cos_angle = np.clip(cos_angle, -1, 1)  # Ensure cos_angle is within the valid range [-1, 1]
                    angle = np.arccos(cos_angle) * 180 / np.pi  # Convert to degrees

                    # Check the angular deviation from 90 or 180 degrees
                    angular_deviation = min(abs(angle - 90), abs(angle - 180))
                    max_angular_deviation = max(max_angular_deviation, angular_deviation)
        else:
            max_angular_deviation = np.nan

        return mean_distance, sd_distance, max_angular_deviation

    def check_octahedral(self, angular_threshold=40) -> bool:
        """
        Function to check if the complex is octahedral or not
        The conditions for a complex to be octahedral are:
            - The complex must have 6 donor atoms
            - The maximum deviation from 90 or 180 degrees in bond angles should be less than `angular_threshold`
        :return: True if the complex is octahedral, False otherwise
        """
        if len(self.donor_indices) != 6:
            return False

        _, _, max_angular_deviation = self.calculate_distortion_parameters()
        is_oct = bool(max_angular_deviation > angular_threshold)    # bool() for json serialisation

        return is_oct

    def count_ligands_with_stoichiometry(self, atoms: list, only_connected=False):
        if isinstance(atoms, str):
            atoms = [atoms]
        atoms = sorted(atoms)

        n = 0
        for lig in self.ligands:
            if not only_connected or lig.was_connected_to_metal:
                lig_atoms = sorted(lig.atomic_props['atoms'])
                if atoms == lig_atoms:
                    n += 1

        return n

    def count_n_unconnected_ligands(self, max_n_atoms: np.inf) -> int:
        """
        Counts the number of unconnected ligands with a maximum number of atoms of `max_n_atoms`.
        :param max_n_atoms: Maximum number of atoms of ligands to count
        :return: Number of ligands with this maximum number of atoms
        """
        n = sum([not lig.was_connected_to_metal and len(lig.atomic_props['atoms']) <= max_n_atoms for lig in self.ligands])
        return n


    def count_atoms_in_ligands(self, atoms: Union[list, str], only_if_connected_to_metal: bool=False, per_ligand:bool=False) -> Union[list, int]:
        """
        Returns the number of n_ligand_instances of the specified elements in the complex.
        :param atoms: list of atoms to count the total n_ligand_instances
        :param only_if_connected_to_metal: If True, ignore unconnected ligands
        :param per_ligand:If True, returns a list of occupancies for each ligand. Otherwise returns the total for all complexes.
        :return: Number of n_ligand_instances of specified elements, either as list per ligand or the total as integer
        """
        if isinstance(atoms, str):
            atoms = [atoms]

        n = []
        for lig in self.ligands:
            if only_if_connected_to_metal and lig.was_connected_to_metal:
                n.append(sum([el in atoms for el in lig.atomic_props['atoms']]))

        if not per_ligand:
            n = sum(n)

        return n

    def get_only_complex_graph_connected_to_metal(self, atom_label: object = 'node_label') -> nx.Graph:
        """
        Returns the graph of only the metal complex without unconnected ligands.
        """
        complex_graph = get_only_complex_graph_connected_to_metal(graph=self.graph, atom_label=atom_label, metal=self.metal)
        return complex_graph

    def has_fragment(self, frag: Union[str, list]) -> bool:
        """
        Checks whether the complex graph has an unconnected fragment with the elements specified in frag.
        :param frag: Either a string specifying a single atom or a list of strings specifying several atoms. E.g. frag='Cl' checks if any unconnected chloride exists. frag='[H, H, O]' checks whether any unconnected water exists. The order in the list doesn't matter but the number of n_ligand_instances of elements does matter.
        """
        if isinstance(frag, str):
            frag = [frag]

        fragments = [sorted(comp) for comp in nx.connected_components(self.graph)]
        el_fragments = [[self.atomic_props['atoms'][i] for i in fragment] for fragment in fragments]

        has_unconnected_fragment = False
        frag = sorted(frag)
        for fragment in el_fragments:
            fragment = sorted(fragment)
            if fragment == frag:
                has_unconnected_fragment = True
                break

        return has_unconnected_fragment

    def count_ligands_containing_only(self,
                                      atoms: Union[str, list],
                                      denticity_range=None,
                                      n_atoms_range=None,
                                      except_elements=None
                                      ) -> int:
        """
        Count how many ligands contain only atoms specified in `atoms`, except for elements in `except elements`.
        :param atoms: atoms to check for.
        :param denticity_range: if specified, count ligands only if it has a n_donors in this range (inclusive)
        :param n_atoms_range: if specified, count ligands only if the number of atoms which are not excluded by `except_element` is in this range (inclusive)
        :param except_elements: Ignore these elements in the ligands when checking the presence of the atoms in the ligand and when checking n_atoms_range.
        :return: Integer
        """
        if denticity_range is None:
            denticity_range = []
        if n_atoms_range is None:
            n_atoms_range = []
        if except_elements is None:
            except_elements = []
        n = 0
        for lig in self.ligands:
            correct_denticity = is_between(lig.n_donors, denticity_range)
            correct_n_atoms = is_between(len([a for a in lig.atomic_props['atoms'] if not a in except_elements]), n_atoms_range)
            if correct_denticity and correct_n_atoms and lig._contains_only(atoms, except_elements=except_elements):
                n += 1

        return n


    def _count_coordinating_atoms_with_distance_to_metal_greater_than(self, distance, element=None, max_n_atoms=np.inf) -> int:
        """
        Returns the number of coordinating atoms of type `element` with a distance to the metal greater than `distance`.
        :param max_n_atoms: Include only ligands with a number of atoms up to this number.
        """
        n = 0
        for lig in self.ligands:
            for idx, el in enumerate(lig.donor_elements):
                correct_atom = el == element or element is None
                if correct_atom and len(lig.atomic_props['atoms']) <= max_n_atoms:
                    if lig._get_atomic_distance_to_original_metal(mode='coordinating')[idx] > distance:
                        n += 1

        return n

    def has_consistent_stoichiometry_with_CSD(self, print_different_elements=False) -> bool:
        try:
            csd_stoichiometry = self.global_props['CSD_stoichiometry']
            pattern = r'[A-Z][a-z]?\d*'
            csd_stoichiometry = ' '.join(re.findall(pattern, csd_stoichiometry))
            csd_comp = Composition(csd_stoichiometry)
        except KeyError:
            warnings.warn('Global property `CSD_stoichiometry` not found. Skip check for consistent stoichiometry with xyz.')
            return True
        except:
            warnings.warn(f'CSD_stoichiometry could not be parsed for complex {self.mol_id}. Skip check for consistent stoichiometry with xyz.')
            return True

        xyz_comp = Composition(self.stoichiometry)
        consistent = csd_comp.almost_equals(xyz_comp)

        if not consistent and print_different_elements:
            xyz_elements = [Element(el).symbol for el in xyz_comp.elements]
            csd_elements = [Element(el).symbol for el in csd_comp.elements]
            different_elements = set(xyz_elements).symmetric_difference(csd_elements)
            if different_elements:
                print(f'Differing elements in stoichiometries of xyz and CSD in complex {self.mol_id}: {list(different_elements)}')
            diff_dict = {el: int(xyz_comp.elements[el]-csd_comp.elements[el]) for el in set(csd_elements+xyz_elements)}
            diff_dict = {key: val for key, val in diff_dict.items() if val != 0}
            print(f'Different elements counts for complex {self.mol_id}: {diff_dict}')

        return consistent

    def has_consistent_stoichiometry_with_smiles(self, smiles: str, ignore_element_count: bool=False, print_warnings: bool=True, only_complex: bool=False) -> bool:
        try:
            mol_graph = pysmiles.read_smiles(smiles, explicit_hydrogen=True, zero_order_bonds=False)

            if only_complex:
                _, frag_atoms = get_graph_fragments(graph=mol_graph, node_label='element')
                atoms = [a for a in frag_atoms if self.metal in a][0]
            else:
                atoms, _ = get_sorted_atoms_and_indices_from_graph(mol_graph, atom_label='element')

            comp = Composition(''.join(atoms))
        except KeyError:
            if print_warnings:
                warnings.warn('Global property `smiles` not found. Skip check for consistent stoichiometry with xyz.')
            return True
        except:
            if print_warnings:
                warnings.warn(f'Smiles could not be parsed for complex {self.mol_id}. Skip check for consistent stoichiometry with xyz.')
            return True

        if only_complex:
            _, frag_atoms = self.get_graph_fragments()
            atoms = [a for a in frag_atoms if self.metal in a][0]
            xyz_comp = Composition(''.join(atoms))
        else:
            xyz_comp = Composition(self.stoichiometry)

        if ignore_element_count:
            # Check only if the same elements occur, not if the elements have the same count.
            xyz_elements = [Element(el).symbol for el in xyz_comp.elements]
            csd_elements = [Element(el).symbol for el in comp.elements]
            different_elements = set(xyz_elements).symmetric_difference(csd_elements)
            consistent = len(different_elements) == 0
        else:
            consistent = comp.almost_equals(xyz_comp)

        return consistent


    def complex_is_biggest_fragment(self, allow_complexes_greater_than: int = np.nan) -> bool:
        """
        Checks whether the complex (fragment with the transition metal) is the fragment with the most atoms. This is useful to check whether the transition metal might not belong to an actual complex but to a counterion.
        :param allow_complexes_greater_than: Always return true for complexes with a higher number of atoms than this parameter
        """
        _, el_fragments = self.get_graph_fragments()

        max_n_other_atoms = 0
        n_complex_atoms = 0
        for atoms in el_fragments:
            is_complex = any([Element(atom).is_metal for atom in atoms])
            if is_complex:
                assert n_complex_atoms == 0, f'There seem to be complexes with more than one transition metal in complex {self.mol_id}.'
                n_complex_atoms = len(atoms)
            else:
                if len(atoms) > max_n_other_atoms:
                    max_n_other_atoms = len(atoms)
        assert n_complex_atoms != 0, f'There seem to be no transition metals in complex {self.mol_id}.'

        if max_n_other_atoms > n_complex_atoms and n_complex_atoms <= allow_complexes_greater_than:
            is_biggest_fragment = False
            # print(f'N complex atoms in {self.mol_id}: {n_complex_atoms}')
        else:
            is_biggest_fragment = True

        return is_biggest_fragment

    def add_additional_complex_information_to_global_props(self):
        info = {}
        if 'partial_charge' in self.atomic_props:
            metal_q = self.atomic_props['partial_charge'][self.atomic_props['atoms'].index(self.metal)]
            info['metal_partial_charge'] = metal_q

        update_dict_with_warning_inplace(
                                            dict_to_update=self.global_props,
                                            dict_with_information=info
                                        )

        return info

    def shift_metal_to_origin(self):
        """
        Actually, this method is outdated not required.
        However, the idea was to shift the molecule to the origin, so that the metal is right in the origin
        which is the first part - the shift

        and afterwards, rotate the the metal, so that one of the functional atoms is algined with the x-axis
        """
        #
        # get shift vector
        metal_symb = identify_metal_in_ase_mol(self.atoms)
        metal_index = self.atomic_props["atoms"].index(metal_symb)

        # do the shifting
        for key in ["x", "y", "z"]:
            self.atomic_props[key] = [value - self.atomic_props[key][metal_index] for value in self.atomic_props[key]]

        """
        # Now we rotate. First we identify the element we want to rotate on
        metal_node = find_node_in_graph_by_label(G=self.graph, label_to_find=metal_symb, expected_hits=1)
        rotation_element_index = list(self.graph.neighbors(metal_node))[0]
        # Next, we obtain the rotation matrix
        rotation_element_vector = np.array([self.atomic_props[key][rotation_element_index] for key in ["x", "y", "z"]])
            # now we can idenftify the desired vector, on which we'd like to rotate the vector
        desired_rotation = np.array([np.linalg.norm(rotation_element_vector), 0, 0])
        # and we can thus find the rotation matrix
        A = R.align_vectors(rotation_element_vector.reshape(-1, 1).T, desired_rotation.reshape(-1, 1).T)[0].as_matrix()
        for index, _ in enumerate(self.atomic_props["x"]):
            location_vec_for_index = np.array([self.atomic_props[key][index] for key in ["x", "y", "z"]])
            # now we rotate
            v_new = location_vec_for_index.dot(A)
            # and we modify the atomic properties
            for i, key in enumerate(["x", "y", "z"]):
                self.atomic_props[key][index] = v_new.tolist()[i] if abs(v_new.tolist()[i]) > 0.001 else 0
        """
        return

    def get_output_info(self) -> dict:

        info = {
            'Complex ID': self.mol_id,
            'Stoichiometry': self.stoichiometry,
            'Metal': self.metal,
            'Metal OS': self.metal_oxi_state,
            'Charge': self.charge,
            'Number of Atoms': self.n_atoms,
            'Fully Connected': self.fully_connected,
            'Coordinating Atoms': self.donor_elements,
            'Number of Donors': self.n_donors,
            'Octahedral': self.is_octahedral,
            'Haptic': any(lig.is_haptic for lig in self.ligands),
            'Beta-Hydrogen': any(lig.n_beta_hydrogens > 0 for lig in self.ligands),
            'Ligand Unique Names': [lig.unique_name for lig in self.ligands],
            **self.global_props,
            }

        return info

    def de_assemble(self, inherit_global_properties: list = ['CSD_code']):
        """
        now only graph based, makeslife waaay easier

        All that is based on the following assumption:
        the nodes are denoted as integers and we shall assume that these integers
        correspond to their index in the atomic properties, i.e. the position
        so, for example for node 1, the x-coordinate can be found by
        self.atomic_props["x"][1] or more general self.atomic_props["x"][node]
        Ich glaube das haelt auf jeden Fall fuer alle selber erzeugten Graphen
        den Dummen Testparamater muss ich mitschleppen, um das gegen meine alte, ultra behinderte Grapherstellung testen
        zu koennen
        """
        if not hasattr(self, "ligands"):
            self.ligands = []

        inherit_global_properties = self._check_input_inherit_global_properties(inherit_global_properties)

        atoms, idc = get_sorted_atoms_and_indices_from_graph(self.graph)
        if 'atoms' in self.atomic_props:
            # if not atoms == self.atomic_props['atoms']:
            #     breakpoint()

            assert atoms == self.atomic_props['atoms'], 'Order of atoms in graph and in atomic_props doesn\'t match.'

        # first we gather some information around the metal in the initial graph
        graph = deepcopy(self.graph)
        metal_in_complex = identify_metal_in_ase_mol(self.atoms)
        metal_node = find_node_in_graph_by_label(G=graph, label_to_find=metal_in_complex, expected_hits=1)
        metal_neighbors = list(graph.neighbors(metal_node))  # all neighbor nodes of the metal
        metal_neighbor_elements = [graph.nodes[i]['node_label'] for i in metal_neighbors]

        for i, el in zip(idc, atoms):
            #graph.nodes[i]['orig_idx'] = i
            #graph.nodes[i]['metal_neighbor'] = i in metal_neighbors
            assert graph.nodes[i]['node_label'] == el, 'atom and index don\'t match.'

        # next we create the ripped graph
        ripped_graph = deepcopy(graph)
        ripped_graph.remove_node(metal_node)

        conn_components = [sorted(comp) for comp in
                           nx.connected_components(ripped_graph)]  # sorting of components very important
        conn_components = [comp for comp in
                           sorted(conn_components,
                                  key=str)]  # important: sort by string of list, makes order of ligands unique

        for component in conn_components:
            # if this set is empty, the ligand has no connection to the metal
            functional_atom_indices = sorted(list(set(component).intersection(set(metal_neighbors))))
            assert max(component) <= len(
                ripped_graph), 'Indices dont make sense. Most likely this is an implementation error where due to the deletion of the metal atom the indices of original and ripped graph dont match.'

            denticity = len(functional_atom_indices)
            if denticity == 0:
                # because in the assembly 0 is the placeholder for the reactant, whereas -1 means this is just a remainder in the .xyz, not
                # connected to the metal at all
                denticity = -1

            # with that it is insanely easy to determine the atomic properties of the ligand
            assert component == sorted(
                component), 'The list of ligand indices is not sorted, but that is assumed in many parts of this project.'
            ligand_atomic_props = {key_: [el for i, el in enumerate(item_) if i in component] for key_, item_ in
                                   self.atomic_props.items()}
            ligand_atomic_props['original_complex_indices'] = component

            ligand_graph = deepcopy(ripped_graph.subgraph(component))
            atoms_lig, idc_lig = get_sorted_atoms_and_indices_from_graph(ligand_graph)

            # problem: the functional_atom_indices are the indices in the full original metal, rather than the ligand only
            # so we have to convert them to the index in the ligand_atomic_props
            local_indices = [component.index(ind) for ind in functional_atom_indices]
            local_elements = [ligand_atomic_props['atoms'][i] for i in local_indices]

            # Doublechecking
            if local_indices != []:  # otherwise error for unconnected ligands
                assert max(local_indices) < len(
                    ligand_graph), 'local_indices make no sense, an index is greater than the number of elements.'
            assert all(
                el in metal_neighbor_elements for el in local_elements), 'Inconsistent elements of the metal neighbors.'
            assert local_indices == sorted(local_indices), 'local_indices is not sorted but should be.'
            assert atoms_lig == ligand_atomic_props['atoms'], 'elements of graph and atomic_props not consistent'
            assert [atoms[i] for i in component] == ligand_atomic_props[
                'atoms'], 'ligand atoms not consistent with original complex atoms.'

            ligand_name, csd = self._ligand_naming(denticity, self.ligands)

            ligand_global_props = {prop: self.global_props[prop] for prop in inherit_global_properties}
            new_lig = Ligand(denticity=denticity,
                             donor_idc=local_indices,
                             atomic_props=ligand_atomic_props,
                             name=ligand_name,
                             graph=ligand_graph,
                             global_props=ligand_global_props,
                             original_metal=Element(metal_in_complex).atomic_number,
                             original_metal_position=self.metal_position,
                             original_metal_os=self.metal_oxi_state
                             )

            self.ligands.append(new_lig)

        if not self.ligands:
            print(f'WARNING: Complex {self.global_props["CSD_code"]} has no ligands extracted.')

    def get_smiles(self, only_core_complex: bool=False) -> str:
        full_smiles = super().get_smiles()
        smiles = full_smiles
        if only_core_complex:
            smiles = [sm for sm in smiles.split('.') if self.metal in sm]
            assert len(smiles) == 1, 'There should be exactly one SMILES string containing the metal.'
            smiles = smiles[0]

        # Assert
        smiles_graph = pysmiles.read_smiles(smiles, explicit_hydrogen=True)
        smiles_hash = get_graph_hash(smiles_graph, node_attr='element')
        core_graph = self.get_only_complex_graph_connected_to_metal()
        core_graph_hash = get_graph_hash(core_graph)
        if not smiles_hash == core_graph_hash and len(core_graph.nodes) < 9999999999:
            view_graph(core_graph, save_path='/Users/timosommer/Downloads/core_graph.png')
            view_graph(smiles_graph, save_path='/Users/timosommer/Downloads/smiles_graph.png', node_label='element')

            a = 1

        # assert smiles_hash == core_graph_hash, 'SMILES and graph hash do not match.'

        return smiles

    def generate_descriptors_of_complex_graph(self, only_core_complex: bool=True, xtb=False):
        """
        Generates descriptors for the complex and its ligands.
        """
        from dev.src11_machine_learning.utils.utilities_ML import get_element_descriptors
        from dev.src11_machine_learning.dataset_preparation.descriptors import RAC
        descriptors = {}

        # Complex descriptors
        descriptors['n_ligands'] = sum(1 for lig in self.ligands if lig.n_donors > 0 or not only_core_complex)

        # Metal center descriptors
        metal_descriptors = get_element_descriptors(self.metal)
        metal_descriptors = {f'metal_{key}': val for key, val in metal_descriptors.items()}
        descriptors.update(metal_descriptors)

        # Coordinating_atom descriptors
        coords_descriptors = [get_element_descriptors(el) for el in self.donor_elements]
        coords_descriptors = {key: np.mean([d[key] for d in coords_descriptors]) for key in coords_descriptors[0].keys()}
        coords_descriptors = {f'coord_{key}': val for key, val in coords_descriptors.items()}
        descriptors.update(coords_descriptors)

        # Graph descriptors
        # RDKit descriptors
        # smiles = self.get_smiles(only_core_complex=only_core_complex)
        # rdkit_descriptors = RDKit_2D([smiles]).compute_2Drdkit().to_dict(orient='records')[0]
        # descriptors.update(rdkit_descriptors)
        # Own RAC descriptors
        graph = self.get_only_complex_graph_connected_to_metal() if only_core_complex else self.graph
        graph_descriptors, labels = RAC().molecule_autocorrelation(mol=graph, return_labels=True)
        graph_descriptors = {f'own_graph_{key}': val for key, val in zip(labels, graph_descriptors)}
        descriptors.update(graph_descriptors)


        # XTB descriptors of ligands
        if xtb:
            xtb_descriptors = []
            for lig in self.ligands:
                if lig.n_donors > 0 or not only_core_complex:
                    xtb_descriptors.append(lig._get_xtb_descriptors())


        return descriptors

    def to_dict(self, include_graph: bool=True) -> dict:
        d = {}

        # Manually initialize special fields
        if include_graph:
            d['graph'] = graph_to_dict_with_node_labels(self.graph)
        d['ligands'] = [lig.to_dict() for lig in self.ligands]

        # Do not output write these fields to the output dictionary, mostly because they are not json serializable
        do_not_output_automatically = ['ligands', 'node_label', 'atoms', 'smiles', 'rdkit_mol', 'graph', 'coordinates', 'hash', 'csd_code', 'graph_with_metal']

        for prop, val in vars(self).items():
            if not prop in do_not_output_automatically:
                d[prop] = val

        return d

    @classmethod
    def from_dict(cls, dict_: dict, **kwargs):
        """
        Reads the ligand from a provided dictionary.
        """
        necessary_props = ["atomic_props", "global_props", "graph_dict"]
        assert set(necessary_props).issubset(set(dict_.keys())), f'Any of the necessary keys {necessary_props} is not present.'

        if 'ligands' in dict_:
            dict_['ligands'] = [Ligand.from_dict(lig) for lig in dict_['ligands']]
            # This sounds stupid but this is because otherwise the RCA_molecule class sets up self.ligands = [] which then collides when the actual ligands which are read in here are added because for safety the code checks that it doesn't overwrite anything.
            has_ligands = False
        else:
            has_ligands = True

        if 'total_q' in dict_:
            dict_['charge'] = dict_['total_q']
            del dict_['total_q']

        other_props = {key: val for key, val in dict_.items() if not key in necessary_props}

        # Optionally add graph if it is present in the dictionary
        if 'graph_dict' in dict_ and not (dict_['graph_dict'] is None):
            graph = graph_from_graph_dict(dict_['graph_dict'])
        else:
            graph = None

        return cls(
            atomic_props=dict_["atomic_props"],
            global_props=dict_["global_props"],
            graph=graph,
            has_ligands=has_ligands,
            other_props=other_props,
            **kwargs
        )


class ComplexDB(BaseDB):
    type = 'Complex'
    def __init__(self, db):
        super().__init__(db=db)

    def check_db_equal(self, db: str):
        db = ComplexDB.from_json(json_=db, type_='Complex')
        return self == db

    @classmethod
    def load_from_json(cls, path: Union[str, Path], n_max: int=None, show_progress: bool=True):
        """
        Load a JSON or JSON Lines file.
        :param path: Path to the JSON or JSON Lines file
        :return: A ComplexDB object
        """
        db = load_complex_db(path=path, n_max=n_max, show_progress=show_progress, molecule='class')
        return cls(db)


    def to_dataframe(self) -> pd.DataFrame:
        """
        Creates a csv file of the database with the most important information.
        """
        data = []
        for name, complex in tqdm(self.db.items(), desc='Create csv file'):
            mol_data = complex.get_output_info()
            mol_data = {key: value for key, value in mol_data.items() if not key.startswith('Unnamed:')}
            data.append(mol_data)
        df = pd.DataFrame(data)

        return df

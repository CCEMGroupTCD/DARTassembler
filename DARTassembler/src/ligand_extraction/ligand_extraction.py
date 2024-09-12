"""
Class for extracting ligands from a database of complexes.
"""
import functools

import jsonlines
import pandas as pd
from copy import deepcopy
from DARTassembler.src.ligand_extraction.DataBase import LigandDB, ComplexDB
from pathlib import Path
from tqdm import tqdm
import json
from typing import Union
import numpy as np
from datetime import datetime
from collections import Counter

from DARTassembler.src.ligand_extraction.DataLoader import DataLoader
from DARTassembler.src.ligand_extraction.Molecule import RCA_Complex, RCA_Ligand
from DARTassembler.src.ligand_extraction.io_custom import NumpyEncoder, iterate_over_json, get_n_entries_of_json_db
from DARTassembler.src.ligand_extraction.utilities_Molecule import unknown_rdkit_bond_orders
from DARTassembler.src.ligand_extraction.utilities_extraction import get_charges_of_unique_ligands, update_ligand_with_charge_inplace
from DARTassembler.src.ligand_extraction.utilities import get_duration_string, series2namedtuple
from DARTassembler.src.constants.testing import CHARGE_BENCHMARKED_COMPLEXES
from DARTassembler.src.constants.constants import odd_n_electrons_warning, \
    similar_molecule_with_diff_n_hydrogens_warning
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

        input_complex_db.to_json(path=self.input_complexes_json, json_lines=True)

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

        if c.global_props['smiles'] is None and not c.has_bond_type(unknown_rdkit_bond_orders):
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

        if c.count_atoms_with_n_bonds(element='C', n_bonds=1) > 0:
            self.excluded_complex_ids['C atom with only 1 bond'].append(c_id)
            return False

        min_dist, _, _ = c.get_atomic_distances_between_atoms()
        if min_dist < 0.5:
            self.excluded_complex_ids['Atoms closer than 0.5A'].append(c_id)
            return False

        min_dist, _, _ = c.get_atomic_distances_between_atoms(skip_elements='H')
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

        if comp.count_coordinating_atoms_with_distance_to_metal_greater_than(distance=1.9, element='O', max_n_atoms=1) > 0:
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
        @param props: list of properties. Must be the exact name of a property of the RCA_Ligand class.
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
        important_ligand_props = ['name', 'stoichiometry', 'n_protons', 'graph_hash', 'denticity', 'has_good_bond_orders',
                                    self.unique_ligand_id, 'original_metal_symbol', 'original_metal_os',
                                    'was_connected_to_metal', 'original_complex_id', 'n_hydrogens', 'original_complex_id',
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
                    comp = RCA_Complex.read_from_mol_dict(dict_=comp_dict)

                    # Normalize the complex
                    comp.remove_node_features_from_molecular_graphs_inplace()
                    comp.normalize_multigraph_into_graph_inplace()

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
                                comp.append_to_file(key=csd_code, writer=complex_writer)    # Write to jsonlines file
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
        # Add `occurrences` column to full ligand df
        self.df_full_ligand_db = pd.merge(self.df_full_ligand_db, self.df_unique_ligand_db['occurrences'], left_on='unique_name', right_index=True, how='left')

        self.print_excluded_complexes()

        return

    # @profile
    def get_unique_ligand_df(self) -> pd.DataFrame:
        """
        Returns a dataframe with all unique ligands from the full ligand db.
        """
        self.grouped_ligands = self.group_same_ligands()
        df_same_graph_hash = self.df_full_ligand_db.loc[:, ['graph_hash', 'denticity', 'graph_hash_with_metal']].groupby('graph_hash').agg(lambda x: x.unique().tolist())
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
            denticities = df_same_graph_hash.loc[graph_hash, 'denticity']
            metals = self.df_full_ligand_db.loc[same_ligand_names, 'original_metal_symbol']
            count_metals = metals.value_counts().sort_values(ascending=False).to_dict()
            n_graph_hashes = len(df_same_graph_hash.loc[graph_hash, 'graph_hash_with_metal'])
            assert not 0 in denticities, 'The denticity for unconnected ligands is assumed to be -1 but here there appears a 0.'
            has_unconnected_ligands = -1 in denticities
            unique_ligand_infos = {
                                    'occurrences': len(same_ligand_names),
                                    'same_graph_denticities': denticities,
                                    'count_metals': count_metals,
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
        self.full_ligand_db.to_json(self.full_ligands_json)

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
            same_ligand_props = self.df_full_ligand_db['has_good_bond_orders'].loc[same_ligands]
            for lig_name, has_good_bond_orders in same_ligand_props.items():
                if has_good_bond_orders:
                    name = lig_name
                    break
        else:
            raise ValueError(
                f'Unknown strategy `{strategy}` to choose the unique ligand representative from all same ligands.')

        return name

    def get_unique_ligand_from_ligand(self, ligand: RCA_Ligand) -> RCA_Ligand:
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
        original_complex_ids = df_same_ligands['original_complex_id'].tolist()
        identical_ligand_info['original_complex_charge'] = self.df_complex_db.loc[original_complex_ids, 'charge'].tolist()
        identical_ligand_info['original_complex_id'] = df_same_ligands['original_complex_id'].tolist()
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

        ulig.n_electrons = ulig.n_protons - ulig.pred_charge
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
                [lig.unique_ligand_information['occurrences'] == 1 for lig in c.ligands])
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
        self.n_pred_charges = self.df_unique_ligand_db['pred_charge'].notna().sum()

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
                    c = RCA_Complex.read_from_mol_dict(dict_=line['value'])
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
            lig.append_to_file(key=lig.name, writer=writer)

        return

    def update_and_save_unique_ligand(self, lig, writer):
        ulig = self.get_unique_ligand_from_ligand(ligand=lig)
        uname = ulig.unique_name
        ulig.append_to_file(key=uname, writer=writer)
        if self.store_database_in_memory:
            self.unique_ligand_db[uname] = ulig

        return

    def update_and_save_complex(self, c, writer):
        # Update global props of complex with some useful information about unique ligands
        df_complex_ligands = self.df_full_ligand_db.loc[self.df_full_ligand_db['original_complex_id'] == c.mol_id, ['unique_name', 'occurrences']]
        n_ligands = len(df_complex_ligands)
        c.global_props['n_ligands'] = n_ligands
        c.global_props['n_unique_ligands'] = len(set(df_complex_ligands['unique_name']))
        n_ligands_occurring_once = sum(df_complex_ligands['occurrences'] == 1)
        c.global_props['n_ligands_occurring_once'] = n_ligands_occurring_once
        c.global_props['frac_ligands_occurring_once'] = n_ligands_occurring_once / n_ligands

        # Write complex data to disk as json
        c.append_to_file(key=c.mol_id, writer=writer)

        return
    def save_databases_to_json(self):
        """
        Saves the databases to jsonlines files. At the same time, additional properties are calculated and stored. This function can read complexes from a temporary jsonlines file if the database is not stored in memory to reduce memory usage. This function does both the calculation of additional properties and the saving of the databases to jsonlines files so that each complex needs to be read only once.
        """

        # Prepare some data for the calculation of additional properties.
        self.all_most_common_n_H = self.df_unique_ligand_db.groupby(['heavy_atoms_graph_hash_with_metal'])[
            'n_hydrogens'].agg(lambda x: Counter(x).most_common(1)[0][0])
        self.all_charges_with_same_graph_hash = self.df_unique_ligand_db.groupby('graph_hash')['pred_charge'].agg(
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

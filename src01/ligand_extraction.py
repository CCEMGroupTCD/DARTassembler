"""
Class for extracting ligands from a database of complexes.
"""
import warnings

import pandas as pd
from copy import deepcopy
from src01.DataBase import LigandDB, ComplexDB
from pathlib import Path
from tqdm import tqdm
import json
from typing import Union
import networkx as nx
import numpy as np
from datetime import datetime
from collections import Counter

from src01.DataLoader import DataLoader
from src01.Molecule import RCA_Complex
from src01.io_custom import load_unique_ligand_db, load_complex_db, load_full_ligand_db, save_unique_ligand_db, \
    save_full_ligand_db, save_complex_db
from src01.utilities_Molecule import unknown_rdkit_bond_orders
from src01.utilities_extraction import unique_ligands_from_Ligand_batch_json_files, update_complex_db_with_ligands, \
    get_charges_of_unique_ligands, update_databases_with_charges, update_ligand_with_charge_inplace
from src01.utilities import sort_dict_recursively_inplace, update_dict_with_warning_inplace, unroll_dict_into_columns, \
    get_duration_string
from constants.testing import CHARGE_BENCHMARKED_COMPLEXES
from constants.constants import odd_n_electrons_warning, unconfident_charge_warning, similar_molecule_with_diff_n_hydrogens_warning
from collections import defaultdict


class LigandExtraction:

    def __init__(self, database_path: str,
                 data_store_path: str,
                 exclude_not_fully_connected_complexes: bool = True,
                 testing: Union[bool, int] = False,
                 graph_strat: str = "default",
                 save_intermediate_databases: bool = False,
                 exclude_charged_complexes: bool = False,
                 only_complexes_with_os: bool = False,
                 only_real_transition_metals: bool = False,
                 unique_ligand_id: str = 'graph_hash_with_metal'
                 ):

        self.ligand_to_unique_ligand = None
        self.unique_ligand_info_props = None
        self.grouped_ligands = None
        self.database_path = ""
        self.data_store_path = ""
        self.exclude_not_fully_connected_complexes = None
        self.testing = None
        self.complex_db = None
        self.basic_ligand_infos = None
        self.test_complexes = None
        self.exclude_charged_complexes = exclude_charged_complexes
        self.only_complexes_with_os = only_complexes_with_os
        self.only_real_transition_metals = only_real_transition_metals
        self.test_complexes = CHARGE_BENCHMARKED_COMPLEXES if not testing == False else []
        self.unique_ligand_id = unique_ligand_id

        self.excluded_complex_ids = defaultdict(list)
        self.save_intermediate_databases = save_intermediate_databases
        self.graph_strat = graph_strat

        self.check_and_set_init_input(
            database_path=database_path,
            data_store_path=data_store_path,
            exclude_not_fully_connected_complexes=exclude_not_fully_connected_complexes,
            testing=testing,
        )

        self.input_complexes_json = Path(self.data_store_path, 'tmQMG.json')
        self.output_complexes_json = Path(self.data_store_path, 'complex_db.json')
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
        Establish and safe the Database (in our case tmQM) as json for simple loading.
        """
        db_dict = DataLoader(database_path_=self.database_path, overwrite=overwrite_atomic_properties).data_for_molDB

        if self.test_complexes:
            db_dict = self.reorder_input_complexes(db_dict=db_dict, first_complexes=self.test_complexes)

        if isinstance(self.testing, list):
            input_complex_db = ComplexDB.from_json(
                json_=db_dict,
                type_="Complex",
                identifier_list=self.testing,
                graph_strategy=self.graph_strat,
                **kwargs
            )
        else:
            input_complex_db = ComplexDB.from_json(
                                                json_=db_dict,
                                                type_="Complex",
                                                max_number=self.testing,
                                                graph_strategy=self.graph_strat,
                                                **kwargs
                                                )
        input_complex_db.to_json(path=self.input_complexes_json)

        return

    def delete_filtered_complexes_from_db(self):
        for c_ids in self.excluded_complex_ids.values():
            for c_id in c_ids:
                try:
                    del self.complex_db.db[c_id]
                except KeyError:
                    pass

    def prefilter_input_complex_db(self):
        """
        Filters input complexes in `self.complex_db` by multiple criteria without needing information about ligands.
        """
        for c_id, c in tqdm(self.complex_db.db.items(), desc='Filter input complexes'):
            if c.n_donors == 0:
                self.excluded_complex_ids['Metal ion'].append(c_id)

            if c.has_fragment(frag='O'):
                self.excluded_complex_ids['Has unconnected O'].append(c_id)

            if c.has_fragment(frag='H'):
                self.excluded_complex_ids['Has unconnected H'].append(c_id)

            if c.has_fragment(frag=['O', 'O']):
                self.excluded_complex_ids['Has unconnected O2'].append(c_id)

            if c.has_fragment(frag=['H', 'O']):
                self.excluded_complex_ids['Has unconnected OH'].append(c_id)

            if not c.global_props['is 3d']:
                self.excluded_complex_ids['Is not 3D'].append(c_id)

            if not c.has_consistent_stoichiometry_with_CSD():
                self.excluded_complex_ids['Inconsistent CSD stoichiometry'].append(c_id)

            # if not c.has_consistent_stoichiometry_with_smiles(smiles=c.global_props['smiles'], ignore_element_count=True, print_warnings=False):
            #     self.excluded_complex_ids['Inconsistent smiles elements'].append(c_id)

            if c.global_props['smiles'] is None and not c.has_bond_type(unknown_rdkit_bond_orders):
                self.excluded_complex_ids['No smiles without bad bonds'].append(c_id)

            # if not c.complex_is_biggest_fragment(allow_complexes_greater_than=10):
            #     self.excluded_complex_ids['Complex is counter ion'].append(c_id)

            if not 'H' in c.atomic_props['atoms']:
                self.excluded_complex_ids['Complex has no H'].append(c_id)

            if not 'C' in c.atomic_props['atoms']:
                self.excluded_complex_ids['Complex has no C'].append(c_id)

            if c.count_atoms_with_n_bonds(element='C', n_bonds=1) > 0:
                self.excluded_complex_ids['C atom with only 1 bond'].append(c_id)

            min_dist, _, _ = c.get_atomic_distances_between_atoms()
            if min_dist < 0.5:
                self.excluded_complex_ids['Atoms closer than 0.5A'].append(c_id)

            min_dist, _, _ = c.get_atomic_distances_between_atoms(skip_elements='H')
            if min_dist < 0.85:
                self.excluded_complex_ids['Heavy atoms closer than 0.85A'].append(c_id)

            if self.only_complexes_with_os and not c.has_metal_os():
                self.excluded_complex_ids['No metal OS'].append(c_id)

            if self.only_real_transition_metals and not c.metal_center_is_transition_metal():
                self.excluded_complex_ids['Metal center not transition metal'].append(c_id)

            if self.exclude_not_fully_connected_complexes and not c.fully_connected:
                self.excluded_complex_ids['Not fully connected'].append(c_id)

            if self.exclude_charged_complexes and c.charge != 0:
                self.excluded_complex_ids['Is charged'].append(c_id)

        # Exclude all invalid complexes.
        self.delete_filtered_complexes_from_db()

        return

    def print_excluded_complexes(self):
        print('Excluded complexes:')
        n_input_complexes = len(self.complex_db)
        for reason, c_ids in self.excluded_complex_ids.items():
            print(f'    - {reason}: {len(c_ids)} ({len(c_ids) / n_input_complexes*100:.2g}%)')

        print(f'  New number of input complexes: {len(self.complex_db)}')

        return

    def postfilter_if_ligands_valid(self, comp: RCA_Complex) -> Union[bool, str]:
        """
        Function for filtering complexes after the ligand extraction. Return False to exclude that complex.
        """
        if not hasattr(comp, 'ligands'):
            raise ValueError(f'The complex {comp.mol_id} has no attribute `ligands`.')
        elif comp.ligands == []:
            raise ValueError(f'The ligand list of the complex {comp.mol_id} is empty.')

        if comp.count_ligands_with_stoichiometry(atoms=['O'], only_connected=True) >= 3:
            return 'More than 3 O ligands'

        if comp.count_ligands_with_stoichiometry(atoms=['N'], only_connected=True) >= 2:
            return 'More than 2 N ligands'

        if comp.count_ligands_with_stoichiometry(atoms=['C']) > 0:
            return 'Ligand which is just C'

        if comp.count_n_unconnected_ligands(max_n_atoms=1) > 5:
            return 'More than 5 unconnected ligands'

        if comp.count_coordinating_atoms_with_distance_to_metal_greater_than(distance=1.9, element='O', max_n_atoms=1) > 0:
            return 'Oxygen ligand more than 1.9A away from metal'

        alkalis = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
        if comp.count_atoms_in_ligands(atoms=alkalis, only_if_connected_to_metal=True) > 0:
            return 'Alkali metal in ligand'

        noble_gases = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'Og']
        if comp.count_atoms_in_ligands(atoms=noble_gases, only_if_connected_to_metal=True) > 0:
            return 'Noble gas in ligand'

        heavy_metals = ['Tl', 'Pb', 'Bi', 'Po', 'Nh', 'Fl', 'Mc', 'Lv']
        if comp.count_atoms_in_ligands(atoms=heavy_metals, only_if_connected_to_metal=True) > 0:
            return 'Heavy metal in ligand'

        # Exclude likely cage structures of metals which are no real ligands.
        metals = ['B', 'Al', 'Ga', 'In', 'Tl', 'Nh', 'Si', 'Ge', 'Sn', 'Pb', 'Fl', 'As', 'Sb', 'Bi', 'Mc', 'Te', 'Po', 'Lv']
        if comp.count_ligands_containing_only(atoms=metals, denticity_range=[1, np.inf], n_atoms_range=[2, np.inf], except_elements=['H']):
            return 'Metal cage structures'

        return True

    def update_unique_ligand_db_with_database_info(self):
        """
        Updates the unique ligand database with information from the whole database.
        """
        needed_props = ['heavy_atoms_graph_hash_with_metal', 'graph_hash', 'pred_charge', 'n_hydrogens']
        data = {uname: {key: val for key, val in ulig.write_to_mol_dict(include_graph_dict=False).items() if key in needed_props} for uname, ulig in self.unique_ligand_db.db.items()}
        df = pd.DataFrame.from_dict(data, orient='index')

        all_most_common_n_H = df.groupby(['heavy_atoms_graph_hash_with_metal'])['n_hydrogens'].agg(lambda x: Counter(x).most_common(1)[0][0])
        all_charges_with_same_graph_hash = df.groupby('graph_hash')['pred_charge'].agg(lambda x: dict(Counter(x)))

        for uname, ulig in self.unique_ligand_db.db.items():
            ulig.same_graph_charges = all_charges_with_same_graph_hash[ulig.graph_hash]
            ulig.n_pred_charges = len(ulig.same_graph_charges)

            most_common_n_H = all_most_common_n_H[ulig.heavy_atoms_graph_hash_with_metal]
            ulig.common_graph_with_diff_n_hydrogens = bool(most_common_n_H != ulig.n_hydrogens)
            if ulig.common_graph_with_diff_n_hydrogens:
                ulig.add_warning(similar_molecule_with_diff_n_hydrogens_warning)

            ulig.n_electrons = ulig.n_protons - ulig.pred_charge
            ulig.odd_n_electron_count = bool(ulig.n_electrons % 2 == 1)
            if ulig.odd_n_electron_count:
                ulig.add_warning(odd_n_electrons_warning)

            ulig.identical_ligand_info['pred_charge'] = [self.full_ligand_db.db[name].pred_charge for name in ulig.identical_ligand_info['name']]
            ulig.identical_ligand_info['pred_charge_is_confident'] = [self.full_ligand_db.db[name].pred_charge_is_confident for name in ulig.identical_ligand_info['name']]

            ulig.has_warnings = bool(len(ulig.warnings) > 0)

        return

    def standardize_input_complex_db(self, complex_db):
        # TODO refactor these so that it outputs what is excluded
        complex_db.remove_node_features_from_molecular_graphs(keep=['node_label'])

        complex_db.normalize_multigraphs_into_simple_graphs()

        return complex_db

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

    def extract_ligands(self,
                        testing: Union[bool, int, list[str]] = False,
                        graph_creating_strategy: str = "default",
                        **kwargs
                        ):

        if isinstance(testing, list) is True:
            self.complex_db = ComplexDB.from_json(json_=str(self.input_complexes_json),
                                                   type_="Complex",
                                                   identifier_list=testing,
                                                   max_number=None,
                                                   graph_strategy=graph_creating_strategy,
                                                   **kwargs
                                                   )
        else:
            # testing in instance bool or integer (as before)
            self.complex_db = ComplexDB.from_json(json_=str(self.input_complexes_json),
                                                   type_="Complex",
                                                   identifier_list=None,
                                                   max_number=testing,
                                                   graph_strategy=graph_creating_strategy,
                                                   **kwargs
                                                   )


        self.prefilter_input_complex_db()
        self.complex_db = self.standardize_input_complex_db(self.complex_db)

        self.basic_ligand_infos = []
        for csd_code, comp in tqdm(self.complex_db.db.items(), desc="Extracting ligands from complexes"):
            comp.de_assemble()

            validity = self.postfilter_if_ligands_valid(comp)
            if validity == True:
                # Make a dataframe with basic ligand infos which doesn't take up much memory but has enough information to compute the grouping in unique ligands.
                ligand_infos = self.get_ligand_class_properties_of_complex(
                                                                            complex_=comp,
                                                                            props=['name', 'graph_hash', 'denticity', 'has_good_bond_orders', self.unique_ligand_id]
                                                                            )
                self.basic_ligand_infos.extend(ligand_infos)
            else:
                self.excluded_complex_ids[validity].append(csd_code)

        self.basic_ligand_infos = pd.DataFrame(self.basic_ligand_infos)

        self.delete_filtered_complexes_from_db()
        self.print_excluded_complexes()

        return

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

        grouped_ligands = self.basic_ligand_infos.groupby(groupby, sort=False).agg(list)

        return grouped_ligands


    @staticmethod
    def choose_unique_ligand_representative_from_all_same_ligands(same_ligands,
                                                                  strategy='good_bond_orders',
                                                                  ) -> str:

        name = None
        if strategy == 'most_common_denticity':

            # Counter ions/ solvent molecules should not count as most common denticity because we are not really interested in them.
            only_counter_ions = all([not lig.was_connected_to_metal for lig in same_ligands.values()])
            if only_counter_ions:
                denticities = [lig.denticity for lig in same_ligands.values()]
            else:
                denticities = [lig.denticity for lig in same_ligands.values() if lig.was_connected_to_metal]

            count_denticities = pd.Series(denticities).value_counts().sort_values(ascending=False)
            most_common_denticity = count_denticities.index[0]
            for name, lig in same_ligands.items():
                if lig.denticity == most_common_denticity:
                    break
        elif strategy == 'first':
            # Just take the first entry.
            name = list(same_ligands.keys())[0]
        elif strategy == 'good_bond_orders':
            # Preferably take ligands which have good bond orders
            name = list(same_ligands.keys())[0]
            for lig_name, lig in same_ligands.items():
                if lig.has_good_bond_orders:
                    name = lig_name
                    break
        else:
            raise ValueError(
                f'Unknown strategy `{strategy}` to choose the unique ligand representative from all same ligands.')

        return name

    def build_unique_ligand_db(self):
        self.grouped_ligands = self.group_same_ligands()
        self.graph_hash_grouped_ligands = self.group_same_ligands(groupby='graph_hash')

        self.unique_ligand_db = {}
        for grouped_ligand in tqdm(self.grouped_ligands.to_dict(orient='index').values(), desc="Building unique ligand db"):
            same_ligands_names = grouped_ligand['name']
            same_ligands = {name: self.full_ligand_db.db[name] for name in same_ligands_names}
            name = self.choose_unique_ligand_representative_from_all_same_ligands(same_ligands=same_ligands)

            unique_ligand = deepcopy(self.full_ligand_db.db[name])

            uname = 'unq_' + name
            unique_ligand.unique_name = uname


            denticities = pd.unique(self.graph_hash_grouped_ligands.loc[unique_ligand.graph_hash]['denticity'])
            metals = [self.full_ligand_db.db[ligand_name].original_metal_symbol for ligand_name in same_ligands_names]

            # Add useful statistical information of all ligands for this unique ligand
            count_metals = pd.Series(metals).value_counts().sort_values(ascending=False).to_dict()
            n_graph_hashes = len(pd.unique(self.graph_hash_grouped_ligands.loc[unique_ligand.graph_hash]['graph_hash_with_metal']))

            assert not 0 in denticities, 'The denticity for unconnected ligands is assumed to be -1 but here there appears a 0.'
            has_unconnected_ligands = -1 in denticities

            identical_ligand_info = defaultdict(list)
            for lig in same_ligands.values():
                identical_ligand_info['name'].append(lig.name)
                identical_ligand_info['original_metal_symbol'].append(lig.original_metal_symbol)
                identical_ligand_info['original_metal_os'].append(lig.original_metal_os)
                identical_ligand_info['original_complex_charge'].append(lig.original_metal_os)
                identical_ligand_info['original_complex_id'].append(lig.global_props['CSD_code'])
            unique_ligand.identical_ligand_info = identical_ligand_info

            unique_ligand_infos = {
                'occurrences': len(same_ligands_names),
                'same_graph_denticities': denticities,
                'count_metals': count_metals,
                'n_same_graph_denticities': len(denticities),
                'n_metals': len(count_metals),
                'n_same_graphs': n_graph_hashes,
                'has_unconnected_ligands': has_unconnected_ligands,
                'all_ligands_metals': metals,
            }
            for prop, val in unique_ligand_infos.items():
                setattr(unique_ligand, prop, val)
            self.unique_ligand_info_props = list(
                unique_ligand_infos.keys())  # for updating the ligands from complex and full ligands db later

            unique_ligand.all_ligand_names = same_ligands_names

            self.unique_ligand_db[uname] = unique_ligand
        self.unique_ligand_db = LigandDB(self.unique_ligand_db)

        # Add useful property for later
        self.ligand_to_unique_ligand = {}
        for uname, ulig in self.unique_ligand_db.db.items():
            for name in ulig.all_ligand_names:
                self.ligand_to_unique_ligand[name] = uname

        return

    @staticmethod
    def update_ligand_with_unique_ligand_information_inplace(lig, ulig, share_properties=None,
                                                             share_global_props=None,
                                                             collect_properties=None):
        if collect_properties is None:
            collect_properties = {}
        if share_global_props is None:
            share_global_props = []
        if share_properties is None:
            share_properties = []

        for prop in share_properties:
            value = deepcopy(getattr(ulig, prop))
            setattr(lig, prop, value)
        update_dict_with_warning_inplace(lig.global_props, deepcopy(ulig.global_props), share_global_props)

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
                                           share_global_props: list = [],
                                           collect_properties: dict = {}
                                           ):
        self.ensure_complex_db()
        self.ensure_unique_ligand_db()

        # Update ligands with unique ligand information
        for c in tqdm(self.complex_db.db.values(), 'Update complex db with unique ligand information'):
            for lig in c.ligands:
                uname = self.ligand_to_unique_ligand[lig.name]
                ulig = self.unique_ligand_db.db[uname]
                self.update_ligand_with_unique_ligand_information_inplace(
                    lig=lig,
                    ulig=ulig,
                    share_properties=share_properties,
                    share_global_props=share_global_props,
                    collect_properties=collect_properties
                )

            # Update global props with some useful information
            c.global_props['n_ligands'] = len(c.ligands)
            c.global_props['n_unique_ligands'] = len(set(lig.unique_name for lig in c.ligands))
            n_ligands_occurring_once = sum(
                [lig.unique_ligand_information['occurrences'] == 1 for lig in c.ligands])
            c.global_props['n_ligands_occurring_once'] = n_ligands_occurring_once
            c.global_props['frac_ligands_occurring_once'] = n_ligands_occurring_once / len(c.ligands)

        return

    def update_full_ligand_db_with_information(self,
                                               share_properties: list = [],
                                               share_global_props: list = [],
                                               collect_properties: dict = {}
                                               ):
        self.ensure_full_ligand_db()
        self.ensure_unique_ligand_db()

        for lig in tqdm(self.full_ligand_db.db.values(), 'Update full ligand db with unique ligand information'):
            uname = self.ligand_to_unique_ligand[lig.name]
            ulig = self.unique_ligand_db.db[uname]
            self.update_ligand_with_unique_ligand_information_inplace(
                lig=lig,
                ulig=ulig,
                share_properties=share_properties,
                share_global_props=share_global_props,
                collect_properties=collect_properties
            )

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
                with open(self.input_complexes_json, "r") as file:
                    n_complexes_in_json = len(json.load(file))

                if not self.testing is None and (self.testing > n_complexes_in_json or self.testing == False):
                    print(
                        f'WARNING: Cannot use existing input json of complexes because it contains less complexes than the testing parameter. Reload xzy, global properties and graph data instead.')
                    self.load_input_data_to_json(overwrite_atomic_properties=overwrite_atomic_properties, **kwargs)
        else:
            self.load_input_data_to_json(overwrite_atomic_properties=overwrite_atomic_properties, **kwargs)

        return

    def get_complex_dict_for_LCS(self):
        """
        Returns a dictionary of all complexes with only the needed properties for the LCS.
        """
        needed_complex_props = ['mol_id', 'metal_oxi_state', 'charge']
        needed_ligand_props = ['stoichiometry', 'name', 'unique_name', 'n_protons']

        charge_complexes = {}
        for c_id, c in self.complex_db.db.items():
            charge_complexes[c_id] = {prop: getattr(c, prop) for prop in needed_complex_props}

            charge_complexes[c_id]['ligands'] = []
            for lig in c.ligands:
                lig_props = {prop: getattr(lig, prop) for prop in needed_ligand_props}
                charge_complexes[c_id]['ligands'].append(lig_props)

        return charge_complexes

    def calculate_ligand_charges(self, max_iterations=None):
        charge_complexes = self.get_complex_dict_for_LCS()
        df_ligand_charges = get_charges_of_unique_ligands(all_complexes=charge_complexes, max_iterations=max_iterations)

        return df_ligand_charges

    def run_ligand_extraction(self,
                              calculate_charges: bool = True,
                              overwrite_atomic_properties: bool = True,
                              use_existing_input_json: bool = True,
                              get_only_unique_ligand_db_without_charges: bool = False,
                              max_charge_iterations: Union[int, None] = None,
                              **kwargs
                              ):
        """
        Runs the entire ligand extraction process from reading in the .xzy files to optionally assigning charges.
        """
        start = datetime.now()

        self.ensure_input_complex_db_exists(overwrite_atomic_properties=overwrite_atomic_properties,
                                            use_existing_input_json=use_existing_input_json,
                                            **kwargs)

        self.extract_ligands(testing=self.testing,
                             graph_creating_strategy=self.graph_strat,
                             **kwargs
                             )
        if self.save_intermediate_databases:
            self.complex_db.to_json(path=self.output_complexes_json, desc='Save intermediate complex db')

        # Attention: The ligands of the full ligand database are just a view on the ligands of the complex db. This is done for memory efficiency and speed. To change this, it is not enough to just set copy to True, but the full ligand db needs to be updated with unique ligand information and with charges.
        self.full_ligand_db = self.build_full_ligand_db(copy=False)
        if self.save_intermediate_databases:
            self.full_ligand_db.to_json(path=self.full_ligands_json, desc='Save intermediate full ligand db')

        self.build_unique_ligand_db()
        if self.save_intermediate_databases:
            self.unique_ligand_db.to_json(path=self.unique_ligands_json, desc='Save intermediate unique ligand db')

        if not get_only_unique_ligand_db_without_charges:
            # Update complex db to include information about the unique ligands for the LCS.
            share_properties = ['unique_name']
            collect_properties = {'unique_ligand_information': self.unique_ligand_info_props}
            self.update_complex_db_with_information(share_properties=share_properties,
                                                    collect_properties=collect_properties)
            if self.save_intermediate_databases:
                self.complex_db.to_json(path=self.output_complexes_json, desc='Save updated complex database')

            if self.save_intermediate_databases:
                self.full_ligand_db.to_json(path=self.full_ligands_json, desc='Save updated full ligand database')

            # Charge assignment using only the linear charge solver (LCS) right now
            if calculate_charges:
                print('\nCHARGE CALCULATION:')
                df_ligand_charges = self.calculate_ligand_charges(max_iterations=max_charge_iterations)
                self.update_databases_with_charges(df_ligand_charges=df_ligand_charges)
        else:
            calculate_charges = False

        self.update_unique_ligand_db_with_database_info()

        self.unique_ligand_db.to_json(self.unique_ligands_json, desc='Save unique ligand db to json')
        self.full_ligand_db.to_json(self.full_ligands_json, desc='Save full ligand db to json')
        self.complex_db.to_json(self.output_complexes_json, desc='Save complex db to json')

        with_charges = 'with charges' if calculate_charges else 'without charges'
        updates = ' Complex db and full ligand db were not updated with unique ligand information.' if get_only_unique_ligand_db_without_charges else ''
        print(f'\nLigand database {with_charges} established successfully!{updates}')

        duration = get_duration_string(start)
        print(f'Duration of extraction: {duration}')

        return

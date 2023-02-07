"""
Class for extracting ligands from a database of complexes.
"""
import pandas as pd
from copy import deepcopy
from src01.DataBase import MoleculeDB
from pathlib import Path
from tqdm import tqdm
import json
from typing import Union
from datetime import datetime

from src01.DataLoader import DataLoader
from src01.io_custom import load_unique_ligand_db, load_complex_db, load_full_ligand_db, save_unique_ligand_db, \
    save_full_ligand_db, save_complex_db
from src01.utilities_extraction import get_charges_of_unique_ligands, update_ligand_with_charge_inplace
from src01.utilities import update_dict_with_warning_inplace


class LigandExtraction:

    def __init__(self, database_path: str,
                 data_store_path: str,
                 exclude_not_fully_connected_complexes: bool = True,
                 testing: Union[bool, int] = False,
                 graph_strat: str = "default"
                 ):

        self.ligand_to_unique_ligand = None
        self.unique_ligand_info_props = None
        self.grouped_ligands = None
        self.database_path = ""
        self.data_store_path = ""
        self.exclude_not_fully_connected_complexes = None
        self.testing = None
        self.complex_db = None
        self.all_hashes = {}

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
            print(f'Path given as `data_store_path` ({data_store_path} doesn\'t exist yet. Making this directory.')
            data_store_path.mkdir(parents=True, exist_ok=True)

        if not (isinstance(testing, int) or isinstance(testing, bool)):
            raise ValueError(f'Input variable `testing` must be int or bool but is {type(testing)}.')

        self.database_path = database_path
        self.data_store_path = data_store_path
        self.exclude_not_fully_connected_complexes = exclude_not_fully_connected_complexes
        self.testing = testing

        return

    def load_input_data_to_json(self, overwrite_atomic_properties: bool = False):
        """
        Establish and safe the Database (in our case tmQM) as json for simple loading.
        """
        db_dict = DataLoader(database_path_=self.database_path, overwrite=overwrite_atomic_properties).data_for_molDB
        tmQMG_DB = MoleculeDB.from_json(json_=db_dict, type_="Complex", max_number=self.testing, graph_strategy=self.graph_strat)

        tmQMG_DB.to_json(path=self.input_complexes_json)

        return

    def filter_input_complex_db(self, complex_db):
        if self.exclude_not_fully_connected_complexes:
            complex_db.filter_not_fully_connected_molecules()
        complex_db.remove_node_features_from_molecular_graphs(keep=['node_label'])
        complex_db.remove_edge_features_from_molecular_graphs()
        complex_db.normalize_multigraphs_into_simple_graphs()

        return complex_db

    def extract_ligands(self,
                        testing: Union[bool, int] = False,
                        graph_creating_strategy: str = "default"
                        ):

        self.complex_db = MoleculeDB.from_json(json_=str(self.input_complexes_json),
                                               type_="Complex",
                                               identifier_list=None,
                                               max_number=testing,
                                               graph_strategy=graph_creating_strategy
                                               )

        #
        # TODO Instead of filtering the whole db, implement to check each complex individually
        self.complex_db = self.filter_input_complex_db(self.complex_db)

        self.all_hashes = {}
        for csd_code, complex_ in tqdm(self.complex_db.db.items(), desc="Extracting ligands from complexes"):
            complex_.de_assemble()

            graph_hashes = {lig.name: lig.graph_hash for lig in complex_.ligands}
            self.all_hashes.update(graph_hashes)

        self.complex_db.to_json(path=self.output_complexes_json)

        return

    def get_full_ligand_db(self):
        full_ligands = {}
        for c in self.complex_db.db.values():
            for lig in c.ligands:
                name = lig.name
                full_ligands[name] = lig
        full_ligand_db = MoleculeDB(full_ligands)
        return full_ligand_db

    def save_full_ligand_db(self):
        full_ligand_db = self.get_full_ligand_db()
        full_ligand_db.to_json(self.full_ligands_json)

        return

    def group_ligands_by_hash(self):
        grouped_ligands_by_hash = {}
        for k, v in self.all_hashes.items():
            grouped_ligands_by_hash.setdefault(v, []).append(k)

        return grouped_ligands_by_hash

    def get_full_ligands_dict(self):
        with open(self.output_complexes_json, 'r') as file:
            complexes = json.load(file)

        ligand_dict = {}
        for c in complexes.values():
            for lig in c['ligands']:
                ligand_dict[lig['name']] = lig

        return ligand_dict

    @staticmethod
    def choose_unique_ligand_representative_from_all_same_ligands(same_ligands,
                                                                  strategy='most_common_denticity'
                                                                  ) -> str:
        #
        name = None
        #
        if strategy == 'most_common_denticity':
            denticities = [lig['denticity'] for lig in same_ligands.values()]
            count_denticities = pd.Series(denticities).value_counts().sort_values(ascending=False)
            most_common_denticity = count_denticities.index[0]
            for name, lig in same_ligands.items():
                if lig['denticity'] == most_common_denticity:
                    break
        elif strategy == 'first':
            name = list(same_ligands.keys())[0]
        else:
            raise ValueError(
                f'Unknown strategy `{strategy}` to choose the unique ligand representative from all same ligands.')

        return name

    def build_unique_ligand_db(self):
        ligand_dict = self.get_full_ligands_dict()

        self.grouped_ligands = self.group_ligands_by_hash()

        unique_ligand_dict = {}
        for same_ligands_names in tqdm(self.grouped_ligands.values(), desc="filling unique ligand dict"):
            same_ligands = {name: ligand_dict[name] for name in same_ligands_names}
            name = self.choose_unique_ligand_representative_from_all_same_ligands(same_ligands=same_ligands)
            uname = 'unq_' + name
            unique_ligand = deepcopy(ligand_dict[name])

            unique_ligand['unique_name'] = uname

            denticities = [ligand_dict[ligand_name]['denticity'] for ligand_name in same_ligands_names]
            metals = [ligand_dict[ligand_name]['original_metal_symbol'] for ligand_name in same_ligands_names]

            # Add useful statistical information of all ligands for this unique ligand
            count_denticities = pd.Series(denticities).value_counts().sort_values(ascending=False).to_dict()
            count_metals = pd.Series(metals).value_counts().sort_values(ascending=False).to_dict()
            chosen_denticity_fraction = sum([dent == unique_ligand['denticity'] for dent in denticities]) / len(
                denticities)

            assert not 0 in count_denticities, 'The denticity for unconnected ligands is assumed to be -1 but here there appears a 0.'
            has_unconnected_ligands = -1 in count_denticities.keys()

            unique_ligand_infos = {
                'occurrences': len(same_ligands_names),
                'count_denticities': count_denticities,
                'count_metals': count_metals,
                'n_denticities': len(count_denticities),
                'n_metals': len(count_metals),
                'chosen_denticity_fraction': chosen_denticity_fraction,
                'has_unconnected_ligands': has_unconnected_ligands
            }
            unique_ligand.update(unique_ligand_infos)
            self.unique_ligand_info_props = list(
                unique_ligand_infos.keys())  # for updating the ligands from complex and full ligands db later

            unique_ligand['all_ligand_names'] = same_ligands_names

            # Delete attribute original metal from unique_ligand since it is confusing and no real attribute of a unique ligand
            del unique_ligand['original_metal']
            del unique_ligand['original_metal_symbol']

            unique_ligand_dict[uname] = unique_ligand

        with open(self.unique_ligands_json, 'w') as file:
            json.dump(unique_ligand_dict, file)
            print(f"Unique ligand database saved to {self.unique_ligands_json}.")

        # Add useful property for later
        self.ligand_to_unique_ligand = {}
        for uname, ulig in unique_ligand_dict.items():
            for name in ulig['all_ligand_names']:
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

        ulig = deepcopy(ulig)

        update_dict_with_warning_inplace(lig, ulig, share_properties)
        update_dict_with_warning_inplace(lig['global_props'], ulig['global_props'], share_global_props)

        # Collect properties from unique ligand in a dictionary in the full ligands.
        for new_prop, old_props in collect_properties.items():
            lig[new_prop] = {prop: ulig[prop] for prop in old_props}

        lig['is_chosen_unique_ligand'] = ulig['name'] == lig['name']

        return

    def update_complex_db_with_information(self, share_properties=None,
                                           share_global_props=None,
                                           collect_properties=None
                                           ):
        if collect_properties is None:
            collect_properties = {}
        if share_global_props is None:
            share_global_props = []
        if share_properties is None:
            share_properties = []

        print('Update complex db with unique ligand information.')
        complexes = load_complex_db(self.output_complexes_json)
        unique_ligands = load_unique_ligand_db(self.unique_ligands_json)

        # Update ligands with unique ligand information
        for c in complexes.values():
            for lig in c['ligands']:
                uname = self.ligand_to_unique_ligand[lig['name']]
                ulig = unique_ligands[uname]
                self.update_ligand_with_unique_ligand_information_inplace(
                    lig=lig,
                    ulig=ulig,
                    share_properties=share_properties,
                    share_global_props=share_global_props,
                    collect_properties=collect_properties
                )

            # Update global props with some useful information
            c['global_props']['n_ligands'] = len(c['ligands'])
            n_ligands_occurring_once = sum(
                [lig['unique_ligand_information']['occurrences'] == 1 for lig in c['ligands']])
            c['global_props']['n_ligands_occurring_once'] = n_ligands_occurring_once
            c['global_props']['frac_ligands_occurring_once'] = n_ligands_occurring_once / len(c['ligands'])

        print(f'Save updated complex database to {self.output_complexes_json}')
        save_complex_db(db=complexes, path=self.output_complexes_json)

        return

    def update_full_ligand_db_with_information(self, share_properties=None,
                                               share_global_props=None,
                                               collect_properties=None
                                               ):
        if collect_properties is None:
            collect_properties = {}
        if share_global_props is None:
            share_global_props = []
        if share_properties is None:
            share_properties = []

        print('Update full ligand db with unique ligand information.')
        full_ligands = load_full_ligand_db(self.full_ligands_json)
        unique_ligands = load_unique_ligand_db(self.unique_ligands_json)

        for lig in full_ligands.values():
            uname = self.ligand_to_unique_ligand[lig['name']]
            ulig = unique_ligands[uname]
            self.update_ligand_with_unique_ligand_information_inplace(
                lig=lig,
                ulig=ulig,
                share_properties=share_properties,
                share_global_props=share_global_props,
                collect_properties=collect_properties
            )
        print(f'Save updated full ligand database to {self.full_ligands_json}')
        save_full_ligand_db(db=full_ligands, path=self.full_ligands_json)

        return

    def update_databases_with_charges(self, df_ligand_charges: pd.DataFrame, ):
        df_ligand_charges = df_ligand_charges.set_index('unique_name')
        charges = df_ligand_charges.to_dict(orient='index')

        unique_ligands = load_unique_ligand_db(path=self.unique_ligands_json)
        not_intersecting_ligands = set(unique_ligands.keys()).symmetric_difference(set(charges.keys()))
        print(f'Number of unique ligands not outputted by the charge assignment: {len(not_intersecting_ligands)}')
        for ulig in unique_ligands.values():
            update_ligand_with_charge_inplace(ulig, charges)
        save_unique_ligand_db(db=unique_ligands, path=self.unique_ligands_json)

        ligands = load_full_ligand_db(self.full_ligands_json)
        for lig in ligands.values():
            update_ligand_with_charge_inplace(lig, charges)
        save_full_ligand_db(db=ligands, path=self.full_ligands_json)

        complexes = load_complex_db(path=self.output_complexes_json)
        for c in complexes.values():
            for lig in c['ligands']:
                update_ligand_with_charge_inplace(lig, charges)
        save_complex_db(db=complexes, path=self.output_complexes_json)

        return

    def ensure_input_complex_db_exists(self, overwrite_atomic_properties: bool, use_existing_input_json: bool):
        if use_existing_input_json:
            if not self.input_complexes_json.exists():
                print(
                    f'WARNING: Cannot use existing input json of complexes because path not found: {self.input_complexes_json}. Reload xzy, global properties and graph data instead.')
                self.load_input_data_to_json(overwrite_atomic_properties=overwrite_atomic_properties)
        else:
            self.load_input_data_to_json(overwrite_atomic_properties=overwrite_atomic_properties)

        return

    def run_ligand_extraction(self,
                              calculate_charges: bool = True,
                              overwrite_atomic_properties: bool = True,
                              use_existing_input_json: bool = True,
                              get_only_unique_ligand_db_without_charges: bool = False
                              ):
        """
        Runs the entire ligand extraction process from reading in the .xzy files to optionally assigning charges.
        """
        start = datetime.now()

        self.ensure_input_complex_db_exists(overwrite_atomic_properties=overwrite_atomic_properties,
                                            use_existing_input_json=use_existing_input_json)

        self.extract_ligands(testing=self.testing,
                             graph_creating_strategy=self.graph_strat
                             )

        self.save_full_ligand_db()
        self.build_unique_ligand_db()

        if not get_only_unique_ligand_db_without_charges:

            # Update complex db to include information about the unique ligands for the LCS.
            share_properties = ['unique_name']
            collect_properties = {'unique_ligand_information': self.unique_ligand_info_props}
            self.update_complex_db_with_information(share_properties=share_properties,
                                                    collect_properties=collect_properties)
            self.update_full_ligand_db_with_information(share_properties=share_properties,
                                                        collect_properties=collect_properties)

            # Charge assignment using only the linear charge solver (LCS) right now
            if calculate_charges:
                print('\nCHARGE CALCULATION:')
                df_ligand_charges = get_charges_of_unique_ligands(all_complexes_path=self.output_complexes_json)
                self.update_databases_with_charges(df_ligand_charges=df_ligand_charges)
                print('Charges successfully assigned.')
        else:
            calculate_charges = False

        with_charges = 'with charges' if calculate_charges else 'without charges'
        updates = ' Complex db and full ligand db were not updated with unique ligand information' if get_only_unique_ligand_db_without_charges else ''
        print(f'\nLigand database {with_charges} established successfully!{updates}')

        end = datetime.now()
        print(f'Duration of extraction: {end - start}')

        return

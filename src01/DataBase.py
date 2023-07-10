import functools
import json
from collections import Counter
from copy import deepcopy
from tqdm import tqdm
import pandas as pd
import numpy as np
from src01.Molecule import RCA_Molecule, RCA_Complex, RCA_Ligand    # important, don't delete, read in by string that's why they appear grey
from src01.utilities_graph import remove_node_features_from_graph, make_multigraph_to_graph, remove_edge_features_from_graph
from src01.utilities import identify_metal_in_ase_mol
from src01.utilities_Molecule import get_all_ligands_by_graph_hashes, group_list_without_hashing
import networkx as nx
from pymatgen.core.periodic_table import Element as Pymatgen_Element
from src01.io_custom import save_json, NumpyEncoder, load_json, load_unique_ligand_db, load_complex_db, iterate_over_json
from scipy.special import comb
from typing import Union
from datetime import datetime
from pathlib import Path
import jsonlines
from memory_profiler import profile
import itertools

from src11_machine_learning.dataset_preparation.descriptors import SOAP_3D


class BaselineDB:
    def __init__(self, dict_: dict):
        """
        Is basically just a dict with extended functionality
        :param dict_: is dict as {identifier: RCA_Mol/RCA_Ligand}
        So the base functionality will assume that we have an RCA_Mol, which assures
        that the functionality is stable even if we use RCA_Ligands (As they are a subclass of RCA_Mol)
        """
        self.db = dict_
        self.names = list(self.db.keys())
        self.reduced_df = self.get_reduced_df()

    def __len__(self):
        return len(self.db)

    def __eq__(self, other):
        return self.db == other.db

    def get_reduced_df(self):
        important_columns = ['name', 'stoichiometry', 'denticity', 'graph_hash_with_metal', 'unique_name', 'pred_charge', 'pred_charge_is_confident']
        data = []
        for name, mol in self.db.items():
            props = {}
            for prop in important_columns:
                try:
                    props[prop] = getattr(mol, prop)
                except AttributeError:
                    pass
            data.append(props)
        df = pd.DataFrame(data)

        return df

    def get_first_entry(self):
        first_key = list(self.db.keys())[0]
        first_item = self.db[first_key]

        return first_item

    def append_DB(self, key, molecule, overwrite=True):
        if key in self.db and overwrite is False:
            print("Key already in DB and no overwrite, hence nothing changed")
        else:
            self.db[key] = molecule

    def get_dict_in_json_format(self, desc: str='Convert db to dict'):
        json_dict = {}
        for key, mol in tqdm(self.db.items(), desc):
            json_dict[key] = mol.write_to_mol_dict()
        return json_dict

    def to_json(self, path, desc: str='Save DB to json', json_lines: bool=False):
        if json_lines:
            with jsonlines.open(path, mode='w', dumps=functools.partial(json.dumps, cls=NumpyEncoder)) as writer:
                for key, mol in tqdm(self.db.items(), desc):
                    data = {'key': key, 'value': mol.write_to_mol_dict()}
                    writer.write(data)
        else:
            d = self.get_dict_in_json_format(desc=desc)
            save_json(d, path=path, indent=4)

        return

    def filter_not_fully_connected_molecules(self):
        deleted_identifiers = []
        for identifier, mol in self.db.items():
            fully_connected = nx.is_connected(mol.graph)
            if not fully_connected:
                deleted_identifiers.append(identifier)

        for identifier in deleted_identifiers:
            del self.db[identifier]

        print(f'Deleted {len(deleted_identifiers)} molecules from input because they were not fully connected.')
        return

    def remove_node_features_from_molecular_graphs(self, keep: list = []):
        """
        Removes all node features from all molecular graphs in the db except the node features specified in keep.
        :param keep: list of node features which will not be removed
        :return: None
        """
        for identifier, mol in self.db.items():
            remove_node_features_from_graph(
                graph=mol.graph,
                keep=keep,
                inplace=True
            )

        print(f'Removed all node features from all graphs except: {", ".join(keep)}')
        return

    def remove_edge_features_from_molecular_graphs(self, keep: list = []):
        """
        Removes all edge features from all molecular graphs in the db except the edge features specified in keep.
        :param keep: list of edge features which will not be removed
        :return: None
        """
        for identifier, mol in self.db.items():
            remove_edge_features_from_graph(
                graph=mol.graph,
                keep=keep,
                inplace=True
            )

        print(f'Removed all edge features from all graphs except: {", ".join(keep)}')
        return

    def normalize_multigraphs_into_simple_graphs(self):
        """
        Makes all molecular multigraphs to graphs.
        :return: None
        """
        for identifier, mol in tqdm(self.db.items(), 'Standardize graphs to simple nx.Graph objects'):
            mol.graph = make_multigraph_to_graph(mol.graph)

        return

class MoleculeDB(BaselineDB):
    type = 'Molecule'
    def __init__(self, dict_):
        super().__init__(dict_=dict_)

    def check_db_equal(self, db: str):
        db = MoleculeDB.from_json(json_=db, type_='Molecule')
        return self == db

    @classmethod
    def from_json(cls,
                  json_,
                  type_: str = "Molecule",
                  graph_strategy: Union[str, None] = None,
                  max_number: Union[int, list[str], None] = None,
                  show_progress: bool = True,
                  **kwargs
                  ):
        """
        :param json_: path to json file
        :param type_: ignored. Only for historical reasons
        :param graph_strategy: strategy to create the graph
        :param max_number: maximum number of molecules to read in
        :param show_progress: show progress bar
        :param kwargs: kwargs for the graph creation
        """
        db = {}
        if isinstance(json_, dict):
            for n, (identifier, mol) in tqdm(enumerate(json_.items()), desc=f"Build {cls.type} Database", disable=not show_progress):
                if isinstance(max_number, int):
                    if n >= max_number:
                        break
                elif isinstance(max_number, (list, tuple, set)):
                    if identifier not in max_number:
                        continue
                db[identifier] = globals()[f"RCA_{cls.type}"].read_from_mol_dict(
                                                                                    dict_=mol,
                                                                                    graph_creating_strategy=graph_strategy,
                                                                                    csd_code=identifier,
                                                                                    **kwargs
                                                                                    )
        else:
            for identifier, mol in tqdm(iterate_over_json(path=json_, n_max=max_number, show_progress=False),
                                        desc=f"Build {cls.type} Database", disable=not show_progress):
                db[identifier] = globals()[f"RCA_{cls.type}"].read_from_mol_dict(
                                                                                dict_=mol,
                                                                                graph_creating_strategy=graph_strategy,
                                                                                csd_code=identifier,
                                                                                **kwargs
                                                                                )

        return cls(db)

class LigandDB(MoleculeDB):
    type = 'Ligand'
    def __init__(self, dict_):
        super().__init__(dict_=dict_)

    @classmethod
    def load_from_json(cls, path: Union[str, Path], n_max: int=None, show_progress: bool=True, only_core_ligands: bool=False):
        """
        Load a JSON or JSON Lines file.
        :param path: Path to the JSON or JSON Lines file
        :return: A LigandDB object
        """
        db = load_unique_ligand_db(path=path, n_max=n_max, show_progress=show_progress, molecule='class')
        if only_core_ligands:
            # Remove all free ligands which were not connected to the metal.
            db = {identifier: lig for identifier, lig in db.items() if lig.denticity > 0}

        return cls(db)


    def check_db_equal(self, db: str) -> bool:
        """
        Checks if two ligand databases are equal.
        :param db: Either the LigandDB itself or the path to a json file
        """
        try:
            db = LigandDB.from_json(json_=db, type_='Ligand')
        except ValueError:
            pass
        return self == db

    def get_SOAP_descriptors(self, r_cut=50.0, n_max=2, l_max=1, crossover=False, average='inner', weighting=None, only_from_donors: bool=False):
        """
        Adds SOAP features to all molecules in the database.
        :return: SOAP descriptors.
        """
        ase_ligands = [lig.mol for lig in self.db.values()]
        soap = SOAP_3D(
            ase_molecules=ase_ligands,
            r_cut=r_cut,
            n_max=n_max,
            l_max=l_max,
            crossover=crossover,
            average=average,
            weighting=weighting
        )

        if only_from_donors:
            positions = [lig.get_atomic_positions(atomic_indices=lig.ligand_to_metal) for lig in self.db.values()]
        else:
            positions = None

        soap_desc = soap.calculate_descriptors(positions=positions)

        # Add SOAP descriptors to molecules
        for i, lig in enumerate(self.db.values()):
            lig.soap = soap_desc[i]

        # Create dataframe with SOAP descriptors
        soap_labels, mol_names = soap.get_descriptor_names(), list(self.db.keys())
        df = pd.DataFrame(data=soap_desc, columns=soap_labels, index=mol_names)

        return df

    @classmethod
    def from_MoleculeDB(cls,
                        molDB: MoleculeDB,
                        metals_of_interest=None,   # important for assembly
                        denticity_numbers_of_interest: list = None, # important for assembly
                        Testing=False):
        """
        molDB: Database of Molecules to extract ligands from
        metals of interest: if None, we want all. (mainly if we want to cut down the extraction set for test reasons)
        denticity: if none we want all (mainly if we want to cut down the extraction set for test reasons)
        """
        if metals_of_interest is None:
            metals_of_interest = []

        if isinstance(denticity_numbers_of_interest, int):
            denticity_numbers_of_interest = [denticity_numbers_of_interest]
        elif denticity_numbers_of_interest is None:
            # then we want all denticities
            denticity_numbers_of_interest = list(range(1, 9999))

        #
        # actual extraction:
        ligand_list = []

        for csd_code, molecule in tqdm(molDB.db.items(), desc="Extracting Ligands from Molecules"):

            mol_metal = identify_metal_in_ase_mol(molecule.mol)
            if mol_metal in metals_of_interest or metals_of_interest == []:
                molecule.de_assemble(Testing=Testing)
                for lig in molecule.ligands:
                    if lig.denticity in denticity_numbers_of_interest:
                        ligand_list.append(lig)
                    else:
                        print('Denticity over 10!')         # TODO remove
                        pass

        return cls(dict_={lig.name: lig for lig in ligand_list})

    # Just in case we need the old format
    def get_lig_db_in_old_format(self):
        """
        sort it by denticity
        """
        lig_dict_old_format = {}

        for lig in self.db.values():

            if lig.denticity not in lig_dict_old_format:
                lig_dict_old_format[lig.denticity] = [lig]
            else:
                lig_dict_old_format[lig.denticity].append(lig)

        return lig_dict_old_format

    # from here on we provide the necessary tools for the duplicant identification
    def hash_check(self):
        """
        returns True, if there exists two ligands with same graph hash, but different stoichiometry
        this means, that the graph hash is not unique on the ligand set,
        and an exact graph comparison is required
        Otherwise, as we assume that for ligands with the same stoichiometric the graph hashes will
        be different (because they are too close for the hashes to be equal)
        and thus graph_hash equal if and only if molecules equal
        """

        all_ligands = list(self.db.values())
        all_ligands_copy = deepcopy(all_ligands)

        while len(all_ligands_copy) > 1:
            lig = all_ligands_copy.pop()
            for lig2 in all_ligands_copy:
                if lig.graph_hash == lig2.graph_hash and lig.stoichiometry != lig2.stoichiometry:
                    return True

        return False

    def exact_comparison(self):

        ligands_by_hash = get_all_ligands_by_graph_hashes(list(self.db.values()))

        grouped_unique_ligands = []
        for graph_hash, ligand_list in tqdm(ligands_by_hash.items(), desc='Compare graphs exact'):
            unique_hash_ligand_list = group_list_without_hashing(ligand_list)
            grouped_unique_ligands.extend(unique_hash_ligand_list)

        return grouped_unique_ligands

    @staticmethod
    def check_property_and_print_if_not_same_for_all_same_ligands(check_props, unique_ligand, ligand):
        for prop in check_props:
            if getattr(unique_ligand, prop) != getattr(ligand, prop):
                print(
                    f'WARNING: Different {prop} for unique ligand {unique_ligand.name} ({getattr(unique_ligand, prop)}) and ligand {ligand.name} ({getattr(ligand, prop)}).')

        return

    def get_unique_ligands_and_set_unique_ligand_name(self, grouped_unique_ligands):
        unique_ligands = []
        for same_ligands in grouped_unique_ligands:

            unique_ligand = same_ligands[0]
            unique_ligand_name = 'unq_' + unique_ligand.name

            for ligand in tqdm(same_ligands, desc="Filter Duplicates"):
                ligand.unique_name = unique_ligand_name
                ligand.n_total_unique_ligands = len(same_ligands)

                check_props = ['denticity', 'graph_hash', 'hash', 'unique_name']
                self.check_property_and_print_if_not_same_for_all_same_ligands(check_props, unique_ligand, ligand)

            unique_ligands.append(deepcopy(unique_ligand))

        return unique_ligands

    # now comes the actual duplicate Filtering
    def filter_duplicates(self) -> dict:

        print('Start filtering duplicates.')
        # all_ligands = list(self.db.values())

        # Calculate graphs and graph_hashes for all ligands and save as attribute.
        # -> no longer required as the graph hashes are already computed during the creation of a ligand
        # all_graph_hashes = [lig.graph_hash for lig in all_ligands]

        if self.hash_check() is False:
            # i.e. no exact comparison required
            ligands_by_hash = get_all_ligands_by_graph_hashes(list(self.db.values()))
            grouped_unique_ligands = [ligand_list for ligand_list in ligands_by_hash.values()]

        else:
            print('Exact comparison required, i.e. comparing ligands by isomorphism.')
            grouped_unique_ligands = self.exact_comparison()

        unique_ligands = self.get_unique_ligands_and_set_unique_ligand_name(grouped_unique_ligands)

        # Get unique ligand dictionary with denticity as output format.
        # Note that this is in the new format of the DB
        unique_ligand_dict = {lig.name: lig for lig in unique_ligands}

        print(f'Number of unique ligands: {len(unique_ligands)}.')
        return unique_ligand_dict

    # lastly we have some comparison methods
    def get_df_of_all_ligands(self):
        """
        Returns a dataframe with name, denticity, CSD code and type of every ligand in the database.
        """
        ligand_props = []
        for lig in self.db.values():
            ligand_props.append({
                'name': lig.name,
                'unique_name': lig.unique_name,
                'csd_code': lig.global_props["CSD_code"] if "CSD_code" in lig.global_props else np.nan,
                'original_metal_symbol': lig.original_metal_symbol if hasattr(lig, 'original_metal_symbol') else np.nan,
                'denticity': lig.denticity,
                'graph_hash': lig.graph_hash,
                'coordinates': str(lig.coordinates),
                'atomic_props': str(lig.atomic_props),
                'n_total_unique_ligands': lig.n_total_unique_ligands,
                # 'hash': l.hash # makes issues in pd.testing.assert_frame_equal, probably because of overflow
            })
        df = pd.DataFrame(ligand_props)
        return df

    def calc_number_of_possible_complexes(self, metals: list[str] = None) -> pd.DataFrame:
        if metals is None:
            metals = ['Cr', 'Mn', 'Fe', 'Ru', 'Co', 'Ni']

        df = []
        for metal in metals:
            df_n_metal_combs = self.calc_number_of_possible_complexes_for_metal(metal)
            df.append(df_n_metal_combs)
        df = pd.concat(df, ignore_index=True)

        return df



    def calc_number_of_possible_complexes_for_metal(self, metal: str, geometries: dict = None) -> pd.DataFrame:
        metal_oxi_states = Pymatgen_Element(metal).common_oxidation_states
        results = []

        # possible geometries for octahedral and square-planar complexes. This list needs to be expanded when adding new geometries.
        if geometries is None:
            geometries = {
                'octahedral': [(3, 2, 1), (4, 1, 1), (5, 1), (3, 3)],
                'square_planar': [(2, 2), (2, 1, 1)]
            }

        for oxi_state in metal_oxi_states:
            target_charge = -oxi_state
            for geometry_name, geometry_list in geometries.items():
                for geometry in geometry_list:
                    count = self.calc_number_of_combinations_of_ligands_for_topology(target_charge=target_charge, geometry=geometry)
                    results.append(
                        {'metal': metal, 'oxi_state': oxi_state, 'geometry': geometry_name, 'denticities': geometry,
                         'count': count})

        return pd.DataFrame(results)

    def calc_number_of_combinations_of_ligands_for_topology(self, target_charge: int, geometry: tuple) -> int:
        """
        Calculates the number of possible ligand combinations for a given target charge and geometry.
        @target_charge: The targeted sum of charges of the ligands.
        @geometry: The topology of the complex, e.g. (3, 2, 1) for a octahedral complex with 3 bidentate, 2 monodentate and 1 tridentate ligand.
        """
        n_ligands = len(geometry)
        geometry = sorted(geometry)

        df = self.reduced_df.query('denticity in @geometry and not pred_charge.isnull()')[['denticity', 'pred_charge']].astype(int)

        df = df.groupby(['denticity', 'pred_charge']).size().reset_index().rename(columns={0: 'count'})
        count = 0
        for ligs in itertools.combinations_with_replacement(list(df.itertuples()), n_ligands):
            correct_denticities = sorted(lig.denticity for lig in ligs) == geometry
            correct_charges = sum(lig.pred_charge for lig in ligs) == target_charge
            if correct_charges and correct_denticities:
                # Group ligands which have the same charge and denticity.
                groups = pd.DataFrame(ligs).groupby(['denticity', 'pred_charge'])['count']
                # Calculate the number of possible combinations for this combination of ligands. If there are multiple ligands with the same charge and denticity, we need to pay attention that we don't count the same combination twice, e.g (lig1, lig2) and (lig2, lig1). That is because we define a complex here just in terms of its set of ligands, without caring about the order of the ligands.
                comb_count = 1
                for _, group in groups:
                    n_same_ligands = len(group)
                    lig_count = group.values[0]
                    if n_same_ligands == 1:
                        comb_count *= lig_count
                    else:
                        # Multiple ligands with same charge and denticity: Avoid double counting.
                        comb_count *= comb(lig_count + n_same_ligands - 1, n_same_ligands, exact=True)
                count += comb_count

        return count


class ComplexDB(MoleculeDB):
    type = 'Complex'
    def __init__(self, dict_):
        super().__init__(dict_=dict_)

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


if __name__ == '__main__':

    db_version = '1.7'
    db_path = f'../data/final_db_versions/unique_ligand_db_v{db_version}.json'
    n_max = None

    db = LigandDB.from_json(json_=db_path, type_='Ligand', max_number=n_max)
    df_metals = db.calc_number_of_possible_complexes()
    print(df_metals)

import json
from copy import deepcopy
from tqdm import tqdm
import pandas as pd
import numpy as np
from src01.Molecule import RCA_Molecule, RCA_Complex, RCA_Ligand    # important, don't delete, read in by string that's why they appear grey
from src01.utilities_graph import remove_node_features_from_graph, make_multigraph_to_graph, remove_edge_features_from_graph
from src01.utilities import identify_metal_in_ase_mol
from src01.utilities_Molecule import get_all_ligands_by_graph_hashes, group_list_without_hashing
import networkx as nx
from src01.io_custom import save_json
from datetime import datetime


class BaselineDB:
    def __init__(self, dict_: dict):
        """
        Is basically just a dict with extended functionality
        :param dict_: is dict as {identifier: RCA_Mol/RCA_Ligand}
        So the base functionality will assume that we have an RCA_Mol, which assures
        that the functionality is stable even if we use RCA_Ligands (As they are a subclass of RCA_Mol)
        """
        self.db = dict_

    def __len__(self):
        return len(self.db)

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

    def to_json(self, path, desc: str='Save DB to json'):
        save_json(self.get_dict_in_json_format(desc=desc), path=path, indent=4)

        return

    @classmethod
    def from_json(cls,
                  json_,
                  type_: str = "Molecule",
                  graph_strategy: str = "default",
                  max_number = None,
                  identifier_list: list = None,
                  **kwargs
                  ):
        """

        :param json_: Either the dict itself or the path to a json file
        :param type_: If we want either a molecule or a Ligand DB
        :param graph_strategy: How we want the graphs to be created (this is only important for molecules,
            because ligand graphs are created by the molecule graphs. For Ligands this will just be dumped as a kwarg)
        :param kwargs: additional arguments for the graph creation (only, if no graphs are present in the .json)
        """
        if isinstance(json_, str):
            if not json_.endswith(".json"):
                print("select a json file as input!")
                return

            with open(json_, "r") as file:
                json_dict = json.load(file)

        elif isinstance(json_, dict):
            json_dict = json_

        else:
            print("Wrong Input format for json_")
            return

        new_dict_ = {}

        if type_ not in ["Ligand", "Molecule", "Complex"]:
            print("Wrong type chosen, will be set to Molecule")
            type_ = "Molecule"

        # print("Start Establishing DB from .json")

        #
        #
        #
        if (max_number == False) or (max_number is None):
            max_number = len(json_dict)

        for i, (identifier, mol_dict) in tqdm(enumerate(json_dict.items()), desc=f"Build {type_}Database"):

            # For testing by identifier list
            if identifier_list is not None:
                if identifier not in identifier_list:
                    continue

            # For testing by Max number
            if i >= max_number:
                break

            new_dict_[identifier] = globals()[f"RCA_{type_}"].read_from_mol_dict(dict_=mol_dict,
                                                                                 graph_creating_strategy=graph_strategy,
                                                                                 csd_code=identifier,
                                                                                 **kwargs
                                                                                 )

        return cls(new_dict_)

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
    def __init__(self, dict_):
        super().__init__(dict_=dict_)


class LigandDB(BaselineDB):
    def __init__(self, dict_):
        super().__init__(dict_=dict_)

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


class ComplexDB(MoleculeDB):
    def __init__(self, dict_):
        super().__init__(dict_=dict_)
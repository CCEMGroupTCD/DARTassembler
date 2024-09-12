import pandas as pd
from pathlib import Path
import json
import os
import glob
from tqdm import tqdm as loading_bar
import numpy as np
import io
import networkx as nx

from DARTassembler.src.ligand_extraction.io_custom import load_json


class DataLoader:
    """
    The idea of this class is to check the input folder structure in the init
    of the dataload.

    And to provide the data in the format for the MoleculeDB
    """

    def __init__(self,
                 database_path_,
                 read_in_graphs_if_exist: bool = True,
                 identifier: str = 'CSD_code',
                 **kwargs
                 ):

        self.path = Path(database_path_)
        self.graph_path = Path(self.path, "graphs")
        self.global_path = Path(self.path, "global_mol_properties.csv")
        self.identifier = identifier
        self.data_for_molDB = {}

        self.has_global = None
        self.has_graphs = None
        self.read_graphs_if_possible = read_in_graphs_if_exist

        self.check_folder_structure()

        self.create_atomic_properties_json(**kwargs)

        self.create_MolDB_dict()

    def check_folder_structure(self):
        """
        We check if the .xyzs were provided in an accurate format
        """
        print("Checking Folder Structure")

        if not os.path.exists(Path(self.path, "atomic_properties")):
            print("No atomic Properties Folder found, no Dataload possible")
            raise FileNotFoundError

        atomic_prop_folder = os.listdir(Path(self.path, "atomic_properties"))
        if not any(str(file).endswith(".xyz") for file in atomic_prop_folder):
            print("No .xyz file found")
            raise FileNotFoundError

        if os.path.exists(self.global_path):
            print("Global Properties found")
            self.has_global = True
        else:
            print("No global properties found")
            self.has_global = False

        if os.path.exists(self.graph_path):
            print("Graphs found")
            self.has_graphs = True
        else:
            print("No graphs found")
            self.has_graphs = False

    def get_global_prop_dict(self):
        # refactor?: maybe use df.to_dict() instead? But, then the None values in the csv are converted to NaN, e.g. for the oxidation state, and this changes things when tested for 5000 complexes. global_property_dict = df.set_index(self.identifier, drop=False).to_dict(orient='index')
        df = pd.read_csv(self.global_path)
        json_dict = json.loads(df.set_index(self.identifier).to_json(orient='split'))

        global_property_dict = {}
        for i, csd_code in enumerate(json_dict["index"]):
            global_property_dict[csd_code] = {prop: json_dict["data"][i][j] for j, prop in
                                              enumerate(json_dict["columns"])}
            global_property_dict[csd_code]["CSD_code"] = csd_code

        return global_property_dict

    # todo: Irgendwie schauen, dass ich das beim einlesen der tmqmG mal teste, vorher ist das eh egal
    def get_graph_dict(self) -> dict:

        graph_dir_list = os.listdir(self.graph_path)
        graph_dict = {}

        for path in loading_bar(graph_dir_list, desc="Loading Graphs from .gml files"):
            if path.endswith("gml"):
                identifier = path.removesuffix(".gml")
                graph_dict[identifier] = nx.read_gml(Path(self.graph_path, path))

        return graph_dict

    def create_MolDB_dict(self):
        """
        From the desired input folder structure we are now going to read in the properties
        in the format, which is required to generate a MoleculeDB
        """
        print("Converting DataFolders into dict for MolDB")

        atomic_props_json = load_json(Path(self.path, "atomic_properties", "atomic_properties.json"))

        if self.has_global:
            global_property_dict = self.get_global_prop_dict()
        else:
            global_property_dict = {}

        if self.has_graphs and self.read_graphs_if_possible:
            graph_dict_dict = self.get_graph_dict()
        else:
            graph_dict_dict = {}

        for key_, item_ in atomic_props_json.items():
            global_props_mol = global_property_dict[key_] if key_ in global_property_dict else None
            graph_dict = graph_dict_dict[key_] if key_ in graph_dict_dict else None

            self.data_for_molDB[key_] = {"atomic_props": item_,
                                         "global_props": global_props_mol,
                                         "graph_dict": graph_dict
                                         }

    # copy-pase von ehemaligen MoleculeDatabase(Timo)
    def read_full_atomic_properties_file(self, path, sep='\n\n', file_nr: int = 0, total_files: int = 0):
        with open(path, 'r') as file:
            full_file = file.read()
        files = full_file.split(sep)

        atomic_props = {}
        for filestr in loading_bar(files, desc=f'Reading atomic properties file {file_nr+1} of {total_files}'):
            mol_id, atoms, values, comment = self.read_single_atomic_properties_filestring(filestr)
            assert not mol_id in atomic_props, f'Molecular id {mol_id} not unique in file{path}.'
            values['atoms'] = atoms
            values['comment'] = comment
            atomic_props[mol_id] = values

        return atomic_props

    # copy-pase von ehemaligen MoleculeDatabase(Timo)
    def read_single_atomic_properties_filestring(self, filestr: str):
        num_atoms, mol_id, col_names, comment = self.read_atomic_properties_file_header(filestr)

        file = io.StringIO(filestr)
        txt_array = np.genfromtxt(file, skip_header=2, dtype=str)
        atoms, value_array = txt_array[:, 0], txt_array[:, 1:].astype(float)

        atoms = [str(a) for a in atoms]  # numpy to regular python list of strings for json compatibility
        col_names.pop(0)

        values = {name: column.tolist() for name, column in zip(col_names, value_array.T)}

        assert num_atoms == len(value_array) and num_atoms == len(
            atoms), 'Number of atoms specified in atomic properties file is not equal to the real number of atoms included.'
        return mol_id, atoms, values, comment

    # copy-pase von ehemaligen MoleculeDatabase(Timo)
    @staticmethod
    def read_atomic_properties_file_header(filestr: str, ATOMIC_PROPERTIES_COLUMN_SEPARATOR='  ===  '):
        lines = filestr.split('\n')
        num_atoms = int(lines[0])

        columns = lines[1].split(ATOMIC_PROPERTIES_COLUMN_SEPARATOR)
        mol_id = columns[0]
        comment = columns[-1]
        col_names = columns[1:-1]

        return num_atoms, mol_id, col_names, comment

    def create_atomic_properties_json(self, overwrite=False):
        """
        kawrgs zB
        ATOMIC_PROPERTIES_COLUMN_SEPARATOR, falls man eignen moechte
        """
        json_safe_path = Path(self.path, "atomic_properties", "atomic_properties.json")
        if os.path.exists(json_safe_path) and not overwrite:
            print("atomic_properties.json already exists and no overwriting selected")
            return

        pattern = Path(self.path, "atomic_properties", '*.xyz')
        all_atomic_props_paths = sorted(glob.glob(str(pattern)))

        all_atomic_props = {}

        for i, atm_path in enumerate(all_atomic_props_paths):
            atomic_props = self.read_full_atomic_properties_file(path=atm_path, file_nr=i, total_files=len(all_atomic_props_paths))

            assert not any([mol_id in all_atomic_props for mol_id in
                            atomic_props.keys()]), 'Molecular ids of atomic property files not unique.'
            all_atomic_props.update(atomic_props)

        #
        # and now the saving process
        print('Start saving atomic properties to json.')

        with open(json_safe_path, 'w') as file:
            json.dump(all_atomic_props, file)

        print(f'Saved atomic properties to {json_safe_path}.')

"""
this script should in theory never be shared, because I read in our customized .mol2 files and convert the 339k of them,
which are of interest for us, into a graph and dump that into a graphs.json
The graphs.json shall then lateron serve as a constant quantity and represent the graphs of the CSD
"""
import os
from src01.utilities_graph import mol2_to_graph, mol2_str_to_graph, view_graph, graph_to_dict_with_node_labels
from tqdm import tqdm
import json

import warnings
warnings.filterwarnings("ignore")


class GraphsFromMol2:
    """
    This is highly customized for the input of the CSD database mol2 files
    and not to be considered a general class
    """
    def __init__(self, dir_path):

        with open("../identifierLists/monometallic_identifier_list.json", "r") as file:
            self.monometallic_csd_entries = json.load(file)

        self.dir_path = dir_path
        self.graph_dict = {}

    @staticmethod
    def get_idenifier_from_mol2_str(mol2_str: str):
        """
        extract the csd code from mol2 files
        """
        b = mol2_str.split("\n")
        for i, line in enumerate(b):
            if line == "@<TRIPOS>MOLECULE":
                return b[i + 1]

    def extract_graphs(self):

        for j, path in enumerate(os.listdir(self.dir_path)):

            with open(f"mol2s/{path}", "r") as file:
                content = file.read()

                file_strs = [str_.strip("\n\n\n").strip("\n") for str_ in content.split("\n\n\n\n")]

                for file_str in tqdm(file_strs, desc=f"Reading file {j} of {len(os.listdir('../mol2s'))}"):

                    identifier = self.get_idenifier_from_mol2_str(file_str)

                    if identifier in self.monometallic_csd_entries:
                        #
                        # view_graph(G)
                        try:
                            # put the graph in desired dict format into the graph dict
                            G = mol2_str_to_graph(str_=file_str)
                            self.graph_dict.update({
                                identifier : graph_to_dict_with_node_labels(G)
                            })
                        except Exception as e:
                            # error_count += 1
                            # print(f"An error has occured: {e}\nerror count: {error_count}")
                            pass

    def safe_graph_dict(self, safe_path: str = None):

        if safe_path is None:
            # will put into the same directory
            safe_path = "CSD_Graphs.json"

        with open(safe_path, "w+") as file:
            json.dump(self.graph_dict, file)
            # dump graphs to json

    def filter_graphs(self, id_list, safe_path: str):
        """
        e.g. to get the graphs for the csd
        """

        new_dict_of_graphs = {}

        for id_, d in self.graph_dict.items():
            if id_ in id_list:
                new_dict_of_graphs[id_] = d

        with open(safe_path, "w+") as file:
            json.dump(new_dict_of_graphs, file)


if __name__ == "__main__":

    Graphs = GraphsFromMol2(dir_path="../mol2s")

    Graphs.extract_graphs()

    Graphs.safe_graph_dict()

"""
Test, which Graph creating methods performs best
"""
import networkx as nx
from src01.DataBase import MoleculeDB, LigandDB
from src01.DataLoader import DataLoader
from src01.constants import metals_in_pse
from pymatgen.core.periodic_table import Element as Pymatgen_Element
import warnings
import numpy as np

import gc

warnings.filterwarnings("ignore")


class GraphTesting:

    def __init__(self,
                 database_path: str = None,
                 run_list: dict = None
                 ):

        self.tmQM_unique_Ligands = None
        self.tmQM_Ligands = None
        self.tmQM_DB = None
        if database_path is None:
            # defaultpath
            self.database_path = '../../database/tmQM'
        else:
            self.database_path = database_path

        self.colorlist = ["r", "b", "g", "yellow", "orange", "black"]

        self.run_list = run_list

        self.results = {}

    def run(self,
            metrics: list[str],
            Testing: bool = True):

        # create empty plots
        self.results = {metric: {} for metric in metrics}

        for identifier, d in self.run_list.items():

            gc.collect()

            # 1. Create Database wih graphs according to run specifications
            if Testing is True:
                self.tmQM_DB = MoleculeDB.from_json(json_=DataLoader(database_path_=self.database_path).data_for_molDB,
                                                    type_="Molecule",
                                                    max_number=1000,
                                                    graph_strategy=d["strat"],
                                                    **d["kwargs"]
                                                    )
            else:
                self.tmQM_DB = MoleculeDB.from_json(json_=DataLoader(database_path_=self.database_path).data_for_molDB,
                                                    type_="Molecule",
                                                    graph_strategy=d["strat"],
                                                    **d["kwargs"]
                                                    )

            # 2. Create Ligand database
            # Note: It is crucial to set the correct numbers of denticities, because otherwise -1 denticitated
            #       ligands are getting filtered out by default
            self.tmQM_Ligands = LigandDB.from_MoleculeDB(molDB=self.tmQM_DB,
                                                         denticity_numbers_of_interest=[-1, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                                                        10]
                                                         )

            # 3. Extract unique Ligands
            self.tmQM_unique_Ligands = LigandDB(self.tmQM_Ligands.filter_duplicates())

            self.evaluate_metrics(metrics, identifier)

    def evaluate_metrics(self, metrics, run_identifier):

        for i, metric in enumerate(metrics):
            self.results[metric][run_identifier] = self.evaluate(name=metric)

    def evaluate(self, name: str):
        """
        we return a value for each metric given a database according to their name
        """

        if name == "number unqiue ligands":
            return len(self.tmQM_unique_Ligands.db.values())
        elif name == "number of not fully connected graphs":
            return len([mol for mol in self.tmQM_DB.db.values() if nx.is_connected(mol.graph) is False])
        elif name == "number of isolated ligands":
            return len([lig for lig in self.tmQM_Ligands.db.values() if lig.denticity == -1])
        else:
            warnings.warn("Metric not implemented")
            pass

    def show_plots(self):
        import matplotlib
        matplotlib.use('TkAgg', force=True)
        import matplotlib.pyplot as plt

        for i, (metric, res) in enumerate(self.results.items()):
            plt.plot(
                [self.run_list[k]["name"] for k in res.keys()],
                res.values(),
                "o",
                alpha=0.5,
                color=self.colorlist[i]
            )
            plt.title(metric)
            plt.xticks(rotation=45)
            plt.show()


if __name__ == "__main__":
    # Setting up the run list
    #
    run_list = {
        0: {
            "name": "nat cutoffs",
            "strat": "default",
            "kwargs": {}
        },
    }

    """
    for corr in np.linspace(-0.1, -1, 10):
        run_list.update(
            {len(run_list): {
                "name": f"corrected nat cutoffs by {corr} A",
                "strat": "default",
                "kwargs": {
                    "skin_": 0.2,
                    "cutoff_corrections_for_metals": {i: corr for i in metals_in_pse}
                }
            }
            }
        )

    for corr in np.linspace(-0.1, -2, 20):
        run_list.update(
            {len(run_list): {
                "name": f"corrected nat cutoffs by {corr} A, skin = 0.3",
                "strat": "default",
                "kwargs": {
                    "skin_": 0.3,
                    "cutoff_corrections_for_metals": {i: corr for i in metals_in_pse}
                }
            }
            }
        )
    """

    run_list.update(
        {len(run_list): {
            "name": f"MolSimplify",
            "strat": "molsimplifyGraphs",
            "kwargs": {}
        }
        }
    )

    '''
    run_list.update(
        {len(run_list): {
            "name": "Pymatgen NN",
            "strat": "pymatgen_NN",
            "kwargs": {}
        }
        }
    )
    '''

    run_list.update(
        {len(run_list): {
            "name": "CSD Graphs",
            "strat": "CSD",
            "kwargs": {}
        }
        }
    )

    #
    #
    # Actually running the comparison
    GT = GraphTesting(run_list=run_list)

    GT.run(metrics=[
        "number unqiue ligands",
        "number of not fully connected graphs",
        "number of isolated ligands"
    ],
        Testing=True
    )

    GT.show_plots()

    print("Done")

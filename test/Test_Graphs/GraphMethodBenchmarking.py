"""
Building up on the new main_ligand_extraction.py,
is going to become the graph testing thing lateron
"""
import networkx as nx
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np

from src01.main_ligand_extraction import main as run
from src01.main_ligand_extraction import select_example_database
from constants.constants import project_path, metals_in_pse


def generate_run_list():
    run_list = {
        0: {
            "name": "nat cutoffs",
            "strat": "default",
            "kwargs": {}
        },
    }

    #
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

    run_list.update(
        {len(run_list): {
            "name": f"MolSimplify",
            "strat": "molsimplifyGraphs",
            "kwargs": {}
        }
        }
    )

    run_list.update(
        {len(run_list): {
            "name": "Pymatgen NN",
            "strat": "pymatgen_NN",
            "kwargs": {}
        }
        }
    )

    run_list.update(
        {len(run_list): {
            "name": "CSD Graphs",
            "strat": "CSD",
            "kwargs": {}
        }
        }
    )

    return run_list


def print_exact_comparison(ex: dict):

    print(f"For the run we obtained a rate of {ex['positive_rate']} correct ligand denticites \n"
          f"while evaluating {ex['total_evaluations']} ligands from the groundtruth")

    return


class GraphTester:

    def __init__(self):

        self.groundtruth = self.prepare_ground_truth()

    @staticmethod
    def prepare_ground_truth(
            test_file_path="test/debug/databases/charge_benchmark/all_ligand_charges_with_high_confidence_v1.3.csv"
    ):

        path = Path(project_path, test_file_path)

        df = pd.read_csv(path)

        return df[["CSD_code", "stoichiometry", "denticity"]].set_index(["CSD_code", "stoichiometry"]).to_dict(orient="index")

    @staticmethod
    def csd_code_from_ligand_name(ligand_name):
        return ligand_name.split("-")[1]

    def run(self,
            graph_strategy,
            selected_DB,
            testing,
            **kwargs
            ):

        if selected_DB.lower() in ["tmqm", "tmqmg"] and graph_strategy == "CSD":
            raise ValueError("CSD Graphs can not be used with tmqm atomic properties")

        #
        database_path, data_store_path = select_example_database(DB=selected_DB)
        #
        db = run(
            database_path_=database_path,
            data_store_path_=data_store_path,
            calculate_charges_=False,
            overwrite_atomic_properties_=False,
            use_existing_input_json_=False,
            exclude_not_fully_connected_complexes_=False,
            get_only_unique_ligand_db_without_charges_=True,
            testing_=testing,
            graph_strat_=graph_strategy,
            **kwargs
        )

        # Complex_DB = db.complex_db: MoleculeDB object
        # Ligand_DB = db.full_ligand_db
        # Full_Ligand_DB = db.unique_ligand_db

        # compute the heuristics
        # which are less important for now at least - only expressive if we make a full run.
        heuristics = {
            "number_of_unconnected_complexes": len([mol for mol in db.complex_db.db.values() if nx.is_connected(mol.graph) is False]),
            "number_of_total_ligands": len(db.full_ligand_db.db),
            "number_of_isolated_ligands": len([lig for lig in db.full_ligand_db.db.values() if lig.denticity == -1]),
            "number_un_ligands_w_denticity": sum([len(ul.count_denticities) for ul in db.unique_ligand_db.db.values()])
        }

        positive, negative = 0, 0

        for name, lig in tqdm(db.full_ligand_db.db.items(), desc="Comparing to groundtruth denticities"):

            csd_code = self.csd_code_from_ligand_name(name)
            stoich = lig.stoichiometry

            predicted_dent_by_graph = lig.denticity

            if (csd_code, stoich) in self.groundtruth:

                if float(predicted_dent_by_graph) == float(self.groundtruth[(csd_code, stoich)]["denticity"]):
                    positive += 1
                else:
                    negative += 1

        exact_denticty_comparison = {
            "positive_rate": positive/(negative+positive),
            "total_evaluations": negative+positive
        }

        return heuristics, exact_denticty_comparison


if __name__ == "__main__":

    GT = GraphTester()

    id_list = list(set([key[0] for key in GT.groundtruth]))

    df = pd.DataFrame(columns=["graph_strat", "rate", "evaluated_entries"])

    runlist = generate_run_list()

    for run_ in runlist.values():

        h, e = GT.run(graph_strategy=run_["strat"],
                      selected_DB="CSD_MM_G",
                      testing=id_list,
                      **run_["kwargs"]
                      )

        df.loc[len(df.index)] = [run_["name"], e["positive_rate"], e["total_evaluations"]]
        # print_exact_comparison(e)

    print(df)
    #
    print("done")

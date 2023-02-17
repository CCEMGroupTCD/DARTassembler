from typing import Union
from copy import deepcopy
import gc
import yaml
from yaml import SafeLoader
from pathlib import Path
from constants.constants import project_path

from src01.main_ligand_extraction import main as run_ligand_extraction
from src01.main_ligand_extraction import select_example_database

from src02_Pre_Assembly_Filtering.Filter import Filter

from src04_Assembly.main import Assembly


class Run:
    """
    This shall facilitate a full run
    """

    def __init__(self,
                 database: str,
                 testing: Union[bool, int] = False,
                 graph_stragy: str = "default",
                 calculate_charges: bool = False,
                 overwrite_atomic_properties: bool = False,
                 use_existing_input_json: bool = False,
                 exclude_not_fully_connected_complexes: bool = False,
                 get_charges: bool = False,
                 filter_setup_path: str = "Filter_setup.yml",
                 batch_yml_path: str = "example_batch.yml",
                 optimisation: bool = True,
                 random_seed: int = 1,
                 gaussian_input_specifications: str = "gaussian_specifications.yml"
                 ):
        """

        testing = False  # if we would like to only do a test run (only works from the second run on)
        graph_strategy = "default"  # the desired graph strategy: default, ase_cutoff, CSD, pymatgen_NN, molsimplifyGraphs
        calculate_charges = True  # if you want to run charge assignment after ligand extraction, takes ~30 min on tmQMg
        overwrite_atomic_properties = True  # if atomic properties json should be overwritten, not really critical
        use_existing_input_json = False  # if the existing input json should be used or the process started from the xzy files
        exclude_not_fully_connected_complexes = False  # only keep complexes which are fully connected
        get_only_unique_ligand_db_without_charges = False  # For graph benchmark useful, reduces runtime because it ignores charge assignment and updating the complex and full ligand db.
        """
        self.database_name = database
        self.testing = testing
        self.graph_strat = graph_stragy
        self.database_path, self.data_store_path = select_example_database(DB=database)
        self.calculate_charges = calculate_charges
        self.overwrite_atomic_properties = overwrite_atomic_properties
        self.use_existing_input_json = use_existing_input_json
        self.exclude_not_fully_connected_complexes = exclude_not_fully_connected_complexes
        # get rid of the weird double negation
        self.get_only_unique_ligand_db_without_charges = not get_charges

        # extraction
        self.db = self.ligand_extraction()
        self.unique_ligand_DB = deepcopy(self.db.unique_ligand_db)

        # free some memory
        del [self.db]
        gc.collect()

        # pre assembly filtering
        self.filter_setup_path = filter_setup_path
        self.assembly_ready_ligands = self.pre_assembly_filtering()

        # assembly
        # set output-path to default and make sure it exists
        self.assembly_output_path = Path(self.data_store_path, "Assembled_Molecules/")
        self.assembly_output_path.mkdir(parents=True, exist_ok=True)

        # run the assembly
        self.assembled_complexes = self.assembly(
            batch_yml_path=batch_yml_path,
            random_seed=random_seed,
            optimisation=optimisation
        )

        self.store_assembled_complexes(
            gaussian_input_specifications=gaussian_input_specifications
        )

        print(f"Run sucessfull: Complexes assembled and stored to {self.assembly_output_path}")

    def ligand_extraction(self):

        return run_ligand_extraction(
            database_path_=self.database_path,
            data_store_path_=self.data_store_path,
            calculate_charges_=self.calculate_charges,
            overwrite_atomic_properties_=self.overwrite_atomic_properties,
            use_existing_input_json_=self.use_existing_input_json,
            exclude_not_fully_connected_complexes_=self.exclude_not_fully_connected_complexes,
            get_only_unique_ligand_db_without_charges_=self.get_only_unique_ligand_db_without_charges,
            testing_=self.testing,
            graph_strat_=self.graph_strat
        )

    def pre_assembly_filtering(self):

        #
        F = Filter(self.data_store_path,
                   filter_yaml_path=self.filter_setup_path,
                   ligand_db=self.unique_ligand_DB
                   )

        F.run_filters()
        F.add_constant_ligands()
        F.add_reactant()
        F.safe()

        return F.filtered_db

    def assembly(self,
                 batch_yml_path: str,
                 optimisation: bool = False,
                 random_seed: int = 1
                 ):

        with open(batch_yml_path, "r") as file:
            batch_inf = yaml.load(file, SafeLoader)

        batch_dict = {
            'Name': batch_inf["Name"],
            "Input_Path": str(Path(self.data_store_path, "Filtering", "filteredLigDB.json")),
            'Output_Path': str(self.assembly_output_path),
            "MAX_num_complexes": batch_inf['MAX_num_complexes'],
            "Isomers": batch_inf["Isomers"],
            "Optimisation_Choice": optimisation,
            "Random_Seed": random_seed
        }

        for i, top in enumerate(batch_inf["Topologies"]):
            batch_dict[f"Topology_{i}"] = top

        for j, metal in enumerate(batch_inf["Metals"]):
            batch_dict[f"Metal_{j}"] = metal

        batch_list = [batch_dict]

        A = Assembly(batch_list)

        c = A.assembly_main()

        return c

    def store_assembled_complexes(self,
                                  gaussian_input_specifications: str
                                  ):

        with open(gaussian_input_specifications, "r") as file:
            gis = yaml.load(file, SafeLoader)

        basis_set_dict = {key: "\n".join(value) for key, value in gis.items()}

        for tmc in self.assembled_complexes:

            tmc.to_json(path=Path(self.assembly_output_path, f"{tmc.name}.json"))

            tmc.to_com(path=Path(self.assembly_output_path, f"{tmc.name}.com"),
                       basis_set_dict=basis_set_dict
                       )

        return




# Standard library imports
from typing import Union, Dict, Any, List, Tuple, Optional

import ase.io
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
import datetime
import logging
import random
import sys

# DART specific imports
from DARTassembler.src.assembler.isomer import AssembledIsomer, AssembledComplex
from DARTassembler.src.assembler.output import AssemblerOutput, BatchAssemblerOutput, ComplexAssemblerOutput
from DARTassembler.src.metalig.utils_molecule import get_standardized_stoichiometry_from_atoms_list
from DARTassembler.src.constants.paths import default_assembler_yml_path
from DARTassembler.src.assembler.ligands import LigandChoice
from DARTassembler.src.modules.modules import BaseModule
from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.misc.io import read_yaml, get_correct_ligand_db_path_from_input, save_json

# Data processing imports
from pathlib import Path
import pandas as pd
import numpy as np

# Warnings
import warnings
warnings.simplefilter("always")

class Assembler(BaseModule):
    """
    This module assembles all isomers of complexes from specified ligand databases and metal centers. It saves all assembled complexes to the output directory as .xyz and .json files, together with a csv file containing information about the assembly run.
    """

    def __init__(self,
                 output_directory: Union[str, Path] = 'DARTassembler',
                 verbosity: int = 2,
                 complex_name_length: int = 8,
                 n_max_ligands: Optional[int] = None,
                 ):
        """
        The main class for the DART workflow. It assembles all isomers of complexes from specified ligand databases and metal centers. Finally, all assembled complexes are saved to the output directory as .xyz and .json files, together with a csv file containing information about the assembly run.
        :param output_directory: Directory where the output files are saved. If it does not exist, it is created.
        :param verbosity: If 0, only errors are printed. If 1, warnings are printed. If 2, info messages are printed. If 3, debug messages are printed.
        :param complex_name_length: The length of the randomly generated complex name. If a duplicate name is generated, the length is increased by one character until a unique name is found.
        :param n_max_ligands: Maximum number of ligands to load from the ligand databases provided in each batch. If None, all ligands are loaded. Useful to speed up test runs.
        """
        super().__init__()
        self.output_directory = Path(output_directory).resolve()
        self.verbosity = verbosity
        self.complex_name_length = complex_name_length
        self.n_max_ligands = n_max_ligands

        self._loaded_ligand_databases = {}  # to avoid reloading the same ligand database multiple times

        # Keep track of the input arguments
        self.init_args = {**locals()}
        self.init_args.pop('self')
        self.init_args.pop('__class__')

        # Set up the output directories
        self.gbl_outcontrol = AssemblerOutput(outdir=self.output_directory)

        # Set up logging
        verbosity2logging = {0: logging.ERROR, 1: logging.WARNING, 2: logging.INFO, 3: logging.DEBUG}
        stream_handler = logging.StreamHandler(stream=sys.stdout)
        stream_handler.setLevel(verbosity2logging[self.verbosity])
        # Print to stdout
        logging.basicConfig(level=verbosity2logging[self.verbosity], format='%(message)s', handlers=[logging.FileHandler(self.gbl_outcontrol.log_path, mode='w'), stream_handler])

    def run_batches(self, batches: list[dict]) -> None:
        """
        Runs the whole assembly for all batches specified in the assembly input file.
        :param batches: List of dictionaries with the batch settings. See the documentation for more information.
        :return None
        """
        self.batches = batches
        for idx, batch in enumerate(self.batches):
            batch['name'] = batch.get('name', f'batch_{idx}')  # Use batch name or generate a default one
        self.n_batches = len(self.batches)
        self.df_info = []
        self.all_tried_complex_names = set()    # Keep track of all tried complex names to avoid duplicates
        self.successfully_assembled_isomer_names = []

        self.batches_args = {**locals()}
        self.batches_args.pop('self')
        self.input_args = {**self.init_args, **self.batches_args}
        self._log_global_info()

        self._check_batch_settings(batches)

        # Save yml file with input arguments to output directory
        self.gbl_outcontrol.save_settings(self.input_args)

        start = datetime.datetime.now()
        for idx, batch_settings in enumerate(self.batches):
            self._log_batch_title_and_settings(batch_settings=batch_settings),
            self._run_batch(batch_idx=idx, **batch_settings)  # run the batch assembly

        self.runtime = datetime.datetime.now() - start  # keep track of the runtime to display later
        self._make_and_save_output_csv()
        self._log_summary()
        self._final_checks()

        return

    @staticmethod
    def _check_batch_settings(batches: list) -> None:
        # Check all batch names are unique
        batch_names = [batch['name'] for batch in batches]
        if not len(batch_names) == len(set(batch_names)):
            raise ValueError(f"DART batch names must be unique. The following batch names are not unique: {batch_names}")
        return

    def _run_batch(self,
                   name: str,
                   batch_idx: int,
                   target_vectors: list[list[list[float]]],
                   metal_centers: Union[str, tuple[str, tuple[float, float, float]]],
                   n_max_complexes: Union[int, str],
                   ligand_db_files: Union[list[str], str] = 'metalig',
                   ligand_archetypes: list[str] = None,
                   ligand_origins: list[tuple[float, float, float]] = None,
                   total_ligand_charges: int = None,
                   monoaxial_optimization: bool = True,
                   permutable_ligands: list[int] = None,
                   force_all_isomers: bool = False,
                   duplicate_tolerance: float = 0.5,
                   clashing_tolerance: float = -0.3,
                   clashing_metal: bool = False,
                   complex_name_suffix: str = '',
                   random_seed: Optional[int] = None,
                   background_file: str = None,
                   background_translation: Optional[Tuple[float, float, float]] = None,
                   ) -> None:
        """
        Run the assembly for one batch.
        :param name: Name of the batch.
        :param batch_idx: Index of the batch in the list of batches.
        :param target_vectors: List of target vectors for the complexes to be assembled. Each target vector is a list of lists of floats/ints.
        :param metal_centers: Metal centers to be used for the complexes. Can be a string (e.g. 'Ru') or a tuple of (metal symbol, (x, y, z)).
        :param n_max_complexes: Maximum number of complexes to be assembled in this batch. If 'all', all complexes are assembled.
        :param total_ligand_charges: Choose ligands so that the total charge of all ligands in the complex is equal to this value. If None, any ligand combination is assembled.
        :param ligand_db_files: List of ligand database files to be used for the batch. If None, defaults to ['metalig'] for each set of target vectors. If 'same_as_previous', the same ligand database as in the previous batch is used. If a single string is provided, it is expanded to a list of the same length as target_vectors.
        :param ligand_origins: Coordinates for each ligand, which will be used as the origin of the ligand rotation in the complex. If None, defaults to the center of all metal center coordinates for each ligand.
        :param complex_name_suffix: Suffix to be added to each complex name. Defaults to ''.
        :param random_seed: Random seed for reproducibility. If None, defaults to the batch_idx so that each batch is deterministic but different from each other.
        :return: None
        """
        # Set random seed for reproducibility. Do this batch-wise so every batch is reproducible independently.
        if random_seed is None:
            random_seed = batch_idx
        random.seed(random_seed)
        np.random.seed(random_seed)

        # Handle defaults
        if isinstance(ligand_db_files, str):    # Expand a single path to a list of paths for each ligand.
            ligand_db_files = [ligand_db_files for _ in target_vectors]

        self.batch_name = name
        self.batch_idx = batch_idx
        self.ligand_db_files = ligand_db_files   # do not resolve path here to keep keywords like 'metalig'
        self.n_max_complexes = n_max_complexes
        self.random_seed = random_seed
        self.total_ligand_charges = total_ligand_charges
        self.target_vectors = target_vectors
        self.ligand_origins = ligand_origins
        self.metal_centers = metal_centers
        self.complex_name_suffix = complex_name_suffix
        self.permutable_ligands = permutable_ligands
        self.clashing_tolerance = clashing_tolerance
        self.clashing_metal = clashing_metal
        self.duplicate_tolerance = duplicate_tolerance
        self.monoaxial_optimization = monoaxial_optimization
        self.ligand_archetypes = ligand_archetypes
        self.force_all_isomers = force_all_isomers
        self.batch_output_path = Path(self.gbl_outcontrol.batch_dir, self.batch_name)
        self.batch_outcontrol = BatchAssemblerOutput(self.batch_output_path)
        self.background_file = background_file
        self.background_translation = background_translation

        # Redirect tqdm to the logging module so that messages appear properly on two different lines
        with (logging_redirect_tqdm()):

            # Load the ligand databases and cache them for later use
            self.ligand_dbs = self._get_ligand_databases()

            # Set up an iterator for the ligand combinations
            ligand_choice = LigandChoice(
                                    ligand_dbs=self.ligand_dbs,
                                    total_ligand_charges=self.total_ligand_charges,
                                    n_max_complexes=self.n_max_complexes,
                                    )
            ligand_combinations = ligand_choice.choose_ligands()

            # Set progress bar with or without final number of assembled complexes
            total = self.n_max_complexes if self.n_max_complexes != 'all' else None
            progressbar = tqdm(desc='Assembling complexes', unit=' complexes', file=sys.stdout, total=total, disable=self.verbosity < 2)

            batch_sum_assembled_complexes = 0  # Number of assembled complexes in this batch
            while ligand_choice.if_make_more_complexes(batch_sum_assembled_complexes):

                # Choose ligands for complex
                try:
                    ligands = next(ligand_combinations)
                except StopIteration:
                    break # If all ligand combinations are exhausted, stop the batch

                complex = AssembledComplex(
                    ligands=ligands,
                    target_vectors=self.target_vectors,
                    ligand_origins=self.ligand_origins,
                    metal_centers=self.metal_centers,
                )
                complex.generate_isomers(
                                            clashing_tolerance= self.clashing_tolerance,
                                            clashing_metal= self.clashing_metal,
                                            duplicate_tolerance= self.duplicate_tolerance,
                                            permutable_ligands = self.permutable_ligands,
                                            monoaxial_optimization = self.monoaxial_optimization,
                                            force_all_isomers=self.force_all_isomers,
                                            complex_name_length=self.complex_name_length,
                                            complex_name_suffix=self.complex_name_suffix,
                                            avoid_names=self.all_tried_complex_names,  # Avoid names of already tried complexes
                )
                # Add the complex name to the set of all tried complex names to avoid duplicates in the next iteration
                self.all_tried_complex_names.add(complex.complex_name)

                self._save_assembled_isomers(complex=complex)

                # Update counters if at least one isomer was successfully assembled for this complex
                if complex.success:
                    batch_sum_assembled_complexes += 1
                    progressbar.update(1)

            progressbar.close()

        return

    def _get_ligand_databases(self) -> list[LigandDB]:
        ligand_databases = []
        for idx, (target_vectors, ligand_db_filepath) in enumerate(zip(self.target_vectors, self.ligand_db_files)):
            if ligand_db_filepath == 'same_as_previous':
                if not idx > 0:
                    raise ValueError("The first ligand database cannot be 'same_as_previous'. Please provide a valid ligand database file path.")
                ligand_databases.append('same_as_previous')
                continue

            if not ligand_db_filepath in self._loaded_ligand_databases:
                ligand_db_filepath = get_correct_ligand_db_path_from_input(ligand_db_filepath)
                self._loaded_ligand_databases[ligand_db_filepath] = LigandDB.from_json(ligand_db_filepath, n_max=self.n_max_ligands, show_progress=self.verbosity > 1)

            # Reduce to the ligands which have the correct n_eff_denticity for the specified target vectors
            n_target_vectors = len(target_vectors)
            database = {name: ligand for name, ligand in self._loaded_ligand_databases[ligand_db_filepath].db.items() if ligand.n_eff_denticities == n_target_vectors}

            # If required, reduce the database to the ligands with the correct archetype
            if self.ligand_archetypes is not None:
                archetype = self.ligand_archetypes[idx]
                database = {name: ligand for name, ligand in database.items() if ligand.archetype == archetype}

            if not database:
                with_archetypes = f' and ligand archetype `{archetype}`' if self.ligand_archetypes is not None else ''
                raise ValueError(f"The provided ligand database contains no ligands with `n_eff_denticities={n_target_vectors}`{with_archetypes}. Please check your input ligand database `{Path(ligand_db_filepath).resolve()}`.")

            database = LigandDB(database)  # Convert to LigandDB object
            ligand_databases.append(database)

        return ligand_databases

    def _save_assembled_isomers(self, complex):
        """
        Save the successfully assembled isomers of the complex to the output directory in the following format:
        {complex.complex_name}.json: Complex data in json format with all isomers and their properties.
        {complex.complex_name}/{complex.complex_name}{isomer_idx}.xyz: XYZ file of the isomer with a comment containing the isomer name, warning, and ligand unique names.
        Additionally, the xyz files are saved to the concatenated xyz files (either passed or failed) of this batch and to the all_xyz file.
        """
        # Save complex json file with all isomer data
        if complex.success:
            complex_dir = Path(self.batch_outcontrol.complex_dir, complex.complex_name)
            complex_json_filepath = complex_dir / f'{complex.complex_name}.json'
            data = complex.to_dict()
            save_json(db=data, path=complex_json_filepath, mkdir=True, indent=4)

        for isomer in complex.isomers:
            success = True if isomer.warning == '' else False
            isomer_name = isomer.isomer_name

            # Add a comment to the xyz file with complex and ligand names for easier identification
            xyz_comment = f'isomer_name: {isomer_name}, warning: {isomer.warning}, ligand_unique_names: ({", ".join(isomer.ligand_info["unique_names"])})'
            xyz_string = isomer.get_xyz_string(comment=xyz_comment)

            # Save xyz of isomer to complex directory
            if success:
                isomer_xyz_filepath = complex_dir / f'{isomer_name}.xyz'
                with open(str(isomer_xyz_filepath), 'w') as xyz_file:
                    xyz_file.write(xyz_string)
                self.successfully_assembled_isomer_names.append(isomer_name)
                if self.background_file is not None:
                    # Combine the DART generated xyz with the extra structure xyz into one new xyz file and save it in the same location

                    extra_structure = ase.io.read(self.background_file, format="xyz")

                    # Ensure we always use a NumPy array for the translation vector
                    translation = np.array(self.background_translation or [0.0, 0.0, 0.0], dtype=float)
                    extra_structure.positions += translation

                    # Combine structures (order doesn’t matter unless you care about atom order)
                    combined_structure = extra_structure + isomer.atoms

                    # Write combined XYZ
                    combined_xyz_filepath = complex_dir / f"{isomer_name}_combined.xyz"
                    ase.io.write(str(combined_xyz_filepath), combined_structure, format="xyz")


            # Save to concatenated xyz files of this batch
            self.batch_outcontrol.save_xyz(xyz_string, success=success, append=True)    # passed/failed xyz files are created automatically
            self.batch_outcontrol.save_file(xyz_string, self.batch_outcontrol.all_xyz_path, append=True)        # save to all_xyz file

            # Save data for csv file.
            isomer_idx = (len(self.successfully_assembled_isomer_names) - 1) if success else None
            self._add_batch_info(isomer=isomer, success=success, isomer_idx=isomer_idx, complex_name=complex.complex_name)

        return

    def _add_batch_info(self, isomer: AssembledIsomer, success, isomer_idx: int, complex_name: str) -> None:
        """
        Add information about the batch to the batch info variable which will be saved to the batch info file.
        """
        elements = isomer.get_metal_symbols()
        metal_stoi = get_standardized_stoichiometry_from_atoms_list(elements)
        data = {
            'success': success,
            'isomer_idx': isomer_idx,
            'isomer_name': isomer.isomer_name,
            'complex_name': complex_name,
            'stoichiometry': isomer.stoichiometry,
            'graph_hash': isomer.graph_hash,
            'warning': isomer.warning,
            'ligand_unique_names': isomer.ligand_info['unique_names'],
            'ligand_archetypes': isomer.ligand_info['archetypes'],
            'ligand_stoichiometries': isomer.ligand_info['stoichiometries'],
            'ligand_charges': isomer.ligand_info['charges'],
            'ligand_donors': isomer.ligand_info['donors'],
            'batch_idx': self.batch_idx,
            'batch_name': self.batch_name,
            'metal_centers': metal_stoi,
            'total_ligand_charges': self.total_ligand_charges,
            'random_seed': self.random_seed,
        }
        self.df_info.append(data)

        return

    @classmethod
    def run_from_yaml(cls, input: Union[str, Path, None]) -> 'Assembler':
        """
        Run the assembler from a yaml file. See the documentation for the input format. The Assembler is run and the output is saved to the specified output directory. The Assembler object is returned.
        :param input: Path to the yml file with the input settings.
        :return: Assembler object.
        """
        if input is None:
            input = default_assembler_yml_path

        options = read_yaml(input)
        batches = options.pop('batches')

        assembler = Assembler(**options)
        assembler.run_batches(batches=batches)

        return assembler

    @classmethod
    def run_from_cli(cls, input: Union[str, Path, None]) -> 'Assembler':
        """
        Run the assembler from a yaml file specified in the command line interface. The Assembler is run and the output is saved to the specified output directory. The Assembler object is returned.
        :param input: Path to the yml file with the input settings.
        :return: Assembler object.
        """
        super()._before_run_from_cli()
        super()._print_cli_input(input=input)
        assembler = cls.run_from_yaml(input=input)
        super()._after_run_from_cli()

        return assembler

    def _make_and_save_output_csv(self) -> None:
        """Save output info csv of all attempts."""
        self.df_info = pd.DataFrame(self.df_info)
        self.df_info['attempt'] = self.df_info.index
        self.df_info = self.df_info[['attempt'] + [col for col in self.df_info.columns if col != 'attempt']]  # Move attempt column to front

        outdf = self.df_info.copy()
        # Make lists in the dataframe to strings for saving to csv
        for col in outdf.columns:
            if isinstance(outdf[col].iloc[0], list):
                outdf[col] = outdf[col].apply(lambda x: f'({", ".join(str(el) for el in x)})' if isinstance(x, list) else x)
        self.gbl_outcontrol.save_run_info_table(outdf)

        return

    def _final_checks(self) -> None:
        """Some final checks."""
        df_test_success = self.df_info[self.df_info['success']]
        batches = df_test_success['batch_idx'].unique()
        for batch in batches:
            df_batch = df_test_success[df_test_success['batch_idx'] == batch]
            # Check for duplicate complex names in the batch
            duplicate_names = df_batch['isomer_name'][df_batch['isomer_name'].duplicated()].values
            assert len(duplicate_names) == 0, f"Duplicate isomer names in batch {batch}: {duplicate_names}. Please report this issue to our GitHub page."

        return

    def _log_summary(self) -> None:
        """Log nice summary per batch."""
        batch_summary_title = '  Summary per batch  '
        logging.info(f'{batch_summary_title:=^80}')
        for batch_idx, batch in enumerate(self.batches):
            batch_name = batch['name']
            df = self.df_info[self.df_info['batch_idx'] == batch_idx]
            logging.info(f"{batch_name}:")
            self._log_success_rate(df)

        # Print total summary of run
        total_summary_title = '  Total summary of DART Assembler run  '
        logging.info(f'{total_summary_title:=^80}')
        self._log_success_rate(self.df_info)
        n_isomers = self.df_info['success'].sum()
        n_complexes = self.df_info[self.df_info['success']]['graph_hash'].nunique()
        logging.info(f"DART Assembler output files saved to directory `{self.output_directory.name}`.")

        # The runtime is printed but not logged, so that slight differences in the runtime do not cause the integration tests to fail.
        if self.verbosity > 1:
            print(f"Total runtime for assembling {n_isomers} isomers (from {n_complexes} complexes): {self.runtime}")

        return

    def _log_global_info(self) -> None:
        """Log global information about the run."""
        logging.info('Starting DART Assembler Module.')
        logging.info(f'Output directory: {self.output_directory.name}')
        plural = 'es' if self.n_batches > 1 else ''  # print plural or singular in next line
        logging.info(f"Running {self.n_batches} batch{plural}...")
        logging.info(f"User-defined global settings:")
        for key, value in self.init_args.items():
            logging.info(f"    {key: <30}{value}")

        return

    @staticmethod
    def _log_success_rate(df):
        """Log success rate of the run."""
        n_total = len(df)
        n_isomers = df['success'].sum()
        n_complexes = df[df['success']]['complex_name'].nunique()

        # Output statistics how many isomers failed each filter
        post_filters = df['warning'].value_counts().to_dict()
        # Merge all warnings that start with 'duplicate' into one
        n_duplicates = sum(count for note, count in post_filters.items() if note.startswith('duplicate'))
        post_filters = {note: count for note, count in post_filters.items() if not note.startswith('duplicate')}
        if n_duplicates > 0:
            post_filters['duplicate'] = n_duplicates
        # Sort the post-filters by the number of occurrences
        post_filters = dict(sorted(post_filters.items(), key=lambda item: item[1], reverse=True))
        post_filter_notes = '\n'.join([f'    - {filter}: {n}' for filter, n in post_filters.items() if not filter == ''])

        logging.info(f"  - {n_total} isomers tried, {n_isomers} isomers (from {n_complexes} complexes) successfully assembled.")
        if post_filter_notes != '':
            logging.info(f"  - {n_total - n_isomers} isomers failed because of filters:")
            logging.info(post_filter_notes)

        return

    @staticmethod
    def _log_batch_title_and_settings(batch_settings: Dict[Any, Any]) -> None:
        """
        Log the title and settings of the batch.
        :param batch_settings: Dictionary with the batch settings.
        :return: None
        """
        batch_title = f'  {batch_settings["name"]}  '
        logging.info(f'{batch_title:=^80}')
        logging.info(f"User-defined settings for {batch_settings['name']}:")
        for key, value in batch_settings.items():
            logging.info(f"    {key: <30}{value}")

        return
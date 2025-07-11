# Standard library imports
from typing import Union, Dict, Any, List, Tuple, Optional
from tqdm import tqdm
import datetime
import logging
import random
import sys

# DART specific imports
from DARTassembler.src.assembler.isomer import Isomer, IsomerFactory
from DARTassembler.src.assembler.output import AssemblerOutput, BatchAssemblerOutput, ComplexAssemblerOutput
from DARTassembler.src.metalig.utils_molecule import get_standardized_stoichiometry_from_atoms_list
from DARTassembler.src.assembler.utils import generate_pronounceable_word
from DARTassembler.src.constants.paths import default_assembler_yml_path
from DARTassembler.src.assembler.ligands import LigandChoice
from DARTassembler.src.modules.modules import BaseModule
from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.misc.io import read_yaml

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
                 concatenate_xyz: bool = True,
                 verbosity: int = 2,
                 same_isomer_names: bool = True,
                 complex_name_length: int = 8,
                 n_max_ligands: Union[int, None] = None,
                 ):
        """
        The main class for the DART workflow. It assembles all isomers of complexes from specified ligand databases and metal centers. Finally, all assembled complexes are saved to the output directory as .xyz and .json files, together with a csv file containing information about the assembly run.
        :param output_directory: Directory where the output files are saved. If it does not exist, it is created.
        :param concatenate_xyz: If True, a single concatenated .xyz file with all assembled complexes is created, for easier browsing in ase.
        :param verbosity: If 0, only errors are printed. If 1, warnings are printed. If 2, info messages are printed. If 3, debug messages are printed.
        :param same_isomer_names: If True, subsequent isomers of the same complex are given the same name as the first isomer, but with a number appended to the end. If False, each isomer gets a completely unique name.
        :param complex_name_length: The length of the randomly generated complex name. Automatically increases by 1 whenever a randomly generated name is duplicated to avoid duplicate names.
        :param n_max_ligands: Maximum number of ligands to load from the ligand databases provided in each batch. If None, all ligands are loaded. Useful to speed up test runs.
        """
        super().__init__()
        self.output_directory = Path(output_directory).resolve()
        self.concatenate_xyz = concatenate_xyz
        self.verbosity = verbosity
        self.same_isomer_names = same_isomer_names
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
        self.all_assembled_isomer_names = []
        self.successfully_assembled_isomer_names = []
        self.assembled_complex_json_paths = []

        self.batches_args = {**locals()}
        self.batches_args.pop('self')
        self.input_args = {**self.init_args, **self.batches_args}
        self._log_global_info()

        self._check_batch_settings(batches)   # todo outcomment when ready

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
            duplicate_names = df_batch['complex_name'][df_batch['complex_name'].duplicated()].values
            assert len(duplicate_names) == 0, f"Duplicate complex names in batch {batch}: {duplicate_names}. Please report this issue to our GitHub page."

        return

    def _log_summary(self) -> None:
        """Log nice summary per batch."""
        batch_summary_title = '  Summary per batch  '
        logging.info(f'{batch_summary_title:=^80}')
        for batch_idx, batch in enumerate(self.batches):
            df = self.df_info[self.df_info['batch_idx'] == batch_idx]
            if df.empty:
                # Exclude batches that do not gen any TMCs else following code will fail
                logging.info(f"Batch[{batch_idx}]: {batch['name']} --> No complexes assembled.")
                continue

            batch_name = df['batch_name'].iloc[0]
            assert batch_name == batch['name']  # todo Use instead of df and remove assert.
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
        logging.info(f'Output directory: {self.output_directory}')
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
        n_complexes = df[df['success']]['graph_hash'].nunique()

        # Output statistics how many isomers failed each filter
        post_filters = df['note'].value_counts()
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

    def _run_batch(self,
                    name: str,
                    batch_idx: int,
                    target_vectors: list[list[list[float]]],
                    metal_centers: Union[str, tuple[str, tuple[float, float, float]]],
                    n_max_complexes: Union[int, str],
                    swap_groups: list[int] = None,
                    filter_duplicate_isomers: bool = True,
                    filter_clashing_structures: bool = True,
                    filter_clashing_structures_cov_radii_buffer: float = 0.0,
                    check_metal_clashes: bool = False,
                    filter_duplicate_isomers_method: str = "fingerprint",
                    filter_duplicate_isomers_grid_size: int = 9,
                    isomer_comparison_mode: str = "max_diff",
                    isomer_comparison_grouping_mode: str = "cluster",
                    isomer_comparison_grouping_cutoff: float = 0.5,
                    opt_mono_rot: bool = True,
                    total_ligand_charges: int = None,
                    ligand_db_files: Union[list[str], str] = 'metalig',
                    ligand_origins: list[tuple[float, float, float]] = None,
                    complex_name_appendix: str = '',
                    random_seed: int = 0,
                    ) -> None:
        """
        Run the assembly for one batch.
        :param name: Name of the batch.
        :param batch_idx: Index of the batch in the list of batches.
        :param target_vectors: List of target vectors for the complexes to be assembled. Each target vector is a list of lists of floats/ints.
        :param metal_centers: Metal centers to be used for the complexes. Can be a string (e.g. 'Ru') or a tuple of (metal symbol, (x, y, z)).
        :param n_max_complexes: Maximum number of complexes to be assembled in this batch. If 'all', all complexes are assembled.
        :param total_ligand_charges: Choose ligands so that the total charge of all ligands in the complex is equal to this value. If None, any ligand combination is assembled.
        :param ligand_db_files: List of ligand database files to be used for the batch. If None, defaults to ['metalig'] for each set of target vectors.
        :param ligand_origins: Coordinates for each ligand, which will be used as the origin of the ligand rotation in the complex. If None, defaults to the center of all metal center coordinates for each ligand.
        :param complex_name_appendix: Appendix to be added to each complex name. Defaults to ''.
        :param random_seed: Random seed for reproducibility. If None, defaults to the batch_idx so that each batch is deterministic but different from each other.
        :return: None
        """
        # Set random seed for reproducibility. Do this batch-wise so every batch is reproducible independently.
        if random_seed is None:
            random_seed = batch_idx
        random.seed(random_seed)

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
        self.complex_name_appendix = complex_name_appendix
        self.swap_groups = swap_groups
        self.filter_duplicate_isomers = filter_duplicate_isomers
        self.filter_clashing_structures = filter_clashing_structures
        self.filter_clashing_structures_cov_radii_buffer = filter_clashing_structures_cov_radii_buffer
        self.check_metal_clashes = check_metal_clashes
        self.filter_duplicate_isomers_method = filter_duplicate_isomers_method
        self.filter_duplicate_isomers_grid_size = filter_duplicate_isomers_grid_size
        self.isomer_comparison_mode = isomer_comparison_mode
        self.isomer_comparison_grouping_mode = isomer_comparison_grouping_mode
        self.isomer_comparison_grouping_cutoff = isomer_comparison_grouping_cutoff
        self.opt_mono_rot = opt_mono_rot
        self.batch_output_path = Path(self.gbl_outcontrol.batch_dir, self.batch_name)
        self.batch_outcontrol = BatchAssemblerOutput(self.batch_output_path)

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
        total = self.n_max_complexes if self.n_max_complexes == 'all' else None
        progressbar = tqdm(desc='Assembling complexes', unit=' complexes', file=sys.stdout, total=total)

        batch_sum_assembled_complexes = 0  # Number of assembled complexes in this batch
        while ligand_choice.if_make_more_complexes(batch_sum_assembled_complexes):

            # Choose ligands for complex
            try:
                ligands = next(ligand_combinations)
            except StopIteration:
                break # If all ligand combinations are exhausted, stop the batch

            factory = IsomerFactory(
                                    ligands=ligands,
                                    target_vectors=self.target_vectors,
                                    ligand_origins=self.ligand_origins,
                                    metal_centers=self.metal_centers,
                                    filter_duplicate_isomers=self.filter_duplicate_isomers,
                                    filter_clashing_structures=self.filter_clashing_structures,
                                    filter_clashing_structures_cov_radii_buffer=self.filter_clashing_structures_cov_radii_buffer,
                                    check_metal_clashes=self.check_metal_clashes,
                                    filter_duplicate_isomers_method=self.filter_duplicate_isomers_method,
                                    filter_duplicate_isomers_grid_size=self.filter_duplicate_isomers_grid_size,
                                    isomer_comparison_mode=self.isomer_comparison_mode,
                                    isomer_comparison_grouping_mode=self.isomer_comparison_grouping_mode,
                                    isomer_comparison_grouping_cutoff=self.isomer_comparison_grouping_cutoff,
                                    swap_groups=self.swap_groups,
                                    opt_mono_rot= self.opt_mono_rot
                                    )
            isomers = factory.generate_isomers()

            # Sort isomers so that the ones without warnings are first, so that the indices of the names are consecutive.
            successful_isomers = [isomer for isomer in isomers if isomer.warning == '']
            unsuccessful_isomers = [isomer for isomer in isomers if isomer.warning != '']
            isomers = successful_isomers + unsuccessful_isomers  # Put successful isomers first, then unsuccessful ones, while keeping the order of isomers within each group.

            # Post-processing of isomers
            logging.debug("Post-processing isomers")
            isomer_idx = 1  # Index for naming the isomer
            for isomer in isomers:
                self._save_assembled_isomer(isomer=isomer, isomer_idx=isomer_idx)
                isomer_idx += 1

            # Update counters if at least one isomer was successfully assembled for this complex
            any_good_isomers = any(isomer.warning == '' for isomer in isomers)
            if any_good_isomers:
                batch_sum_assembled_complexes += 1
                progressbar.update(1)
            logging.debug("Leaving post-processing")

        progressbar.close()

        return

    def _get_ligand_databases(self) -> list[LigandDB]:
        ligand_databases = []
        for target_vectors, ligand_db_filepath in zip(self.target_vectors, self.ligand_db_files):
            if not ligand_db_filepath in self._loaded_ligand_databases:
                self._loaded_ligand_databases[ligand_db_filepath] = LigandDB.from_json(ligand_db_filepath, n_max=self.n_max_ligands)

            # Reduce to the ligands which have the correct n_eff_denticity for the specified target vectors
            n_target_vectors = len(target_vectors)
            database = {name: ligand for name, ligand in self._loaded_ligand_databases[ligand_db_filepath].db.items() if ligand.n_eff_denticities == n_target_vectors}
            database = LigandDB(database)  # Convert to LigandDB object

            if not database.db:
                raise ValueError(f"No ligands found in database {ligand_db_filepath} with the same `n_eff_denticities` as the number of specified target vectors ({n_target_vectors}). Please check your input ligand database `{ligand_db_filepath}`.")
            ligand_databases.append(database)

        return ligand_databases

    @staticmethod
    def _extract_ligand_origins_and_vectors(geometry: List[dict]) -> Tuple[List[Dict[str, List[float]]], List[List[float]]]:
        """
        Extracts and processes ligand origins and vectors from the geometry.

        :param geometry: List of fragments from batch["geometry"]
        :return: (target_vectors, ligand_origins)
        """
        axis_map = {
            "x": [1.0, 0.0, 0.0],
            "y": [0.0, 1.0, 0.0],
            "z": [0.0, 0.0, 1.0],
            "-x": [-1.0, 0.0, 0.0],
            "-y": [0.0, -1.0, 0.0],
            "-z": [0.0, 0.0, -1.0],
        }

        def _parse_vector(value, key: str, name: str) -> List[float]:
            if isinstance(value, str):
                vec = axis_map.get(value.strip().lower())
                if vec is None:
                    raise ValueError(
                        f"Fatal Error: Invalid string '{value}' in '{key}' for ligand '{name}'. "
                        f"Expected one of: {', '.join(axis_map.keys())}."
                    )
                return vec
            elif isinstance(value, list) and len(value) == 3:
                if not all(isinstance(v, (int, float)) for v in value):
                    raise ValueError(f"Fatal Error: '{key}' in ligand '{name}' must be a 3-element list of numbers.")
                return [float(v) for v in value]
            else:
                raise ValueError(f"Fatal Error: '{key}' in ligand '{name}' must be a symbolic axis string or a 3-element list.")

        target_vectors = []
        ligand_origins = []

        for fragment in geometry:
            assert len(fragment) == 1, f"{len(fragment)} keys present for {fragment}, expected 1"
            name = next(iter(fragment))
            if name.startswith("ligand"):
                spec = fragment[name]
                ligand_origins.append(spec["origin"])

                vector_dict = {
                    key: _parse_vector(value, key=key, name=name)
                    for key, value in spec.items()
                    if key.startswith("vector")
                }

                target_vectors.append(vector_dict)

        return target_vectors, ligand_origins

    def _save_assembled_isomer(self, isomer: Isomer, isomer_idx: int):
        """
        Save the successfully assembled complex to the output files.
        """
        success = True if isomer.warning == '' else False

        name = self._get_unique_complex_name(complex=isomer, isomer_idx=isomer_idx)

        total_ligand_charges = sum(isomer.ligand_info['charges'])  # Don't take the global property self.total_ligand_charges in case it is None
        isomer.global_props = {  # Overwrite potentially existing global properties with the new ones
            'complex_name': name,
            'stoichiometry': isomer.stoichiometry,
            'total_ligand_charges': total_ligand_charges,
            'graph_hash': isomer.graph_hash,
            'batch_name': self.batch_name,
        }

        # Add a comment to the xyz file with complex and ligand names for easier identification
        xyz_comment = f'complex_name: {name}, ligand_unique_names: ({", ".join(isomer.ligand_info["unique_names"])}), note: {isomer.warning}'
        xyz_string = isomer.get_xyz_string(comment=xyz_comment)

        # Save to complex directory
        if success:
            complex_dir = Path(self.batch_outcontrol.complex_dir, name)
            complex_outcontrol = ComplexAssemblerOutput(complex_dir)
            complex_outcontrol.save_all_complex_data(complex=isomer)
            complex_outcontrol.save_structure(xyz_string)

            # Keep track of number and names of assembled complexes
            self.successfully_assembled_isomer_names.append(name)
            self.assembled_complex_json_paths.append(complex_outcontrol.data_path)

        # Save to concatenated xyz file of this batch
        if self.concatenate_xyz:
            self.batch_outcontrol.save_xyz(xyz_string, success=success, append=True)

        # Save data for csv file.
        complex_idx = (len(self.successfully_assembled_isomer_names) - 1) if success else None
        self._add_batch_info(complex=isomer, success=success, complex_idx=complex_idx)

        self.all_assembled_isomer_names.append(name)

        return

    def _get_unique_complex_name(self, complex, isomer_idx, decimals=6) -> str:
        """
        Returns a unique new name for the complex. If the complex is a subsequent isomer, the name is based on the first isomer.
        """
        if isomer_idx > 1 and self.same_isomer_names:  # subsequent isomers
            # Fix subsequent isomer to always have the same name as the first isomer, but counting up.
            n_digits_last_isomer = len(str(isomer_idx - 1))
            n_digits_appendix = len(self.complex_name_appendix)
            n_digits_remove = n_digits_last_isomer + n_digits_appendix
            last_isomers_name = self.all_assembled_isomer_names[-1]
            last_isomers_stem = last_isomers_name[:-n_digits_remove]
            # Check that we can reconstruct the last isomers name.
            assert last_isomers_name == last_isomers_stem + str(isomer_idx - 1) + self.complex_name_appendix, f'The complex name seems to work different than implemented.'
            # Construct the new isomers name after the same rules as above.
            name = last_isomers_stem + str(isomer_idx) + self.complex_name_appendix
            assert not name in self.all_assembled_isomer_names, f"Complex name {name} already exists in the assembled complex names list even though it is a subsequent isomer."
        else:
            # Generate new name for new complex.
            complex_name_length = self.complex_name_length
            while True:  # emulate a do-while loop
                # Get a random name for the complex
                if self.same_isomer_names:
                    hash_string = complex.graph_hash
                else:
                    xyz = complex._get_xyz_array()
                    sorted_indices = np.lexsort((xyz[:, 2], xyz[:, 1], xyz[:, 0]), axis=0)
                    xyz = np.round(xyz, decimals=decimals)  # round to 6 decimals to get rid of numerical noise
                    xyz = xyz[sorted_indices]
                    elements = [el for _, el in sorted(zip(sorted_indices, complex.atomic_props['atoms']))]  # sort elements according to xyz
                    hash_string = str(elements) + str(xyz)  # make hash string

                # Generate a pronounceable word from the hash
                name = generate_pronounceable_word(length=complex_name_length, seed=hash_string)

                # If the name is based on the graph hash AND there are multiple isomers, add a number to the end of each name, otherwise start without a number.
                if self.same_isomer_names:  # Names based on graph hash
                    assert isomer_idx == 1, f'Isomer idx that is {isomer_idx} should be 1 here because subsequent isomers are handled differently.'
                    name = name + str(1)

                # Add the specified appendix to the name
                name += self.complex_name_appendix

                # If the name is already used, redo name generation with one more character. For the next complex, it starts with the original character length again.
                if name in self.all_assembled_isomer_names:
                    complex_name_length += 1
                    continue
                else:
                    break  # name is unique, break the loop

        return name

    def _add_batch_info(self, complex: Isomer, success, complex_idx: int) -> None:
        """
        Add information about the batch to the batch info variable which will be saved to the batch info file.
        """
        elements = complex.get_metal_symbols()
        metal_stoi = get_standardized_stoichiometry_from_atoms_list(elements)
        data = {
            'success': success,
            'complex_idx': complex_idx,
            'complex_name': complex.global_props['complex_name'],
            'stoichiometry': complex.global_props['stoichiometry'],
            'graph_hash': complex.global_props['graph_hash'],
            'note': complex.warning,
            'ligand_unique_names': complex.ligand_info['unique_names'],
            'ligand_geometries': complex.ligand_info['geometries'],
            'ligand_stoichiometries': complex.ligand_info['stoichiometries'],
            'ligand_charges': complex.ligand_info['charges'],
            'ligand_donors': complex.ligand_info['donors'],
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
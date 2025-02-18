import shutil
import sys
import datetime
from tqdm import tqdm
import random
import warnings
import pandas as pd
from pathlib import Path
from typing import Union
import ase
from copy import deepcopy
import numpy as np

from DARTassembler.src.constants.Paths import projectpath

warnings.simplefilter("always")
import logging
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') # Disable rdkit warnings

from DARTassembler.src.assembly.utilities_assembly import generate_pronounceable_word
from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from DARTassembler.src.assembly.ligands import LigandChoice
from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand
from DARTassembler.src.ligand_extraction.io_custom import load_json
from DARTassembler.src.assembly.Assembly_Input import LigandCombinationError
from DARTassembler.src.assembly.Assembly_Output import AssemblyOutput, BatchAssemblyOutput, ComplexAssemblyOutput
from dev.Assembler_revision_jan_2025.assembler.utilities import AssembledIsomer
from DARTassembler.src.ligand_extraction.utilities_Molecule import get_standardized_stoichiometry_from_atoms_list





class DARTAssembly(object):

    def __init__(
                    self,
                    output_directory: Union[str, Path] = 'DARTassembler',
                    concatenate_xyz: bool = True,
                    verbosity: int = 2,
                    same_isomer_names: bool = True,
                    complex_name_length: int = 8,
                    n_max_ligands: Union[int,None] = None,
                    pre_delete: bool = False
                 ):
        self.output_directory = Path(output_directory).resolve()
        self.concatenate_xyz = concatenate_xyz
        self.verbosity = verbosity
        self.same_isomer_names = same_isomer_names
        self.complex_name_length = complex_name_length
        self.n_max_ligands = n_max_ligands
        self.pre_delete = pre_delete

        # Initialize some necessary variables for later
        self.df_info = None
        self.random_seed = None

        # Keep track of the input arguments
        self.init_args = {**locals()}
        self.init_args.pop('self')

        # Delete the output directory if it already exists to avoid mixing of old and new data
        if self.pre_delete:
            shutil.rmtree(self.output_directory, ignore_errors=True)

        # Set up the output directories
        self.gbl_outcontrol = AssemblyOutput(outdir=self.output_directory)

        # Set up logging
        verbosity2logging = {0: logging.ERROR, 1: logging.WARNING, 2: logging.INFO, 3: logging.DEBUG}
        stream_handler = logging.StreamHandler(stream=sys.stdout)
        stream_handler.setLevel(verbosity2logging[self.verbosity])
        # Print to stdout
        logging.basicConfig(level=verbosity2logging[self.verbosity], format='%(message)s', handlers=[logging.FileHandler(self.gbl_outcontrol.log_path, mode='w'), stream_handler])


    def run_batches(self, batches: list[dict]) -> None:
        """
        Runs the whole assembly for all batches specified in the assembly input file.
        @param batches: list of dictionaries with the batch settings. See the documentation for more information.
        @:return None
        """
        # Save yml file with input arguments to output directory
        self.batches_args = {**locals()}
        self.batches_args.pop('self')
        self.input_args = {**self.init_args, **self.batches_args}
        self.gbl_outcontrol.save_settings(self.input_args)

        self.batches = batches
        self.n_batches = len(self.batches)
        self.df_info = []
        self.assembled_complex_names = []
        self.assembled_complex_json_paths = []
        self.last_ligand_db_path = None     # to avoid reloading the same ligand database in the next batch

        self._log_global_info()

        start = datetime.datetime.now()
        for idx, batch_settings in enumerate(self.batches):
            # Set batch settings for the batch run
            self.batch_name = batch_settings['name']
            self.ligand_db_files = batch_settings['ligand_db_files']
            self.max_num_complexes = batch_settings['max_num_complexes']
            self.random_seed = batch_settings['random_seed']
            self.total_charge = batch_settings['total_charge']
            self.total_metal_oxidation_state = batch_settings['total_metal_oxidation_state']
            self.target_vectors = batch_settings['target_vectors']
            self.ligand_origins = batch_settings['ligand_origins']
            self.metal_centers = batch_settings['metal_centers']
            self.complex_name_appendix = batch_settings['complex_name_appendix']
            self.geometry_modifier_filepath = batch_settings['geometry_modifier_filepath']

            self.batch_output_path = Path(self.gbl_outcontrol.batch_dir, self.batch_name)
            self.batch_outcontrol = BatchAssemblyOutput(self.batch_output_path)
            self.batch_idx = idx

            self._log_batch_title_and_settings(batch_settings)
            self._run_batch()  # run the batch assembly

        self.runtime = datetime.datetime.now() - start  # keep track of the runtime to display later
        self._make_and_save_output_csv()
        self._log_summary()
        self._final_checks()

        return

    def _make_and_save_output_csv(self) -> None:
        """Save output info csv of all attempts."""
        self.df_info = pd.DataFrame(self.df_info)
        self.df_info['attempt'] = self.df_info.index
        self.df_info = self.df_info[['attempt'] + [col for col in self.df_info.columns if col != 'attempt']]  # Move attempt column to front
        self.gbl_outcontrol.save_run_info_table(self.df_info)

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
            batch_name = df['batch_name'].iloc[0]
            logging.info(f"Batch {batch_idx} ('{batch_name}'):")
            self._log_success_rate(df)

        # Print total summary of run
        total_summary_title = '  Total summary of DART Assembler run  '
        logging.info(f'{total_summary_title:=^80}')
        self._log_success_rate(self.df_info)
        n_isomers = self.df_info['success'].sum()
        n_complexes = self.df_info['graph_hash'].nunique()
        logging.info(f"DART Assembler output files saved to directory `{self.output_directory.name}`.")

        # The runtime is printed but not logged, so that slight differences in the runtime do not cause the integration tests to fail.
        if self.verbosity > 1:
            print(f"Total runtime for assembling {n_isomers} isomers (from {n_complexes} complexes): {self.runtime}")

        logging.info('Done! All complexes assembled. Exiting DART Assembler Module.')

        return

    def _log_global_info(self) -> None:
        """Log global information about the run."""
        logging.info('Starting DART Assembler Module.')
        logging.info(f'Output directory: `{self.output_directory.name}`')
        plural = 'es' if self.n_batches > 1 else ''                # print plural or singular in next line
        logging.info(f"Running {self.n_batches} batch{plural}...")
        logging.info(f"User-defined global settings:")
        for key, value in self.init_args.items():
            logging.info(f"    {key: <30}{value}")

        return

    def _log_success_rate(self, df):
        """Log success rate of the run."""
        n_total = len(df)
        n_isomers = df['success'].sum()
        n_complexes = df['graph_hash'].nunique()

        # Output statistics how many isomers failed each filter
        post_filters = df['note'].value_counts()
        post_filter_notes = '\n'.join([f'    - {filter}: {n}' for filter, n in post_filters.items() if not filter == ''])

        logging.info(f"  - {n_total} isomers tried, {n_isomers} isomers (from {n_complexes} complexes) successfully assembled.")
        if post_filter_notes != '':
            logging.info(f"  - {n_total - n_isomers} isomers failed because of filters:")
            logging.info(post_filter_notes)

        return

    def _log_batch_title_and_settings(self, batch_settings: dict):
        batch_title = f'  Batch {self.batch_idx}: {self.batch_name}  '
        logging.info(f'{batch_title:=^80}')
        logging.info(f"User-defined settings for batch {self.batch_idx}:")
        for key, value in batch_settings.items():
            logging.info(f"    {key: <30}{value}")

    def _run_batch(self):
        """
        Run the assembly for one batch.
        """
        # Set random seed for reproducibility. Do this batch-wise so every batch is reproducible independently.
        random.seed(self.random_seed)

        # Here we load the ligand database and avoid reloading the same ligand database if it is the same as the last one
        if self._check_if_reload_database():
            self.ligand_db = self._get_ligand_db()

        # Set up an iterator for the ligand combinations
        ligand_choice = LigandChoice(
                                database=self.ligand_db,
                                metal_oxidation_state=self.total_metal_oxidation_state,
                                total_complex_charge=self.total_charge,
                                max_num_assembled_complexes=self.max_num_complexes,
                                )
        ligand_combinations = ligand_choice.choose_ligands()

        # Set progress bar with or without final number of assembled complexes
        total = self.max_num_complexes if self.max_num_complexes == 'all' else None
        progressbar = tqdm(desc='Assembling complexes', unit=' complexes', file=sys.stdout, total=total)

        batch_sum_assembled_complexes = 0  # Number of assembled complexes in this batch
        while ligand_choice.if_make_more_complexes(batch_sum_assembled_complexes):

            # Choose ligands for complex
            try:
                ligands = next(ligand_combinations)
            except StopIteration:
                break # If all ligand combinations are exhausted, stop the batch

            isomers, warnings = AssembledIsomer.from_ligands_and_metal_centers(
                                                                        ligands=ligands,
                                                                        target_vectors=self.target_vectors,
                                                                        ligand_origins=self.ligand_origins,
                                                                        metal_centers=self.metal_centers
                                                                        )

            # Post-processing of isomers
            logging.debug("Post-processing isomers...")
            isomer_idx = 1  # Index for naming the isomer
            for isomer, warning in zip(isomers, warnings):
                self._save_assembled_isomer(isomer=isomer, ligands=ligands, isomer_idx=isomer_idx, note=warning)
                isomer_idx += 1

            # Update counters if at least one isomer was successfully assembled for this complex
            any_good_isomers = any(warning == '' for warning in warnings)
            if any_good_isomers:
                batch_sum_assembled_complexes += 1
                progressbar.update(1)
            logging.debug("Leaving post-processing")

        progressbar.close()

        return

    def _modify_ligand_geometry(self, geometry_modifier_path: Union[str, Path], building_blocks: dict):
        old_mol, new_mol = ase.io.read(geometry_modifier_path, index=':', format='xyz')
        coordinates_for_matching = old_mol.get_positions()
        new_coordinates = new_mol.get_positions()
        match_atoms = old_mol.get_chemical_symbols()

        new_stk_ligand_building_blocks_list = {}
        for idx, ligand in building_blocks.items():
            coordinates_of_ligand = ligand.get_position_matrix()
            ligand_atoms = [DART_Element(atom.get_atomic_number()).symbol for atom in ligand.get_atoms()]

            new_ligand_coordinates = deepcopy(coordinates_of_ligand)
            for lig_coord_idx, (lig_coord, lig_atom) in enumerate(zip(coordinates_of_ligand, ligand_atoms)):
                if lig_atom == 'Hg':
                    continue

                for match_coord_idx, (match_coord, match_atom) in enumerate(zip(coordinates_for_matching, match_atoms)):
                    if lig_atom == match_atom and np.allclose(lig_coord, match_coord, atol=1e-5):
                        # Match of coordinates and elements found. Replace coordinates with new ones.
                        new_ligand_coordinates[lig_coord_idx] = new_coordinates[match_coord_idx]
                        break

            new_ligand = ligand.with_position_matrix(new_ligand_coordinates)
            new_stk_ligand_building_blocks_list[idx] = new_ligand

        return new_stk_ligand_building_blocks_list

    def _save_assembled_isomer(self, isomer: AssembledIsomer, ligands: list[RCA_Ligand], isomer_idx: int, note: str):
        """
        Save the successfully assembled complex to the output files.
        """
        success = True if note == '' else False

        name = self._get_unique_complex_name(complex=isomer, isomer_idx=isomer_idx)
        isomer.global_props['complex_name'] = name

        # Enrich the isomer with some information
        isomer.global_props['stoichiometry'] = isomer.stoichiometry
        isomer.global_props['charge'] = self.total_charge
        isomer.global_props['total_oxi_state'] = self.total_metal_oxidation_state
        isomer.global_props['graph_hash'] = isomer.graph_hash
        isomer.global_props['batch_name'] = self.batch_name

        # Save to complex directory
        if success:
            complex_dir = Path(self.batch_outcontrol.complex_dir, name)
            complex_outcontrol = ComplexAssemblyOutput(complex_dir)
            complex_outcontrol.save_all_complex_data(complex=isomer)
            complex_outcontrol.save_ligand_info(ligands)  # Save the ligand info
            complex_outcontrol.save_structure(isomer.get_xyz_file_format_string())

            # Keep track of number and names of assembled complexes
            self.assembled_complex_names.append(name)
            self.assembled_complex_json_paths.append(complex_outcontrol.data_path)

        # Save to concatenated xyz file of this batch
        self.batch_outcontrol.save_xyz(isomer.get_xyz_file_format_string(), success=success, append=True)

        # Save data for csv file.
        complex_idx = ( len(self.assembled_complex_names) - 1 ) if success else None
        self._add_batch_info(complex=isomer, success=success, note=note, complex_idx=complex_idx)

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
            last_isomers_name = self.assembled_complex_names[-1]
            last_isomers_stem = last_isomers_name[:-n_digits_remove]
            # Check that we can reconstruct the last isomers name.
            assert last_isomers_name == last_isomers_stem + str(isomer_idx - 1) + self.complex_name_appendix, f'The complex name seems to work different than implemented.'
            # Construct the new isomers name after the same rules as above.
            name = last_isomers_stem + str(isomer_idx) + self.complex_name_appendix
            assert not name in self.assembled_complex_names, f"Complex name {name} already exists in the assembled complex names list even though it is a subsequent isomer. This should be impossible."
        else:
            # Generate new name for new complex.
            complex_name_length = self.complex_name_length
            while True:     # emulate a do-while loop
                # Get a random name for the complex
                if self.same_isomer_names:
                    hash_string = complex.graph_hash
                else:
                    xyz = complex.mol.get_xyz_as_array()
                    sorted_indices = np.lexsort((xyz[:, 2], xyz[:, 1], xyz[:, 0]), axis=0)
                    xyz = np.round(xyz, decimals=decimals)  # round to 6 decimals to get rid of numerical noise
                    xyz = xyz[sorted_indices]
                    elements = [el for _, el in sorted(zip(sorted_indices, complex.mol.get_elements_list()))] # sort elements according to xyz
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
                if name in self.assembled_complex_names:
                    complex_name_length += 1
                    continue
                else:
                    break   # name is unique, break the loop

        return name

    def _add_batch_info(self, complex: AssembledIsomer, success, note: str, complex_idx: int) -> None:
        """
        Add information about the batch to the batch info variable which will be saved to the batch info file.
        """
        elements = complex.get_metal_center_atoms().get_chemical_symbols()
        metal_stoi = get_standardized_stoichiometry_from_atoms_list(elements)
        data = {
            'success': success,
            'complex_idx': complex_idx,
            'complex_name': complex.global_props['complex_name'],
            'stoichiometry': complex.global_props['stoichiometry'],
            'graph_hash': complex.global_props['graph_hash'],
            'note': note,
            'ligand_unique_names': complex.ligand_info['unique_names'],
            'ligand_geometries': complex.ligand_info['geometries'],
            'ligand_stoichiometries': complex.ligand_info['stoichiometries'],
            'ligand_charges': complex.ligand_info['charges'],
            'ligand_donors': complex.ligand_info['donors'],
            'batch_idx': self.batch_idx,
            'batch_name': self.batch_name,
            'metal_centers': metal_stoi,
            'total_oxi_state': self.total_metal_oxidation_state,
            'charge': self.total_charge,
            'random_seed': self.random_seed,
        }
        self.df_info.append(data)

        return

    def _check_if_reload_database(self) -> bool:
        """
        Check if the ligand database needs to be reloaded because any of the ligand json files have changed. Only for performance.
        """
        if isinstance(self.ligand_db_files, Path) and isinstance(self.last_ligand_db_path, Path):
            reload_database = self.last_ligand_db_path.resolve() != self.ligand_db_files.resolve()
        elif isinstance(self.ligand_db_files, list) and isinstance(self.last_ligand_db_path, list):
            reload_database = len(self.ligand_db_files) != len(self.last_ligand_db_path) or any([last_path.resolve() != path.resolve() for path, last_path in zip(self.ligand_db_files, self.last_ligand_db_path)])
        else:
            reload_database = True

        return reload_database

    def _get_ligand_db(self) -> Union[LigandDB, list[LigandDB]]:
        """
        Load the ligand database from the json files.
        @return: ligand database or list of ligand databases.
        """
        multiple_db = isinstance(self.ligand_db_files, list)
        if multiple_db:
            ligand_db = []
            for path in self.ligand_db_files:
                ligand_db.append(LigandDB.load_from_json(path))
                if len(ligand_db[-1]) == 0:
                    raise LigandCombinationError(
                        f"No ligands found in the ligand database {path}. Please check your ligand database files.")
        else:
            ligand_db = LigandDB.load_from_json(self.ligand_db_files)
            if len(ligand_db) == 0:
                raise LigandCombinationError(
                    f"No ligands found in the ligand database {self.ligand_db_files}. Please check your ligand database files.")
        self.last_ligand_db_path = self.ligand_db_files

        return ligand_db


if __name__ == '__main__':

    monopath = projectpath('dev/DART_refactoring_to_v1_1_0/data/assembler/example_ligand_db/monodentate.jsonlines')
    bidentpath = projectpath('dev/DART_refactoring_to_v1_1_0/data/assembler/example_ligand_db/bidentate.jsonlines')
    haptic_path = projectpath('dev/DART_refactoring_to_v1_1_0/data/assembler/example_ligand_db/haptic.jsonlines')

    # The input options as a dictionary for the refactored version of the assembler.
    options = {
        "output_directory": "data/assembler/data_output",
        "concatenate_xyz": True,
        "verbosity": 2,
        "same_isomer_names": True,
        "complex_name_length": 8,
        "n_max_ligands": 100,
        "batches": [{
                    "name": "First_batch",
                    "total_metal_oxidation_state": 2,
                    "total_charge": 0,
                    "ligand_db_files": [bidentpath, monopath, haptic_path],
                    "target_vectors": [
                                        [[1, 0, 0], [0, 1, 0]],
                                        [-1, 0, 0],
                                        [0, -1, 0],
                                        ],
                    "ligand_origins": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                    "metal_centers": [
                                        [['Ru', [0, 0, 0]]],
                                        [['Ru', [0, 0, 0]]],
                                        [['Ru', [0, 0, 0]]]
                                        ],
                    "max_num_complexes": 100,
                    "random_seed": 0,
                    "geometry_modifier_filepath": None,
                    "complex_name_appendix": 'test',
            }],
    }

    # out_json = load_json('/Users/timosommer/PhD/projects/DARTassembler/dev/DART_refactoring_to_v1_1_0/data/ABAFIFUQ_Ru_OH_data.json')

    batches = options.pop('batches')
    dart = DARTAssembly(**options, pre_delete=True)
    dart.run_batches(batches=batches)

    # Check if the code can read in a complex json file
    json_path = dart.assembled_complex_json_paths[0]
    isomer = AssembledIsomer.from_json(json_path)

    # out_json2 = load_json('/Users/timosommer/PhD/projects/DARTassembler/dev/DART_refactoring_to_v1_1_0/data/assembler/data_output/batches/First_batch/complexes/ADIYOWOZ1test/ADIYOWOZ1test_data.json')

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    old_dir = dart.output_directory.parent / Path('benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=dart.output_directory, old_dir=old_dir)
        test.compare_all()
        print('Test for installation passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')


    print('Done')
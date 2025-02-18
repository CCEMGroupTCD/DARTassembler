import shutil
import sys
from unittest.mock import patch
import datetime

from DARTassembler.src.assembly.utilities_assembly import generate_pronounceable_word, format_topologies
from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.assembly.Monkeypatch_stk import MONKEYPATCH_STK_SmartsFunctionalGroupFactory
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from DARTassembler.src.assembly.Assemble import PlacementRotation
from DARTassembler.src.assembly.ligands import LigandChoice
from DARTassembler.src.assembly.Isomer import BuildIsomers
from DARTassembler.src.assembly.Optimise import OPTIMISE
from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand
from DARTassembler.src.ligand_extraction.io_custom import load_json
from DARTassembler.src.assembly.Post_Filter import PostFilter
from tqdm import tqdm
import random
import warnings
import pandas as pd
from DARTassembler.src.assembly.TransitionMetalComplex import TransitionMetalComplex as TMC
from pathlib import Path
from typing import Union
from DARTassembler.src.assembly.Assembly_Input import AssemblyInput, LigandCombinationError, _isomers
from DARTassembler.src.assembly.Assembly_Output import AssemblyOutput, BatchAssemblyOutput, _gbl_optimization_movie, \
    ComplexAssemblyOutput, append_global_concatenated_xyz
import ase
from copy import deepcopy
import numpy as np

from DARTassembler.src.ligand_extraction.io_custom import load_unique_ligand_db
from DARTassembler.src.ligand_extraction.utilities_Molecule import get_standardized_stoichiometry_from_atoms_list
from dev.Assembler_revision_jan_2025.assembler.utilities import AssembledIsomer

warnings.simplefilter("always")
import logging
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') # Disable rdkit warnings


class DARTAssembly(object):

    def __init__(self,
                 assembly_input_path: Union[str, Path] = 'assembly_input.yml',
                 pre_delete: bool = False
                 ):

        # Save input paths
        self.assembly_input_path = assembly_input_path

        # Read in global settings and check if input looks correct
        self.settings = AssemblyInput(path=self.assembly_input_path)
        self.verbose = self.settings.verbose
        self.optimization_movie = self.settings.optimization_movie
        self.concatenate_xyz = self.settings.concatenate_xyz
        self.output_path = self.settings.Output_Path
        self.complex_name_length = self.settings.complex_name_length
        self.batches = self.settings.Batches
        self.n_batches = len(self.batches)
        self.same_isomer_names = self.settings.same_isomer_names

        if pre_delete:
            shutil.rmtree(self.output_path, ignore_errors=True)

        self.df_info = None
        self.gbl_outcontrol = AssemblyOutput(outdir=self.output_path)

        # Initialize some necessary variables
        self.random_seed = None
        self.topology_similarity = None
        self.metal_list = None

        # Set up logging
        verbosity2logging = {0: logging.ERROR, 1: logging.WARNING, 2: logging.INFO, 3: logging.DEBUG}
        stream_handler = logging.StreamHandler(stream=sys.stdout)
        stream_handler.setLevel(verbosity2logging[self.verbose])
        # Print to stdout
        logging.basicConfig(level=verbosity2logging[self.verbose], format='%(message)s', handlers=[logging.FileHandler(self.gbl_outcontrol.log_path, mode='w'),
                              stream_handler])


    def run_all_batches(self):
        """
        Runs the whole assembly for all batches specified in the assembly input file.
        """
        start = datetime.datetime.now()
        self.df_info = []
        self.assembled_complex_names = []
        self.last_ligand_db_path = None     # to avoid reloading the same ligand database in the next batch

        logging.info('Starting DART Assembler Module.')
        logging.info(f'Output directory: `{self.output_path.name}`')
        plural = 'es' if self.n_batches > 1 else ''                # print plural or singular in next line
        logging.info(f"Running {self.n_batches} batch{plural}...")
        for idx, batch_settings in enumerate(self.batches):
            # Set batch settings for the batch run
            self.batch_name, self.ligand_json, self.max_num_assembled_complexes, self.generate_isomer_instruction,\
            self.optimisation_instruction, self.random_seed, self.total_charge, metal_list, self.topology_similarity,\
            self.complex_name_appendix, self.geometry_modifier_filepath, bidentate_rotator, \
                                 = self.settings.check_and_return_batch_settings(batch_settings)

            self.batch_output_path = Path(self.gbl_outcontrol.batch_dir, self.batch_name)
            self.batch_idx = idx
            self.batch_outcontrol = BatchAssemblyOutput(self.batch_output_path)
            self.metal_type = metal_list[0]
            self.metal_ox_state = metal_list[1]
            self.build_options = {'bidentate_rotator': bidentate_rotator}
            self.multiple_db = isinstance(self.ligand_json, list)
            if self.generate_isomer_instruction == 'Generate All':
                self.multiple_isomers = True
            elif self.generate_isomer_instruction == 'Generate Lowest Energy':
                self.multiple_isomers = False
            else:
                raise ValueError(f"{_isomers} must be either 'Generate All' or 'Generate Lowest Energy', but is {self.generate_isomer_instruction}.")

            self.print_batch_title_and_settings(batch_settings)
            self.run_batch()  # run the batch assembly

        self.runtime = datetime.datetime.now() - start

        # Save output info csv of all attempts
        self.df_info = pd.DataFrame(self.df_info)
        self.df_info['attempt'] = self.df_info.index
        self.df_info = self.df_info[['attempt'] + [col for col in self.df_info.columns if col != 'attempt']]  # Move attempt column to front
        self.gbl_outcontrol.save_run_info_table(self.df_info)

        # Save yaml file with input settings
        self.gbl_outcontrol.save_settings(self.settings.global_settings)

        # Print nice summary per batch
        batch_summary_title = '  Summary per batch  '
        logging.info(f'{batch_summary_title:=^80}')
        for batch_idx, batch in enumerate(self.batches):
            df = self.df_info[self.df_info['batch_idx'] == batch_idx]
            batch_name = df['batch_name'].iloc[0]
            logging.info(f"Batch {batch_idx} ('{batch_name}'):")
            self.print_success_rate(df)

        # Print total summary of run
        total_summary_title = '  Total summary of DART Assembler run  '
        logging.info(f'{total_summary_title:=^80}')
        self.print_success_rate(self.df_info)
        n_success = self.df_info['success'].sum()
        logging.info(f"DART Assembler output files saved to directory `{self.output_path.name}`.")
        print(f"Total runtime for assembling {n_success} complexes: {self.runtime}")    # print to leave runtime out of log file for integration tests
        logging.info('Done! All complexes assembled. Exiting DART Assembler Module.')

        # Some final tests
        df_test_success = self.df_info[self.df_info['success']]
        batches = df_test_success['batch_idx'].unique()
        for batch in batches:
            df_batch = df_test_success[df_test_success['batch_idx'] == batch]
            # Check for duplicate complex names in the batch
            duplicate_names = df_batch['complex_name'][df_batch['complex_name'].duplicated()].values
            assert len(duplicate_names) == 0, f"Duplicate complex names in batch {batch}: {duplicate_names}. Please report this issue to our GitHub page."

        return

    def print_success_rate(self, df):
        n_success = df['success'].sum()
        n_total = len(df)
        post_filters = df['note'].value_counts()
        successful_assembly_notes = ['no optimization', 'optimized']
        post_filter_notes = '\n'.join([f'    - {filter}: {n}' for filter, n in post_filters.items() if filter not in successful_assembly_notes])

        logging.info(f"  - {n_total} complexes tried, {n_success} complexes successfully assembled.")
        if post_filter_notes != '':
            logging.info(f"  - {n_total - n_success} complexes failed because of filters:")
            logging.info(post_filter_notes)

        return

    def print_batch_title_and_settings(self, batch_settings: dict):
        batch_title = f'  Batch {self.batch_idx}: {self.batch_name}  '
        logging.info(f'{batch_title:=^80}')
        logging.info(f"User-defined settings for batch {self.batch_idx}:")
        for key, value in batch_settings.items():
            logging.info(f"    {key: <30}{value}")

    def run_batch(self):

        random.seed(int(self.random_seed))  # Set random seed for reproducibility

        # Here we load the ligand database and avoid reloading the same ligand database if it is the same as the last one
        if self.check_if_reload_database():
            self.ligand_db = self.get_ligand_db()

        Topology, Similarity = format_topologies(self.topology_similarity)

        # Choose ligands for each complex
        choice = LigandChoice(
            database=self.ligand_db,
            topology=Topology,
            instruction=Similarity,
            metal_oxidation_state=int(self.metal_ox_state),
            total_complex_charge=self.total_charge,
            max_num_assembled_complexes=self.max_num_assembled_complexes,
        )
        ligand_combinations = choice.choose_ligands()

        # Set progress bar with or without final number of assembled complexes
        if self.max_num_assembled_complexes == 'all':
            # If we don't know the final number of assembled complexes, we don't set the total number of iterations for the progress bar
            progressbar = tqdm(desc='Assembling complexes', unit=' complexes', file=sys.stdout)
        else: # self.max_num_assembled_complexes is an integer
            # If we know the final number of assembled complexes, we set the total number of iterations for the progress bar
            progressbar = tqdm(total=self.max_num_assembled_complexes, desc='Assembling complexes', unit=' complexes', file=sys.stdout)

        batch_sum_assembled_complexes = 0  # Number of assembled complexes produced.
        while choice.if_make_more_complexes(batch_sum_assembled_complexes):

            # 1. Choose Ligands for Complex
            try:
                ligands = next(ligand_combinations)
            except StopIteration:
                break # If all ligand combinations are exhausted, stop the batch
            ligands = list(ligands.values())

            # 2. Detect certain conditions which hinder complex assembly, e.g. tridentate non-planar
            complex_can_be_assembled = self.check_if_complex_can_be_assembled(ligands)
            if not complex_can_be_assembled:
                continue

            ligand_origins = [[0,0,0] for _ in ligands]
            ligand_target_vectors = [
                [[1, 0, 0], [0, 1, 0]],
                [-1, 0, 0],
                [0, -1, 0],
            ]
            isomers, warnings = AssembledIsomer.from_ligands_and_metal_centers(
                                                                        ligands=ligands,
                                                                        target_vectors=ligand_target_vectors,
                                                                        ligand_origins=ligand_origins,
                                                                        metal_centers=self.metal_type
                                                                        )

            # Post-processing of isomers
            logging.debug("Post-processing isomers...")
            isomer_idx = 1  # Index for naming the isomer
            for isomer, warning in zip(isomers, warnings):
                success = True if warning == '' else False
                self.save_assembled_isomer(isomer=isomer, success=success, ligands=ligands, isomer_idx=isomer_idx, note=warning)
                isomer_idx += 1

                if success:  # New complex successfully built.
                    batch_sum_assembled_complexes += 1
                    progressbar.update(1)
            logging.debug("Leaving post process")

        progressbar.close()

        return

    def modify_ligand_geometry(self, geometry_modifier_path: Union[str, Path], building_blocks: dict):
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

    def relax_and_check_structure(self, complex, building_blocks, ligands):
        """
        This function takes a complex, optionally relaxes it via forcefield and then checks if it is a good complex.
        """
        complex_is_good = PostFilter(isomer=complex,
                                     metal_centre=self.metal_type,
                                     metal_offset=-0.9,  # The more negative the number the more crap gets through
                                     ligand_atom_offset=0.20,  # default is 0.5
                                     building_blocks=building_blocks,
                                     ).closest_distance()

        if not complex_is_good:
            logging.debug("!!!Warning!!! -> None detect in optimiser -> Returning None")
            ff_movie = None
            note = 'clashing ligands'
        elif self.optimisation_instruction == False:
            ff_movie = None
            note = 'no optimization'
        else:  # Optimise via forcefield and re-check structure
            complex, building_blocks, ff_movie = OPTIMISE(
                isomer=complex,
                ligand_list=ligands,
                building_blocks=building_blocks,
                nsteps=50,
            ).Optimise_STK_Constructed_Molecule(return_ff_movie=True)

            complex_is_good = PostFilter(isomer=complex,
                                         metal_centre=self.metal_type,
                                         building_blocks=building_blocks,
                                         metal_offset=-0.9,
                                         ligand_atom_offset=0.5,
                                         ).post_optimisation_filter()
            note = 'broken bonds' if not complex_is_good else 'optimized'

        return complex, complex_is_good, ff_movie, note

    def save_assembled_isomer(self, isomer: AssembledIsomer, success: bool, ligands: list, isomer_idx: int, note: str):
        """
        Save the successfully assembled complex to the output files.
        """
        name = self.get_complex_name(complex=isomer, isomer_idx=isomer_idx)
        isomer.global_props['name'] = name

        # Save to complex directory
        if success:
            complex_dir = Path(self.batch_outcontrol.complex_dir, name)
            complex_outcontrol = ComplexAssemblyOutput(complex_dir)
            complex_outcontrol.save_all_complex_data(complex=isomer)

        # Save to concatenated xyz file of this batch
        self.batch_outcontrol.save_xyz(isomer.get_xyz_file_format_string(), success=success, append=True)

        # Save data for csv file
        complex_idx = len(self.assembled_complex_names)
        self.add_batch_info(success=success, reason=note, ligands=ligands, complex_idx=complex_idx, complex_name=name, complex_graph_hash=isomer.graph_hash)

        # Keep track of number and names of assembled complexes
        self.assembled_complex_names.append(name)

        return

    def save_failed_assembled_complex(self, complex: TMC, ligands: dict, note: str):
        """
        Save the successfully assembled complex to the output files.
        """
        self.add_batch_info(success=False, reason=note, ligands=ligands)

        # Save to concatenated xyz file of this batch
        xyz_string = complex.mol.get_xyz_file_format_string()
        self.batch_outcontrol.save_failed_xyz(xyz_string, append=True)

        return

    def get_complex_name(self, complex, isomer_idx, decimals=6) -> str:
        """
        Returns the name of the complex.
        """
        # Fix subsequent isomer to always have the same name as the first isomer, but counting up.
        if isomer_idx > 1 and self.same_isomer_names and self.multiple_isomers:  # subsequent isomers
            n_digits_last_isomer = len(str(isomer_idx - 1))
            n_digits_appendix = len(self.complex_name_appendix)
            n_digits_remove = n_digits_last_isomer + n_digits_appendix
            last_isomers_name = self.assembled_complex_names[-1]
            last_isomers_stem = last_isomers_name[:-n_digits_remove]
            # Check that we can reconstruct the last isomers name.
            assert last_isomers_name == last_isomers_stem + str(isomer_idx - 1) + self.complex_name_appendix, f'The complex name seems to work different than implemented.'
            # Now construct the new isomers name after the same rules as above.
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
                if self.same_isomer_names and self.multiple_isomers:  # Names based on graph hash
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

    def add_batch_info(self, success, ligands: list[RCA_Ligand], reason: str = '', complex_idx=None, complex_name: str=None, complex_graph_hash: str=None):
        """
        Add information about the batch to the batch info variable which will be saved to the batch info file.
        """
        ligand_names = tuple(ligand.unique_name for ligand in ligands)
        ligand_stoichiometries = tuple(ligand.stoichiometry for ligand in ligands)
        ligand_charges = tuple(ligand.pred_charge for ligand in ligands)
        ligand_donors = tuple('-'.join(sorted(ligand.local_elements)) for ligand in ligands)
        atoms = [self.metal_type] + [atom for ligand in ligands for atom in ligand.atomic_props['atoms']]
        stoichiometry = get_standardized_stoichiometry_from_atoms_list(atoms)

        data = {
            "success": success,
            "complex_idx": complex_idx,
            'complex_name': complex_name,
            "stoichiometry": stoichiometry,
            'graph_hash': complex_graph_hash,
            "note": reason,
            "ligand_names": ligand_names,
            "ligand_stoichiometries": ligand_stoichiometries,
            "ligand_charges": ligand_charges,
            "ligand_donors": ligand_donors,
            "batch_idx": self.batch_idx,
            "batch_name": self.batch_name,
            "metal": self.metal_type,
            "oxi_state": self.metal_ox_state,
            "total_charge": self.total_charge,
            "optimization": self.optimisation_instruction,
            "isomers": self.generate_isomer_instruction,
            "random_seed": self.random_seed,
        }
        self.df_info.append(data)

        return

    def check_if_complex_can_be_assembled(self, ligands: list[RCA_Ligand]):
        """
        Check if the complex can be assembled based on certain conditions to avoid errors.
        """
        # Non-planar tridentate and tetradentate ligands cannot be assembled
        dent_names = {3: 'tridentate', 4: 'tetradentate'}
        for ligand in ligands:
            if ligand.denticity in [3, 4]:
                is_planar = ligand.if_donors_planar(with_metal=True)
                if not is_planar:
                    self.add_batch_info(success=False, reason=f'non-planar {dent_names[ligand.denticity]}', ligands=ligands)
                    return False

        # Hydride ligand cannot be assembled
        hydride_found = False
        for ligand in ligands:
            if ((ligand.atomic_props['atoms'][0] == "H") or (ligand.atomic_props['atoms'][0] == "Se")) and (
                    len(ligand.atomic_props['atoms']) == 1):
                hydride_found = True
            else:
                pass
        if hydride_found:
            self.add_batch_info(success=False, reason='hydride', ligands=ligands)
            return False

        # Haptic ligands cannot be assembled
        for ligand in ligands:
            if ligand.has_neighboring_coordinating_atoms:
                self.add_batch_info(success=False, reason='haptic ligand', ligands=ligands)
                return False

        return True

    def check_if_reload_database(self):
        """
        Check if the ligand database needs to be reloaded because any of the ligand json files have changed. Only for performance.
        """
        if isinstance(self.ligand_json, Path) and isinstance(self.last_ligand_db_path, Path):
            reload_database = self.last_ligand_db_path.resolve() != self.ligand_json.resolve()
        elif isinstance(self.ligand_json, list) and isinstance(self.last_ligand_db_path, list):
            reload_database = len(self.ligand_json) != len(self.last_ligand_db_path) or any([last_path.resolve() != path.resolve() for path, last_path in zip(self.ligand_json, self.last_ligand_db_path)])
        else:
            reload_database = True

        return reload_database

    def get_ligand_db(self) -> Union[dict, list[dict]]:
        """
        Load the ligand database from the json files.
        @return: ligand database or list of ligand databases. The db are in the format {denticity: {charge: [ligand, ligand, ...]}}
        """
        if self.multiple_db:
            ligand_db = []
            for path in self.ligand_json:
                ligand_db.append(LigandDB.load_from_json(path).get_lig_db_in_old_format())
                if len(ligand_db[-1]) == 0:
                    raise LigandCombinationError(
                        f"No ligands found in the ligand database {path}. Please check your ligand database files.")
        else:
            ligand_db = LigandDB.load_from_json(self.ligand_json).get_lig_db_in_old_format()
            if len(ligand_db) == 0:
                raise LigandCombinationError(
                    f"No ligands found in the ligand database {self.ligand_json}. Please check your ligand database files.")
        self.last_ligand_db_path = self.ligand_json

        return ligand_db


if __name__ == '__main__':

    # The input options as a dictionary for the refactored version of the assembler.
    options = {
        "output_directory": "data/assembler/data_output",
        "ffmovie": False,
        "concatenate_xyz": True,
        "verbosity": 2,
        "same_isomer_names": True,
        "complex_name_length": 8,
        "batches": [{
                    "name": "First_batch",
                    "metal_center": "Fe",
                    "metal_oxidation_state": 2,
                    "total_charge": 0,
                    "geometry": "2-1-1",
                    "ligand_db_file": "test_metalig",
                    "max_num_complexes": 100,
                    "isomers": "all",
                    "random_seed": 0,
                    "forcefield": False,
                    "bidentate_rotator": "auto",
                    "geometry_modifier_filepath": None,
                    "complex_name_appendix": None,
            }],
    }

    # Set the path to the assembly input file
    assembly_input_path = 'assembler.yml'

    out_json = load_json('/Users/timosommer/PhD/projects/DARTassembler/dev/DART_refactoring_to_v1_1_0/data/ABAFIFUQ_Ru_OH_data.json')

    dart = DARTAssembly(assembly_input_path=assembly_input_path, pre_delete=True)
    dart.run_all_batches()

    out_json2 = load_json('/Users/timosommer/PhD/projects/DARTassembler/dev/DART_refactoring_to_v1_1_0/data/assembler/data_output/batches/First_batch/complexes/ACIJESIQ1/ACIJESIQ1_data.json')

    # #%% ==============    Doublecheck refactoring    ==================
    # from dev.test.Integration_Test import IntegrationTest
    # old_dir = dart.output_path.parent / Path('benchmark_data_output')
    # if old_dir.exists():
    #     test = IntegrationTest(new_dir=dart.output_path, old_dir=old_dir)
    #     test.compare_all()
    #     print('Test for installation passed!')
    # else:
    #     print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')


    print('Done')
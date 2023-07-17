import shutil
from unittest.mock import patch
from src05_Assembly_Refactor.Monkeypatch_stk import MONKEYPATCH_STK_SmartsFunctionalGroupFactory
from src01.DataBase import LigandDB
from src05_Assembly_Refactor.Assemble import PlacementRotation
from src05_Assembly_Refactor.ligands import ChooseRandomLigands, ChooseIterativeLigands
from src05_Assembly_Refactor.Isomer import BuildIsomers
from src05_Assembly_Refactor.Optimise import OPTIMISE
from src05_Assembly_Refactor.Post_Filter import PostFilter
import random
import warnings
from constants.Paths import project_path
from TransitionMetalComplex import TransitionMetalComplex as TMC
from pathlib import Path
from typing import Union
from src05_Assembly_Refactor.Assembly_Input import AssemblyInput
from src05_Assembly_Refactor.Assembly_Output import AssemblyOutput, BatchAssemblyOutput, _gbl_optimization_movie, _gbl_concatenated_xyz, ComplexAssemblyOutput

warnings.simplefilter("always")


class DARTAssembly(object):

    def __init__(self,
                 assembly_input_path: Union[str, Path] = 'assembly_input.yml'
                 ):
        """
        :param assembly_input_path: Path to the assembly input file.
        """

        # Save input paths
        self.assembly_input_path = assembly_input_path

        # Read in global settings and check if input looks correct
        self.settings = AssemblyInput(path=self.assembly_input_path, ensure_empty_output_dir=True)
        self.verbose = self.settings.verbose
        self.optimization_movie = self.settings.optimization_movie
        self.concatenate_xyz = self.settings.concatenate_xyz
        self.output_path = self.settings.Output_Path
        self.batches = self.settings.Batches

    def run_all_batches(self):
        """
        Runs the whole assembly for all batches specified in the assembly input file.
        """
        sum_assembled_complexes = 0
        for idx, batch_settings in enumerate(self.batches):
            batch_name, ligand_json, max_num_assembled_complexes, generate_isomer_instruction, optimisation_instruction, \
        random_seed, total_charge, metal_list, topology_similarity = self.settings.check_and_return_batch_settings(batch_settings)
            batch_output_path = Path(self.output_path, batch_name)
            batch = BatchAssembly(
                                    batch_name=batch_name,
                                    ligand_json=ligand_json,
                                    max_num_assembled_complexes=max_num_assembled_complexes,
                                    generate_isomer_instruction=generate_isomer_instruction,
                                    optimisation_instruction=optimisation_instruction,
                                    random_seed=random_seed,
                                    total_charge=total_charge,
                                    metal_type=metal_list[0],
                                    metal_ox_state=metal_list[1],
                                    metal_spin=metal_list[2],
                                    topology_similarity=topology_similarity,
                                    output_path=batch_output_path,
                                    assembly_input_path=self.assembly_input_path,
                                    batch_idx = idx,
                                    )
            sum_assembled_complexes += batch.run(complex_starting_index=sum_assembled_complexes)

        return

class BatchAssembly(object):

    def __init__(
                    self,
                    batch_name: str,
                    ligand_json: Union[str, Path],
                    max_num_assembled_complexes: int,
                    generate_isomer_instruction: str,
                    optimisation_instruction: str,
                    random_seed: int,
                    total_charge: int,
                    metal_type: str,
                    metal_ox_state: str,
                    metal_spin: int,
                    topology_similarity: float,
                    output_path: Union[str,Path],
                    assembly_input_path: Union[str, Path] = None,
                    batch_idx: int = None,
                    ):
        self.output_path = Path(output_path)

        self.batch_name = batch_name
        self.ligand_json = ligand_json
        self.max_num_assembled_complexes = max_num_assembled_complexes
        self.generate_isomer_instruction = generate_isomer_instruction
        self.optimisation_instruction = optimisation_instruction
        self.random_seed = random_seed
        self.total_charge = total_charge
        self.metal_type = metal_type
        self.metal_ox_state = metal_ox_state
        self.metal_spin = metal_spin
        self.topology_similarity = topology_similarity
        self.assembly_input_path = assembly_input_path
        self.batch_idx = batch_idx
        self.global_settings_controler = AssemblyInput(path=self.assembly_input_path)

        self.output_control = BatchAssemblyOutput(self.output_path)


    def run(self, complex_starting_index: int = 0):

        print(f"\n\n\n\n**********Batch: {self.batch_name}**********\n\n\n\n")

        # Here we load the ligand database
        F = LigandDB.from_json(json_=self.ligand_json,
                               type_="Ligand")  # We initiate the database in the RandomComplexAssembler
        RCA = PlacementRotation(database=F)

        ########################################################################
        # This section of code needs to be commented out if you do not want to
        # generate every combination of ligands but would rather build in a
        # random manner

        """IterativeLigands = ChooseIterativeLigands(database=RCA.ligand_dict,
                                                  top_list=topology_list,
                                                  charge=Total_Charge,
                                                  metal=metal_list,
                                                  random_seed=Random_Seed)
        IterativeLigands.obtain_all_combinations()"""
        ########################################################################

        j = 0  # Assembly iteration we are on
        Sum_Assembled_Complexes = 0  # Number of assembled complexes produced
        # Note: j is not always equal to Sum_Assembled_Complexes
        while Sum_Assembled_Complexes < self.max_num_assembled_complexes:
            print(
                f"###############################__Attempting_Assembly_of_Complex_#_{str(self.random_seed + j)}__###############################")
            #
            #
            # 1. Choose Random seed
            random.seed(int(self.random_seed) + j)

            #
            #
            # 3. Random Choice of Inputs
            Topology, Similarity = RCA.format_topologies(self.topology_similarity)

            #
            #
            # 3a. Choose Ligands Randomly
            # This section of code needs to be uncommented if you want complexes assembled randomly
            ligands = ChooseRandomLigands(database=RCA.ligand_dict,
                                          topology=Topology,
                                          instruction=Similarity,
                                          metal_oxidation_state=int(self.metal_ox_state),
                                          total_complex_charge=self.total_charge,
                                          max_attempts=1000000).choose_ligands()

            #
            #
            # 3b. Choose Ligands Iteratively
            # This section of code needs to be commented if you want complexes assembled randomly
            """ligands = IterativeLigands.choose_ligands(iteration=Random_Seed + j)"""

            #
            #
            # 4. Detect Tridentate Non-Planar
            if not RCA.planar_check_(ligands):
                warnings.warn(
                    "!!!Warning!!! -> The program has encountered a non_planar tridentate ligand, The program as of yet can not assemble with this ligand, skipping to next assembly")
                j += 1
                continue
            else:
                pass

            hydride_found = False
            for ligand in ligands.values():
                if ((ligand.atomic_props['atoms'][0] == "H") or (ligand.atomic_props['atoms'][0] == "Se")) and (
                        len(ligand.atomic_props['atoms']) == 1):
                    hydride_found = True
                else:
                    pass

            if hydride_found:
                warnings.warn(
                    "!!!Warning!!! -> The program has encountered a hydride. The program as of yet can not assemble with this ligand, skipping to next assembly.")
                j += 1
                continue

            #
            #
            # 5. Obtain rotated building blocks
            # Here we pass in our ligands and get out our stk building blocks
            with patch('stk.SmartsFunctionalGroupFactory', new=MONKEYPATCH_STK_SmartsFunctionalGroupFactory):   # This is a monkey patch to fix an inconvenience in stk where every molecule is sanitized with rdkit, which throws errors for some of our molecules
                stk_ligand_building_blocks_list, denticities = RCA.convert_ligand_to_building_block_for_complex(
                    ligands=ligands,
                    topology=Topology,
                    metal=self.metal_type)

            #
            #
            # 6. Generate Isomers
            Isomers = BuildIsomers(topology=self.topology_similarity,
                                   building_blocks_list=stk_ligand_building_blocks_list,
                                   metal_input=self.metal_type,
                                   charge_input=self.metal_ox_state,
                                   denticity_list=denticities,
                                   return_all_isomers=self.generate_isomer_instruction,
                                   opt_choice=self.optimisation_instruction,
                                   ligand_list=ligands)

            Assembled_Complex_list, Building_Block_list = Isomers.Generate()

            #
            #
            # 7. Post-Process
            # Post process includes error detection and optimization
            Post_Process_Complex_List = []
            print("entering post process")
            for complex, building_blocks in zip(Assembled_Complex_list, Building_Block_list):
                complex, complex_is_good, ff_movie = self.relax_and_check_structure(complex, building_blocks, ligands)

                tmc = TMC(compl=complex, ligands=ligands, metal=self.metal_type, metal_charge=int(self.metal_ox_state), spin=self.metal_spin)
                if complex_is_good:               # New complex successfully built.
                    self.save_successfully_assembled_complex(tmc, ff_movie, complex_idx=Sum_Assembled_Complexes+complex_starting_index)
                    Post_Process_Complex_List.append(complex)
                    Sum_Assembled_Complexes += 1
                else:
                    self.save_failed_assembled_complex(complex=tmc)
            print("Leaving post process")

            #
            #
            # 8. Format Outputs
            # Todo: this function should not be used anymore. It is only used for the old output format.
            RCA.output_controller_(list_of_complexes_wih_isomers=Post_Process_Complex_List,
                                   ligands=ligands,
                                   metal=self.metal_type,
                                   metal_ox_state=int(self.metal_ox_state),
                                   metal_multiplicity=self.metal_spin,
                                   view_complex=False,
                                   write_gaussian_input_files=False,
                                   output_directory=self.output_path, # todo: refactor to only use output_path
                                   )

            if j == 138746543956439563475683496736:
                exit()
            j += 1

        return Sum_Assembled_Complexes

    def relax_and_check_structure(self, complex, building_blocks, ligands):
        """
        This function takes a complex, optionally relaxes it via forcefield and then checks if it is a good complex.
        """
        complex_is_good = PostFilter(isomer=complex,
                                     metal_centre=self.metal_type,
                                     metal_offset=-0.9,  # The more negative the number the more crap gets through
                                     ligand_atom_offset=0.5,
                                     building_blocks=building_blocks,
                                     ).closest_distance()

        if not complex_is_good:
            warnings.warn("!!!Warning!!! -> None detect in optimiser -> Returning None")
            ff_movie = None
        elif self.optimisation_instruction == False:
            ff_movie = None
        else:   # Optimise via forcefield and re-check structure
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

        return complex, complex_is_good, ff_movie

    def save_successfully_assembled_complex(self, complex: TMC, ff_movie: str, complex_idx: int):
        """
        Save the successfully assembled complex to the output files.
        """
        # Save the force field movie
        if ff_movie is not None:
            self.output_control.save_passed_ff_movie(ff_movie)

            # Save to global optimization movie. Todo: Remove this once it is not needed anymore
            global_optimization_movie_path = Path(self.output_path.parent, _gbl_optimization_movie)
            with open(global_optimization_movie_path, "a") as f:
                f.write(ff_movie)

        # Save to concatenated xyz file of this batch
        xyz_string = complex.mol.get_xyz_file_format_string()
        self.output_control.save_passed_xyz(xyz_string, append=True)

        # Save to complex directory
        complex_name = self.get_complex_name(complex_idx)
        complex_dir = Path(self.output_control.complex_dir, complex_name)
        output_controler = ComplexAssemblyOutput(complex_dir)
        output_controler.save_all_complex_data(
                                                complex=complex,
                                                complex_idx=complex_idx,
                                                xyz_structure=xyz_string,
                                                ff_movie=ff_movie,
                                                assembly_input_path=self.assembly_input_path,
                                                batch_idx=self.batch_idx,
                                                )

        return

    def save_failed_assembled_complex(self, complex: TMC):
        """
        Save the successfully assembled complex to the output files.
        """
        # Save to concatenated xyz file of this batch
        xyz_string = complex.mol.get_xyz_file_format_string()
        self.output_control.save_failed_xyz(xyz_string, append=True)

        return

    def get_complex_name(self, complex_idx) -> str:
        """
        Returns the name of the complex based on the complex index.
        """
        return f"complex_{complex_idx}"



    # This is not working atm, it gives different results!!!
    # @classmethod
    # def from_assembly_input(cls, assembly_input_file: Union[str,Path], batch_idx: int, output_path: Union[str,Path]):
    #     global_settings = AssemblyInput(path=assembly_input_file)
    #     batch_settings = global_settings.Batches[batch_idx]
    #
    #     batch_name, ligand_json, max_num_assembled_complexes, generate_isomer_instruction, optimisation_instruction, \
    #     random_seed, total_charge, metal_list, topology_similarity = global_settings.check_and_return_batch_settings(batch_settings)
    #     metal_type = metal_list[0]
    #     metal_ox_state = str(metal_list[1])
    #     metal_spin = int(metal_list[2])
    #
    #     return cls(
    #                 batch_name=batch_name,
    #                 ligand_json=ligand_json,
    #                 max_num_assembled_complexes=max_num_assembled_complexes,
    #                 generate_isomer_instruction=generate_isomer_instruction,
    #                 optimisation_instruction=optimisation_instruction,
    #                 random_seed=random_seed,
    #                 total_charge=total_charge,
    #                 metal_type=metal_type,
    #                 metal_ox_state=metal_ox_state,
    #                 metal_spin=metal_spin,
    #                 topology_similarity=topology_similarity,
    #                 output_path=output_path
    #                 )



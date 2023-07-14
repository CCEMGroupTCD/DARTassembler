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
from pathlib import Path
from typing import Union
from src05_Assembly_Refactor.Assembly_Input import AssemblyInput

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
        self.settings = AssemblyInput(path=self.assembly_input_path)
        self.verbose = self.settings.verbose
        self.optimization_movie = self.settings.optimization_movie
        self.concatenate_xyz = self.settings.concatenate_xyz
        self.output_path = self.settings.Output_Path
        self.batches = self.settings.Batches
        # Set global paths
        self.optimization_movie_path = project_path().extend('tmp', 'opt_movie.xyz')
        self.concat_xyz_path = Path(self.output_path, 'INTEGRATION_TEST.xyz')

        # Set up files and directories
        self.initialize_assembly()

    def initialize_assembly(self):
        """
        Initializes the assembly by setting up files and directories. At the moment, this only includes deleting files from previous runs that we do not want to append to.
        """
        # Clean up files from previous runs to which we do not want to append
        if self.optimization_movie:
            self.optimization_movie_path.unlink(missing_ok=True)
        if self.concatenate_xyz:
            self.concat_xyz_path.unlink(missing_ok=True)

        return

    def run_all_batches(self):
        """
        Runs the whole assembly for all batches specified in the assembly input file.
        """
        for batch_settings in self.batches:
            batch_name, ligand_json, max_num_assembled_complexes, generate_isomer_instruction, optimisation_instruction,\
            random_seed, total_charge, metal_list, topology_similarity = self.settings.check_and_return_batch_settings(batch_settings)
            concat_xyz_path = self.concat_xyz_path if self.concatenate_xyz else None
            batch = BatchAssembly(
                                    batch_name=batch_name,
                                    ligand_json=ligand_json,
                                    max_num_assembled_complexes=max_num_assembled_complexes,
                                    generate_isomer_instruction=generate_isomer_instruction,
                                    optimisation_instruction=optimisation_instruction,
                                    random_seed=random_seed,
                                    total_charge=total_charge,
                                    metal_list=metal_list,
                                    topology_similarity=topology_similarity,
                                    output_path=self.output_path,
                                    concat_xyz_path=concat_xyz_path
                                    )
            batch.run()

        return

class BatchAssembly(object):

    def __init__(
                    self,
                    batch_name: str,
                    ligand_json: Union[str, Path],
                    max_num_assembled_complexes: int,
                    generate_isomer_instruction: str,
                    optimisation_instruction: bool,
                    random_seed: int,
                    total_charge: int,
                    metal_list: list,
                    topology_similarity: str,
                    output_path: Union[str,Path],
                    concat_xyz_path: Union[str,Path,None] = None,
                    ):
        self.output_path = Path(output_path)
        self.concat_xyz_path = Path(concat_xyz_path) if concat_xyz_path is not None else None

        self.batch_name = batch_name
        self.ligand_json = ligand_json
        self.max_num_assembled_complexes = max_num_assembled_complexes
        self.generate_isomer_instruction = generate_isomer_instruction
        self.optimisation_instruction = optimisation_instruction
        self.random_seed = random_seed
        self.total_charge = total_charge
        self.metal_list = metal_list
        self.topology_similarity = topology_similarity

    def run(self):

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
            metal_type = self.metal_list[0]
            metal_ox_state = str(self.metal_list[1])
            metal_spin = int(self.metal_list[2])

            #
            #
            # 3a. Choose Ligands Randomly
            # This section of code needs to be uncommented if you want complexes assembled randomly
            ligands = ChooseRandomLigands(database=RCA.ligand_dict,
                                          topology=Topology,
                                          instruction=Similarity,
                                          metal_oxidation_state=int(metal_ox_state),
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
            with patch('stk.SmartsFunctionalGroupFactory', new=MONKEYPATCH_STK_SmartsFunctionalGroupFactory):   # This is a monkey patch to fix an inconvenience in stk
                stk_ligand_building_blocks_list, denticities = RCA.convert_ligand_to_building_block_for_complex(
                    ligands=ligands,
                    topology=Topology,
                    metal=metal_type)

            #
            #
            # 6. Generate Isomers
            Isomers = BuildIsomers(topology=self.topology_similarity,
                                   building_blocks_list=stk_ligand_building_blocks_list,
                                   metal_input=metal_type,
                                   charge_input=metal_ox_state,
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
                complex_built = PostFilter(isomer=complex,
                                           metal_centre=metal_type,
                                           metal_offset=-0.9,  # The more negative the number the more crap gets through
                                           ligand_atom_offset=0.5,
                                           building_blocks=building_blocks,
                                           instruction=self.optimisation_instruction).closest_distance()

                output = OPTIMISE(isomer=complex_built,
                                  ligand_list=ligands,
                                  building_blocks=building_blocks,
                                  instruction=self.optimisation_instruction).Optimise_STK_Constructed_Molecule()

                complex_built = output[0]
                building_blocks_list_opt = output[1]

                complex_normal_built = PostFilter(isomer=complex_built,
                                                  metal_centre=metal_type,
                                                  building_blocks=building_blocks_list_opt,
                                                  metal_offset=-0.9,
                                                  ligand_atom_offset=0.5,
                                                  instruction=self.optimisation_instruction).post_optimisation_filter()
                Post_Process_Complex_List.append(complex_normal_built)

                if complex_normal_built is not None:
                    Sum_Assembled_Complexes += 1
            print("Leaving post process")

            #
            #
            # 8. Format Outputs

            RCA.output_controller_(list_of_complexes_wih_isomers=Post_Process_Complex_List,
                                   ligands=ligands,
                                   metal=metal_type,
                                   metal_ox_state=int(metal_ox_state),
                                   metal_multiplicity=metal_spin,
                                   view_complex=False,
                                   concatonate_xyz_name=self.concat_xyz_path,
                                   write_gaussian_input_files=False,
                                   output_directory=self.output_path, # todo: refactor to only use output_path
                                   frames=1)

            if j == 138746543956439563475683496736:
                exit()
            j += 1



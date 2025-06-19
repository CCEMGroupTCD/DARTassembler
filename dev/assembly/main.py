from DARTassembler.src.constants.paths import project_path
from DARTassembler.src.assembler.ligands import LigandChoice
import random
import yaml

from utils import *


def main(input_dict: dict):
    """
    This script executes a complete DART assembler run.
    :param: input_dict: A dictionary containing the input parameters for the DART assembler run.
    :return: None: The function does not return anything, but it generates assembled isomers based on the input parameters.
    """

    # Loop through batches
    for batch in input_dict["batches"]:
        # Create a BatchInput object
        batch = BatchInput(batch=batch)

        # Set the random seed
        random.seed(batch.random_seed)

        # Generate ligand combinations [lig1, lig2], [lig3, lig4], ...
        ligand_combinations = LigandChoice(database=[ligand.ligand_db for ligand in batch.ligands],
                                           metal_oxidation_state=batch.total_metal_oxidation_state,
                                           total_complex_charge=batch.total_charge,
                                           max_num_assembled_complexes=batch.max_num_complexes).choose_ligands()

        # Loop through ligand_combinations. Stop if max_num_complexes < idx
        for idx, ligands in enumerate(ligand_combinations):
            if batch.max_num_complexes < idx:
                break
            # Todo chuck functionality at end of IsomerFactory
            # Todo

            # ------------------------------- #
            # 1. Initial Isomer Generation
            # ------------------------------- #
            isomers_assembled = IsomerFactory.from_batch_input(batch_input=batch, ligands=ligands).generate()

            # ------------------------------- #
            # 2. Mono-ligands Optimization
            # ------------------------------- #
            isomers_mono_opt = AxialOptFactory.from_batch_input(batch_input=batch, isomers=isomers_assembled).opt_mono_rotation()

            # ------------------------------- #
            # 3. Remove redundant isomers
            # ------------------------------- #
            isomers_unique = ReduceIsomers.from_batch_input(batch_input=batch, isomers=isomers_mono_opt).reduce_isomers()

            # ------------------------------- #
            # 4. Post-assembly filter
            # ------------------------------- #
            isomers_filtered = ClashFilter.from_batch_input(batch_input=batch, isomers=isomers_unique).process_isomers()

            # ------------------------------- #
            # 5. Save the Output
            # ------------------------------- #
            # Todo: This is an incredibly Janky Save method that needs to be refactored ASAP
            test = AssemblerSave.from_batch_input(batch_input=batch, isomers=isomers_filtered)
            test.run_batch(batch_settings=batch)



    return None


if __name__ == '__main__':
    print("#===> Executing the DART assembler main script <===#\n")

    # Input file
    input_file = project_path().extend(*'dev/assembly/assembly_input/example.yaml'.split('/'))

    # Open the input file
    with open(input_file, "r") as yaml_file:
        yaml_dict = yaml.safe_load(yaml_file)

    # Execute the main function with the input dictionary
    main(input_dict=yaml_dict)

    print("\n#===> DART assembler main script finished <===#")

#########################################################################################
# This file contains the classes and methods that are used to process the input data    #
# and generate the assembled transition metal complexes                                 #
#########################################################################################
from DARTassembler.src.assembly.ligands import LigandChoice
from utilities import BatchInput, AssemblyComplex, AtomsCombiner
from ase.visualize import view
from pathlib import Path
import random
import yaml


if __name__ == "__main__":

    # Load input instructions from the input YAML file
    input_file = Path("assembly_input.yml")
    with open(input_file, "r") as yaml_file:
        yaml_dict = yaml.safe_load(yaml_file)

    # Loop through all the batches from the input file
    for batch in yaml_dict["batches"]:

        # The BatchInput object is used to store all the information about the assembly input
        assembly_input: BatchInput = BatchInput(batch)

        # The random seed is set to ensure reproducibility of the results
        random.seed(assembly_input.random_seed)

        # todo remove the list casting and implement tighter integration with the code base
        # todo in the input file we need to identify which ligands are connected and which ligands are not connected also which ligands can be swapped and which ligands can not be swapped
        # todo need to redo ligand choice to not rely on the denticity of the ligands but rather the number of vectors supplied by the user
        ligand_db_list = [ligand.ligand_db for ligand in assembly_input.ligands]
        ligand_CN_list = [ligand.temp_dent for ligand in assembly_input.ligands]
        ligand_similarity_list = [i + 1 for i in range(len(assembly_input.ligands))]  # Essentially all ligands are currently treated as "unique"

        # Now we will create the ligand choice object which will allow us to choose ligands under a number of different constraints
        ligand_choice = LigandChoice(database=ligand_db_list,
                                     metal_oxidation_state=assembly_input.total_metal_oxidation_state,
                                     total_complex_charge=assembly_input.total_charge,
                                     max_num_assembled_complexes=assembly_input.max_num_complexes)
        ligand_combinations = ligand_choice.choose_ligands()

        # Now we will loop through each of the unique ligand combinations and generate all the isomers for each combination
        all_geom = []
        for idx, ligand_combination in enumerate(ligand_combinations):

            # Create the assembly input for the current ligand combination
            ChemBuild = AssemblyComplex.from_batch_input(batch_input=assembly_input, ligands=ligand_combination)

            # Generate all the isomers for the current ligand combination
            isomers = ChemBuild.generate()


            # Add extra components
            for isomer in isomers:
                combined_struct = AtomsCombiner.from_batch_input(batch_input=assembly_input, base_atoms=isomer).combine()
                all_geom.append(combined_struct)


            if idx >= assembly_input.max_num_complexes:
                break  # If we have reached the maximum number of complexes we will break the loop

        # For the time being we will view each of the geometries for testing using a for loop and view
        a = input("Press q to quit or any other key to continue")
        for geom in all_geom:
            if a == "q":
                break
            view(geom)






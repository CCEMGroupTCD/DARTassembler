#########################################################################################
# This file contains the classes and methods that are used to process the input data    #
# and generate the assembled transition metal complexes                                 #
#########################################################################################
from DARTassembler.src.assembly.ligands import LigandChoice
from utilities import BatchInput, AssemblyComplex
from ase.visualize import view
from pathlib import Path
import random
import yaml


if __name__ == "__main__":

    # we will open the input file and read the instructions
    input_file = Path("assembly_input.yml")
    with open(input_file, "r") as yaml_file:
        yaml_dict = yaml.safe_load(yaml_file)

    # Next we should loop through all the batches and generate the instructions object
    for batch in yaml_dict["batches"]:

        # We should first create the batch input object
        assembly_input: BatchInput = BatchInput(batch)

        # Set the random seed
        random.seed(assembly_input.random_seed)

        # To make your life Timo hopefully a little easier Timo I have formatted the Assembly inputs as lists in case you don't want to use the BatchInput object
        ligand_db_list = [ligand.ligand_db.get_lig_db_in_old_format() for ligand in assembly_input.ligands]
        ligand_CN_list = [ligand.temp_dent for ligand in assembly_input.ligands] # todo this is a temporary fix. Ideally I would like to use the number of vectors supplied by the user as a proxy for the "topology" in the LigandChoice object
        ligand_target_vectors = [ligand.vectors for ligand in assembly_input.ligands]
        ligand_origin = [ligand.origin for ligand in assembly_input.ligands]
        ligand_similarity_list = [i+1 for i in range(len(assembly_input.ligands))]  # Essentially all ligands are currently treated as "unique"
        metal_type_list = [metal.metal_type for metal in assembly_input.metals]
        metal_origin_list = [metal.coord for metal in assembly_input.metals]

        # Now we will create the ligand choice object which will allow us to choose ligands under a number of different constraints
        ligand_choice = LigandChoice(database=ligand_db_list,
                                     topology=ligand_CN_list,
                                     instruction=ligand_similarity_list,
                                     metal_oxidation_state=assembly_input.total_metal_oxidation_state,
                                     total_complex_charge=assembly_input.total_charge,
                                     max_num_assembled_complexes=assembly_input.max_num_complexes)
        ligand_combinations = ligand_choice.choose_ligands()

        # Now we will loop through each of the unique ligand combinations and generate all the isomers for each combination
        all_geom = []
        for idx, ligand_combination in enumerate(ligand_combinations):
            ChemBuild = AssemblyComplex(ligands=ligand_combination,
                                        target_vectors=ligand_target_vectors,
                                        ligand_origins=ligand_origin,
                                        metal_types=metal_type_list,
                                        metal_origins=metal_origin_list,
                                        monometallic=False)
            if idx >= assembly_input.max_num_complexes:
                break   # If we have reached the maximum number of complexes we will break the loop

            isomers = ChemBuild.get_isomers()
            #isomers = ReduceIsomers(isomers, rssd_threshold=0.01).get_unique_isomers()
            all_geom.extend(isomers)

        # For the time being we will view each of the geometries for testing using a for loop and view
        for geom in all_geom:
            #input("Press Enter to view the next geometry")
            view(geom)


# TODO LIST for Cian
# 1. todo create a more simple version of the input yml file where users can specify certain cardinal axes instead of coords i.e. x instead of [1,0,0]
# 2. todo coordinate with post filters to ensure that the output is reasonable
# 3. todo test multiple geometries
# 4. todo implement isomer instruction
# 5. todo it is weird to me the Ligand choice generator does not yield the amount inputted by the user (max_num_assembled_complexes)
# 6. todo at least in the 3-3 case too many isomers are generated 4 instead of 2
# 7. todo "cis-trans" exchange isomers in the 1-1-1-1-1-1 case (come back to this)
# 8. todo what about the 2-2-1-1 mirror image case i.e. [Co(en)2Cl2]+ (https://crunchchemistry.co.uk/isomerism-in-transition-metal-complexes/)
# 9. todo optimize the TMCs




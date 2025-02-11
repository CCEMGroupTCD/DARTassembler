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
# 1. todo create a more simple version of the input yml file where users can specify certain cardinal axes instead of coords i.e. x instead of [1,0,0]  -> Timo: let's think about this sometime together since this is more part of the DART workflow, not the assembler
# 2. todo coordinate with post filters to ensure that the output is reasonable -> Timo: agree, we can just copy the code from DART, that seemed to work great. I think that should have a high priority.
# 3. todo test multiple geometries  -> Timo: agree, especially one interesting geometry would be to do a 5-1-1 geometry with a pentagonal planar ligand and two monodentate ligands, since we need to get the number of isomers correct here and to test out a non-octahedral geometry.
# 4. todo implement isomer instruction  -> Timo: I'm not sure what you mean with this, could you clarify?
# 5. todo it is weird to me the Ligand choice generator does not yield the amount inputted by the user (max_num_assembled_complexes)    -> Timo: agree, will be fixed soon (use the eff_denticities from the number of vectors supplied by the user)
# 6. todo at least in the 3-3 case too many isomers are generated 4 instead of 2 -> Timo: agree. I think your idea which checks which complex isomers are equivalent sounds great, this should have a high priority because it will solve a lot of issues.
# 7. todo "cis-trans" exchange isomers in the 1-1-1-1-1-1 case (come back to this) -> Timo: agree that this is an interesting option in DART, even though I think other projects have higher priority
# 8. todo what about the 2-2-1-1 mirror image case i.e. [Co(en)2Cl2]+ (https://crunchchemistry.co.uk/isomerism-in-transition-metal-complexes/)
# 9. todo optimize the TMCs -> Timo: unfortunately, if we want to use an MIT license as Max said before, we can use neither openbabel nor xtb. We can use rdkit, but the ff there is so bad I don't think that should be even an option. I would suggest we setup DART such that it generates reasonable starting structures and leave later optimization to the user, as we do in the OER project.

# Todo: Suggestions by Timo for other projects with high priority:
# High priority:
# 1. slight refactoring of input of the AssemblyComplex() class, since the code doesn't actually need dictionaries but lists, so it would be good to reflect this in the actual input as well. This should be a quick fix.
# 2. Automatically adjust bond lengths of the rotated ligands depending on the metal type. For now, I suggest to keep this simple, e.g. simply taking the parent metal of the ligand (e.g. Fe) and the current metal for the assembly (e.g. Ru) and then adjusting the mean bond length of the ligand by shift=(DART_Element('Fe').covalent_radius - DART_Element('Ru').covalent_radius) and the shift all the ligand atoms in the direction of the mean of the donor vectors. This would still make the ligand always rotatable around the metal center, but in the correct distance. In a later stage, we could potentially also implement more complex things here, like looking at other ligand instances with the desired metal and os, but this information is currently not in the metalig database, since the metalig only records some basic information about all ligand instances (specifically parent complex id, parent metal os, parent metal type, parent complex charge) but not all the ligand instances themselves, because that would make the file much larger.
# Medium priority:
# 2. monodentate rotater using ff (we can use the openbabel ff for now and if we do decide to get rid of it, maybe use the rdkit ff, which should be ok since it's only singlepoint calculations and only needs to be rough)


# # #=========== Idea for how to set up the assembler input ===========#
# ligands = [...]
# target_vectors = [...]
# ligand_origins = [...]
# ru1 = ase.Atom('Ru', [1, 0, 0])
# ru2 = ase.Atom('Ru', [-1, 0, 0])
# metal_centers = [  # metal centers for a bi-metallic complex with one bridging ligand
#     [ru1],
#     [ru1, ru2],
#     [ru2]
# ]
# ChemBuild = AssemblyComplex(
#     ligands=ligands,
#     target_vectors=target_vectors,
#     metal_centers=metal_centers,
#     ligand_origins=ligand_origins,
# )
# isomers = ChemBuild.get_isomers()



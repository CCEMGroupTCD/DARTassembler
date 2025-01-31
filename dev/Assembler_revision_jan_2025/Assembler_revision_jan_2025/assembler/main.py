from pathlib import Path
import yaml
from utilities import BatchInput
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from DARTassembler.src.assembly.ligands import LigandChoice
from DARTassembler.src.assembly.ligand_geometries import align_donor_atoms, assign_geometry
from DARTassembler.src.ligand_filters.Ligand_Filters import LigandFilters
from DARTassembler.src.constants.Paths import project_path



if __name__ == "__main__":

    # we will open the input file and read the instructions
    input_file = Path("input.yml")
    with open(input_file, "r") as yaml_file:
        yaml_dict = yaml.safe_load(yaml_file)

    # Next we should loop through all the batches and generate the instructions object
    for batch in yaml_dict["batches"]:
        assembly_input = BatchInput(batch)

        # We should now load the ligand database for each of the ligand sites


        # Based on charge compatability we should then choose ligands for the assembly


        # We should the align the donor atoms and assemble the final complex


        print("breakpoint")


from pathlib import Path
from DARTassembler.src.assembly.ligand_geometries import align_donor_atoms, assign_geometry
from DARTassembler.filter_ligands import filter_ligands
from DARTassembler.src.constants.Paths import project_path

import yaml
with open("input.yml", "r") as yaml_file:
    yaml_dict = yaml.safe_load(yaml_file)

print(yaml_dict)


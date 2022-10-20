from src.Molecule import RCA_Molecule
from datetime import date


class TransitionMetalComplex:
    """
    Here we define how we want the final object to look like
    to bring us in a good position for the post assembly filtering
    """

    def __init__(self, mol: RCA_Molecule, metal_symbol: str, ligands: dict, xyz_str: str):
        """
        self.molecule: basic functionality, the underlying molecule can be plotted and graphed
        self.name: the name of the ligand
        self.ligands: stores a dict of the ligands for the complex
        self.functional_groups: Quick access to the functional groups
        """
        self.molecule = mol
        self.ligands = ligands
        self.functional_groups = {key_: lig.get_assembly_dict()["type"] for key_, lig in self.ligands.items()}

        self.name = self.assemble_name(metal=metal_symbol)
        self.xyz_str = xyz_str

        self.date_of_creation = date.today()

    def assemble_name(self, metal):
        """
        this methods encodes our naming scheme
        """
        name = metal
        for i, ligand in self.ligands.items():
            try:
                name += f"_{ligand.name}"
            except AttributeError:
                # Ligand has no name assigned
                name += f"_{i}dentLig"

        return name

    def print_to_xyz(self, path):
        """
        create a xyz file for the complex at the given path
        """
        with open(f"{path}/{self.name}.xyz", "w+") as file:
            file.write(self.xyz_str)

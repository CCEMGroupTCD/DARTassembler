# todo: tbd

from src.Molecule import RCA_Molecule
from datetime import datetime


class TransitionMetalComplex:
    """
    Here we define how we want the final object to look like
    to bring us in a good position for the post assembly filtering
    """

    def __init__(self):
        self.molecule: RCA_Molecule
        self.functional_groups: dict

        self.name: str
        self.xyz_str: str

        self.date_of_creation: datetime
    def print_to_xyz(self, path):
        """
        create an xyz file for the complex at the given path
        """
        pass



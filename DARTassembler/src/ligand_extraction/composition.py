from typing import Union, List
import re

class Composition:
    """
    Class to represent a chemical composition. Similar to pymatgen's Composition class.
    """

    def __init__(self, comp: Union[str, List[str]]):
        if isinstance(comp, str):
            # Check if the string has any digits to decide the parsing method
            if any(char.isdigit() for char in comp):
                self.elements = self._parse_composition(comp)
            else:
                self.elements = self._parse_concatenated_atoms(comp)
        elif isinstance(comp, list):
            self.elements = self._parse_atom_list(comp)
        else:
            raise ValueError("Invalid input type. Only string or list of atoms is accepted.")

    def _parse_composition(self, comp_str: str) -> dict:
        """
        Parse the composition string into a dictionary of elements and their counts.

        :param comp_str: A string representation of the composition.
        :return: Dictionary of elements and their counts.
        """
        # Remove any whitespace
        clean_str = comp_str.replace(" ", "")

        # Regular expression to match element symbols (e.g., Fe, O) and their counts (e.g., 2, 3.5)
        # Elements start with a capital letter followed by 0 or more lowercase letters.
        # Counts are one or more digits possibly followed by a dot and more digits.
        pattern = re.compile(r'([A-Z][a-z]*)(\d*\.?\d*)')

        elements = {}

        for match in pattern.findall(clean_str):
            element, count = match
            # Default count to 1.0 if not provided
            count = float(count) if count else 1.0
            if element in elements:
                elements[element] += count
            else:
                elements[element] = count

        return elements

    def _parse_concatenated_atoms(self, atom_str: str) -> dict:
        """
        Parse a concatenated string of atom symbols.

        :param atom_str: A concatenated string of atom symbols.
        :return: Dictionary of elements and their counts.
        """
        # Regular expression to match element symbols without any digits
        pattern = re.compile(r'([A-Z][a-z]*)')

        elements = {}

        for match in pattern.findall(atom_str):
            if match in elements:
                elements[match] += 1.0
            else:
                elements[match] = 1.0

        return elements

    def _parse_atom_list(self, atom_list: List[str]) -> dict:
        """
        Parse the list of atoms into a dictionary of elements and their counts.

        :param atom_list: A list of atom symbols.
        :return: Dictionary of elements and their counts.
        """
        elements = {}
        for atom in atom_list:
            if atom in elements:
                elements[atom] += 1.0
            else:
                elements[atom] = 1.0

        return elements

    def almost_equals(self, other, rtol: float = 0.1, atol: float = 1e-8) -> bool:
        """
        Returns true if compositions are equal within a tolerance.

        Args:
            other (Composition): Other composition to check
            rtol (float): Relative tolerance
            atol (float): Absolute tolerance
        """
        sps = set(list(self.elements.keys()) + list(other.elements.keys()))
        # sps = set(self.elements + other.elements)
        for sp in sps:
            a = self.elements[sp] if sp in self.elements else 0
            b = other.elements[sp] if sp in other.elements else 0
            tol = atol + rtol * (abs(a) + abs(b)) / 2
            if abs(b - a) > tol:
                return False
        return True

    def get_stoichiometry(self, omit_1=False) -> str:
        """
        Returns the stoichiometry of the composition as a dictionary.
        """
        stoichiometry = ''.join([f"{element}{count:.0f}" if count != 1 or not omit_1 else element for element, count in self.elements.items()])
        return stoichiometry


if __name__ == "__main__":
    comps = ['Cl', 'CH', 'C1Br1', 'C2H5Ir1P']
    for comp in comps:
        stoi = Composition(comp).get_stoichiometry(omit_1=False)
        print(f"{comp} --> {stoi}")
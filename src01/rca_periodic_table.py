"""
This is currently not needed anywhere.
"""

from pymatgen.core.periodic_table import Element as Pymatgen_Element


class RCA_Element():

    def __init__(self, el):
        if isinstance(el, str):
            self.element = Pymatgen_Element(el)
        elif isinstance(el, int):
            self.element = Pymatgen_Element.from_Z(el)
        else:
            raise ValueError(f'Unknown element input {el}.')

        self.atomic_number = self.element.Z
        self.symbol = self.element.symbol
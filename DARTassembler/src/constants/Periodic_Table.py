"""
This file contains data of the periodic table.
It is mostly provided by mendeleev (but much faster than using mendeleev.element()) and in parts by pymatgen.
"""
# todo:
#  - add electronic ionization energies for each TMC oxidation state from NIST: https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
#  - add electronic DFT data for each electronic level for each element: https://math.nist.gov/DFTdata/atomdata/tables/ptable.html or https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7

import ast
import pandas as pd
from DARTassembler.src.constants.Paths import element_data_path
from typing import Union

# Load periodic table data from mendeleev and store it in dictionaries for fast access
ptable = pd.read_csv(element_data_path, converters={'pmg_common_oxidation_states': ast.literal_eval})
ptable_dict = ptable.set_index('symbol', drop=False).to_dict(orient='index') # dictionary of the periodic table with element symbol as key
all_atomic_symbols = list(ptable_dict.keys()) # list of all atomic symbols
atomic_number_to_symbol = ptable.set_index('atomic_number')['symbol'].to_dict() # dictionary to convert atomic number to symbol

class DART_Element(object):

    def __init__(self, el: Union[str, int]):

        self.symbol = self.get_element_symbol_from_input(el)
        self.atomic_number = ptable_dict[self.symbol]['atomic_number']
        self.Z = self.atomic_number # alias for atomic number
        self.covalent_radius = ptable_dict[self.symbol]['covalent_radius_pyykko']       # in pm
        self.covalent_radius_angstrom = self.covalent_radius / 100                      # in Angstrom
        self.atomic_mass = ptable_dict[self.symbol]['pmg_atomic_mass']
        self.common_oxidation_states = ptable_dict[self.symbol]['pmg_common_oxidation_states']
        self.is_transition_metal = ptable_dict[self.symbol]['pmg_is_transition_metal']
        self.is_metal = ptable_dict[self.symbol]['pmg_is_metal']

    @staticmethod
    def get_element_symbol_from_input(el: Union[str, int]) -> str:
        """
        Convert input to element symbol, e.g. 'H' for hydrogen.
        """
        try:
            # if el is an integer, convert it to the corresponding element symbol
            el = atomic_number_to_symbol[el]
        except KeyError:
            # if el is symbol, check if it is in the periodic table and otherwise raise an error
            try:
                ptable_dict[el]
            except KeyError:
                raise KeyError(f"Element {el} not found in periodic table")

        return str(el)






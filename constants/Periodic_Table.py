"""
This file contains the periodic table with data as provided by mendeleev, but much faster than using mendeleev.element().
"""
from mendeleev.fetch import fetch_table
from typing import Union

# Load periodic table data from mendeleev and store it in dictionaries for fast access
ptable = fetch_table("elements")
ptable_dict = ptable.set_index('symbol', drop=False).to_dict(orient='index') # dictionary of the periodic table with element symbol as key
all_atomic_symbols = list(ptable_dict.keys()) # list of all atomic symbols
atomic_number_to_symbol = ptable.set_index('atomic_number')['symbol'].to_dict() # dictionary to convert atomic number to symbol



class DART_Element(object):

    def __init__(self, el: Union[str, int]):

        self.symbol = self.get_element_symbol_from_input(el)

        self.atomic_number = ptable_dict[self.symbol]['atomic_number']
        self.covalent_radius = ptable_dict[self.symbol]['covalent_radius_pyykko']

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






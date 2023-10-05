"""
This file is for generating a table with properties of pymatgen elements.
"""

from pymatgen.core.periodic_table import Element as pmg_Element
from mendeleev.fetch import fetch_table
import pandas as pd

# Load periodic table data from mendeleev and store it in dictionaries for fast access
ptable = fetch_table("elements")
ptable_dict = ptable.set_index('symbol', drop=False).to_dict(orient='index') # dictionary of the periodic table with element symbol as key
all_elements = list(ptable_dict.keys()) # list of all atomic symbols

for el in all_elements:
    pmg_el = pmg_Element(el)
    assert pmg_el.Z == ptable_dict[el]['atomic_number']
    ptable_dict[el]['pmg_atomic_mass'] = pmg_el.atomic_mass
    ptable_dict[el]['pmg_common_oxidation_states'] = pmg_el.common_oxidation_states
    ptable_dict[el]['pmg_is_transition_metal'] = pmg_el.is_transition_metal
    ptable_dict[el]['pmg_is_metal'] = pmg_el.is_metal

ptable = pd.DataFrame.from_dict(ptable_dict, orient='index')
ptable.to_csv('../constants/element_data_table.csv')

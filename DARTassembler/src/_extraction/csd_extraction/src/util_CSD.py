"""
Some CSD API specific utility
"""


def get_molecule_from_csd(CSD_code: str, csd_reader):
    """
    searches through the csd, returning the entry with the respective csd code as a molecule object
    """
    for i, el in enumerate(csd_reader.molecules()):
        if el.identifier == CSD_code:
            return el


def get_entry_from_csd(CSD_code: str, csd_reader):
    """
    searches through the csd, returning the entry with the respective csd code as an entry object
    """
    for i, el in enumerate(csd_reader):
        if el.identifier == CSD_code:
            return el

"""
This is the script for creating a table for entries in the CSD with certain properties
"""

from ccdc import io
import pandas as pd
from tqdm import tqdm


def CSD_2_CSV(safe_to_csd: bool = True) -> pd.DataFrame:
    """
    We read through all entries of the CSD and store the desired information
    """

    csd_reader = io.EntryReader("CSD")
    d = []

    for i, entry in tqdm(enumerate(csd_reader.entries()), desc="Reading Entries from the csd"):
        d.append({
            #
            "CSD_code": entry.identifier,
            "CSD_iupac_name": entry.chemical_name,
            "date": entry.deposition_date,
            #
            # molecule properties
            "CSD_stoichiometry": entry.molecule.formula,
            "smiles": entry.molecule.smiles,
            "CSD_formal_charge": entry.molecule.formal_charge,
            "is 3d": entry.molecule.is_3d
        })

    df = pd.DataFrame(data=d)

    if safe_to_csd is True:
        df.to_csv("../data/CSD_full.csv")

    return df


if __name__ == "__main__":

    csd_df = CSD_2_CSV(safe_to_csd=True)



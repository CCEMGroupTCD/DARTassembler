from ccdc import io
from src.util_CSD import get_molecule_from_csd, get_entry_from_csd


if __name__ == "__main__":
    """
    very simple example, how you get access to a certain entry of the CSD, either as mol or entry object,
    easy setup to play around
    """

    # init the csd reader
    csd_reader = io.EntryReader("CSD")

    mol = get_molecule_from_csd(CSD_code="ABAFOZ", csd_reader=csd_reader)

    entry = get_entry_from_csd(CSD_code="ABAFOZ", csd_reader=csd_reader)


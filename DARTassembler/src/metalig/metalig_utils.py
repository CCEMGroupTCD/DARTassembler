from pathlib import Path
from typing import Union

from DARTassembler.src.constants.Paths import default_ligand_db_path, test_ligand_db_path


def get_correct_ligand_db_path_from_input(path) -> Union[Path]:
    """
    Returns either a valid path to a ligand database file or None if the path is not valid.
    """
    path = str(path)
    if path.lower() in ('', 'none', 'null', 'default', 'metalig'):
        path = default_ligand_db_path
        assert Path(path).is_file(), f"Default ligand database file '{path}' not found."
    elif path.lower() in ('test_metalig', 'test'):
        path = test_ligand_db_path
        assert Path(path).is_file(), f"Test ligand database file '{path}' not found."
    else:
        try:
            path = Path(path)
        except TypeError:
            raise ValueError(f"Invalid ligand database path string '{path}'.")
        if not path.is_file():
            raise FileNotFoundError(f"Ligand database filepath not found: '{path}'.")

    return Path(path)

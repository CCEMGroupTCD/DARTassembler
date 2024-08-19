# Ignore warnings
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

# Check if MetaLigDB exists, else uncompress it from zip
from .src.ligand_extraction.io_custom import check_if_MetaLig_exists_else_uncompress_from_zip
check_if_MetaLig_exists_else_uncompress_from_zip(delete_zip=False)

# Import all DART modules for easy access via the CLI
from .filter_ligands import filter_ligands
from .assemble_complexes import assemble_complexes
from .make_ligand_db_csv import save_dbinfo
from .concat import concatenate_ligand_databases
from .test_installation import run_installation_test
from .configs import get_default_config_files_saved

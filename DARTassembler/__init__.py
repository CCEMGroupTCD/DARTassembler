__version__ = '1.0.1'
# Ignore warnings
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

# Check if MetaLigDB exists, else uncompress it from zip
from .src.ligand_extraction.io_custom import check_if_MetaLig_exists_else_uncompress_from_zip
check_if_MetaLig_exists_else_uncompress_from_zip(delete_zip=False)

# Import all DART modules for easy access via the CLI
from .ligandfilters import ligandfilters
from .assembler import assembler
from .dbinfo import dbinfo
from .concat import concat
from .installtest import installtest
from .configs import configs

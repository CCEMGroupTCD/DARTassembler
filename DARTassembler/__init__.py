__version__ = '1.0.4'
# Ignore warnings
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

# Set openbabel environment variable BABEL_DATADIR if possible. Fixes #3.
import os
from pathlib import Path
from openbabel import __file__ as openbabel_file
from openbabel import pybel
pybel.ob.obErrorLog.StopLogging()  # Remove Openbabel warnings
openbabel_path = Path(openbabel_file).parent
ff_files = ['UFF.prm']  # forcefield parameter files
babel_datadir = None
for root, _, files in os.walk(openbabel_path):
    files_in_dir = all(file in files for file in ff_files)
    if files_in_dir:
        babel_datadir = root
        break
# If babel_datadir is not found using the above algorithm, openbabel was most likely installed by conda and hopefully already set. Should not matter though because we install it via pip, this is more for backwards compatibility with the conda installation on Mac, where the code was developed.
if babel_datadir is not None:
    os.environ['BABEL_DATADIR'] = babel_datadir

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

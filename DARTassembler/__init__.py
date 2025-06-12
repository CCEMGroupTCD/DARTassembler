__version__ = '1.0.4'
# Ignore warnings
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
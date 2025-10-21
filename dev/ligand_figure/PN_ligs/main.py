from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.metalig.mol import Ligand
from tqdm import tqdm
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['svg.fonttype'] = 'none'
import numpy as np


if __name__ == "__main__":
    db = LigandDB.from_json(path="/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/DART/examples/Pd_Ni_Cross_Coupling/generate_complexes/output/ligand_databases/P_N_ligands/ligand_db_P_N_donors.jsonlines").db

    print("done")

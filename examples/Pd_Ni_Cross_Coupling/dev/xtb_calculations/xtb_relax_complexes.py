import json

from ase.io import read, write
from ase.visualize import view
from ase.constraints import FixedPlane
import pandas as pd
from tblite.ase import TBLite
from ase.optimize import BFGS
import warnings
from tqdm import tqdm
from copy import deepcopy
from tblite.interface import Calculator
import os
import numpy as np
from ase.optimize.sciopt import SciPyFminCG
from ase.constraints import FixAtoms
from DARTassembler.src.constants.Periodic_Table import DART_Element
from pathlib import Path
from datetime import datetime

from DARTassembler.src.ligand_extraction.io_custom import save_json
from dev.src11_machine_learning.utils.utilities_ML import get_xtb_descriptors


def is_assembled_complex_dir(dir: str):
    complex_name = Path(dir).name
    structure_file = Path(dir, f'{complex_name}_structure.xyz')
    return structure_file.exists()

def relax_complex(input_file, method, force_criterion, verbosity):
    atoms = read(input_file)
    calc = TBLite(method=method, verbosity=verbosity)
    atoms.calc = calc

    opt = BFGS(atoms, trajectory=traj_file, logfile=log_file)
    opt.run(fmax=force_criterion)

    return calc, opt

def get_own_homo_lumo(atoms):
    """
    Function to calculate homo and lumo from scratch in tblite. Gives very similar but not exactly the same values as morfeus.
    """
    try:
        numbers = atoms.get_atomic_numbers()
        bohr = 1.8897259886
        positions = atoms.get_positions() * bohr # input coordinates in bohr units!!!
        xtb = Calculator("GFN2-xTB", numbers, positions, charge=0)
        xtb.set('verbosity', 0)
        res = xtb.singlepoint()

        occs = res.get('orbital-occupations')
        ens = res.get('orbital-energies')
        # df = pd.DataFrame(data=[occs, ens]).T

        homo_index = int(np.ceil(sum(occs) / 2)) - 1
        homo = ens[homo_index] * 27
        lumo_index = homo_index + 1
        lumo = ens[lumo_index] * 27
    except Exception as e:
        print(f'Complex failed homo/lumo calculation: {e}')
        homo = np.nan
        lumo = np.nan

    return homo, lumo



if __name__ == "__main__":

    assembled_complexes_dir = 'data/Pd_Ni_all_assembled_complexes'

    method = 'GFN2-xTB'
    force_criterion = 1e-3  # lower is more accurate
    verbosity = 0
    # todo: add accuracy to input




    calcdirs = [calcdir for calcdir, _, _ in os.walk(assembled_complexes_dir) if is_assembled_complex_dir(calcdir)]
    calcdirs = calcdirs#[2:4]   # TODO: DEBUGGING

    n_calcs = len(calcdirs)
    global_start = datetime.now()
    df = []
    for calcdir in tqdm(calcdirs, desc='Relaxing complexes'):
        start = datetime.now()

        cname = Path(calcdir).name
        input_file = Path(calcdir, f'{cname}_structure.xyz')

        # Set up output file paths
        complex_path = str(Path(calcdir, f'{cname}'))
        traj_file = complex_path + '_xtbtrajectory.traj'
        log_file = complex_path + '_xtblog.txt'
        relaxed_structure_file = complex_path + '_xtbstructure.xyz'
        data_file = complex_path + '_xtbdata.json'

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ResourceWarning)
            calc, opt = relax_complex(input_file, method, force_criterion, verbosity)

        atoms = opt.atoms

        metal_idx = [atom.index for atom in atoms if DART_Element(atom.symbol).is_metal]
        assert len(metal_idx) == 1, f"Expected 1 metal centre, found {len(metal_idx)}."
        metal_idx = metal_idx[0]

        # Write output relaxed structure
        relaxed_structure = opt.atoms   # todo: remove forces and charges from output structure, output only coordinates
        write(relaxed_structure_file, relaxed_structure)

        xtb_decriptors = get_xtb_descriptors(relaxed_structure_file)

        # Write output data of xtb calculation
        add_data = {'complex': cname, 'n_opt_steps': opt.nsteps, 'method': method, 'force_criterion': force_criterion, **xtb_decriptors}
        dict = {**add_data, **calc.results}
        save_json(dict, data_file)

        # Record data
        forces = np.linalg.norm(calc.results['forces'], axis=1)     # magnitudes of forces per atom
        df.append({
                    'complex': cname,
                    'duration': str(datetime.now() - start),
                    'n_opt_steps': opt.nsteps,
                    'n_atoms': len(opt.atoms),
                    'energy': calc.results['energy'],
                    'free_energy': calc.results['free_energy'],
                    'max_force': forces.max(),
                    'metal_force': forces[metal_idx],
                    'max_charge': calc.results['charges'].max(),
                    'min_charge': calc.results['charges'].min(),
                    'metal_charge': calc.results['charges'][metal_idx],
                    'dipole_x': calc.results['dipole'][0],
                    'dipole_y': calc.results['dipole'][1],
                    'dipole_z': calc.results['dipole'][2],
                    **xtb_decriptors,
                    })


    df = pd.DataFrame(df)
    overview_csv = Path(assembled_complexes_dir, 'xtb_relaxations.csv')
    df.to_csv(overview_csv, index=False)

    global_duration = datetime.now() - global_start
    print(f"Relaxation of {n_calcs} complexes took {global_duration}.")
    print('Done!')
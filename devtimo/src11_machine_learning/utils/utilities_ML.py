from pymatgen.core.periodic_table import Element as Pymatgen_Element
from io import StringIO
import tempfile
from pathlib import Path
import numpy as np

def hartree2ev(energy: float) -> float:
    return energy * 27.2114

def get_element_descriptors(el: str) -> dict:
    el = Pymatgen_Element(el)
    descriptors = {
        'Z': el.Z,
        'row': el.row,
        'group': el.group,
        'atomic_mass': el.atomic_mass,
        'atomic_radius': el.atomic_radius,
        'electron_affinity': el.electron_affinity,
        'electronegativity': el.X,
        'min_oxidation_state': el.min_oxidation_state,
        'max_oxidation_state': el.max_oxidation_state,
        'ionization_energy': el.ionization_energy
            }

    return descriptors

def get_xtb_descriptors(xyz: str, charge: int=0, n_unpaired: int=None) -> dict:
    """
    Returns a dictionary with global XTB descriptors for a given molecule.
    @param xyz: path to an xyz file or xyz string or a tuple/list with elements and coordinates
    """
    import morfeus
    from xtb.interface import XTBException

    if isinstance(xyz, list) or isinstance(xyz, tuple):
        elements, coordinates = xyz
    else:
        try:
            elements, coordinates = morfeus.read_xyz(xyz)
        except (FileNotFoundError, OSError):
            with tempfile.NamedTemporaryFile() as tmp:
                with open(tmp.name, 'w') as f:
                    f.write(xyz)
                elements, coordinates = morfeus.read_xyz(tmp.name)

    try:
        xtb = morfeus.XTB(elements, coordinates, charge=charge, n_unpaired=n_unpaired)
        dipole_x, dipole_y, dipole_z = xtb.get_dipole()
        descriptors = {
                        'ionization_potential': xtb.get_ip(),
                        'electron_affinity': xtb.get_ea(),
                        'homo': hartree2ev(xtb.get_homo()),
                        'lumo': hartree2ev(xtb.get_lumo()),
                        'dipole_x': dipole_x,
                        'dipole_y': dipole_y,
                        'dipole_z': dipole_z
                        }
    except (ModuleNotFoundError, XTBException):
        print('XTB calculation failed for some reason. Return NaNs.')
        descriptors = {
                        'ionization_potential': np.nan,
                        'electron_affinity': np.nan,
                        'HOMO': np.nan,
                        'LUMO': np.nan,
                        'Dipole_x': np.nan,
                        'Dipole_y': np.nan,
                        'Dipole_z': np.nan
                        }

    return descriptors
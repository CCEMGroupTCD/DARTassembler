import warnings
from pathlib import Path
import os
from typing import Union, Tuple, List
import numpy as np
import cclib
import ase

from DARTassembler.src.assembly.TransitionMetalComplex import TransitionMetalComplex
from DARTassembler.src.constants.Periodic_Table import DART_Element

class GaussianOutput(object):

    def __init__(self, dirpath: Union[str, Path]):
        self.dirpath = Path(dirpath).resolve()
        self.name = self.dirpath.name
        self.fchk_file = self.get_fchk_file()

        # This is the parser used for the Gaussian output file. For info how to extract which data, see https://cclib.github.io/data.html
        self.parser = self.load_gaussian_parser()

        self.relaxed_coords = self.parser.atomcoords[-1]
        self.elements = [DART_Element(Z).symbol for Z in self.parser.atomnos]

    def get_fchk_file(self) -> Path:
        fchk_files = list(self.dirpath.glob('*.fchk'))
        assert len(fchk_files) == 1, f'Found multiple fchk files in directory {self.dirpath}!'
        return Path(fchk_files[0])

    def load_gaussian_parser(self):
        return cclib.io.ccread(self.fchk_file)

    def get_gaussian_relaxed_complex(self, data_json: Union[str, Path] = None) -> TransitionMetalComplex:
        coords = self.parser.atomcoords[-1]
        elements = [DART_Element(Z).symbol for Z in self.parser.atomnos]
        complex = TransitionMetalComplex.from_json(data_json, xyz=(elements, coords))

        return complex

    def get_homo_index(self) -> int:
        homo_index = self.parser.homos[0]
        return homo_index

    def get_homo(self) -> float:
        """
        Returns the HOMO in eV.
        """
        homo_index = self.get_homo_index()
        homo_energy = self.parser.moenergies[0][homo_index]
        return homo_energy

    def get_lumo_index(self) -> int:
        """
        Returns and validates the LUMO index by checking if its energy is higher than the HOMO energy.
        """
        try:
            # Get HOMO index and its energy
            homo_index = self.get_homo_index()
            homo_energy = self.parser.moenergies[0][homo_index]

            # Start looking for the LUMO index
            lumo_index_candidate = homo_index + 1  # The orbital just above the HOMO as a starting point

            # Loop to find the first orbital with energy greater than HOMO
            mo_energies = self.parser.moenergies[0]
            for i in range(lumo_index_candidate, len(mo_energies)):
                if mo_energies[i] > homo_energy:
                    if i > lumo_index_candidate:
                        warnings.warn(f'Degenerate orbitals at HOMO level recognized. The returned LUMO index is not the one after the HOMO index in the list, but the first one which has a higher energy.')
                    return i  # Found the valid LUMO index

            # If loop finishes without finding a valid LUMO
            warnings.warn(f"No valid LUMO found in calculation {self.dirpath}")
            return None

        except AttributeError as e:
            print(f'Encountered error {e} in calculation {self.dirpath}')
            return None

    def get_lumo(self) -> float:
        """
        Returns the LUMO in eV.
        """
        lumo_index = self.get_lumo_index()
        if lumo_index is not None:
            lumo_energy = self.parser.moenergies[0][lumo_index]
        else:
            lumo_energy = np.nan
        return lumo_energy

    def get_homo_lumo_hlgap(self) -> Tuple[float, float, float]:
        """
        Returns the HOMO, LUMO and HOMO-LUMO gap in eV.
        """
        try:
            # Get HOMO and LUMO energies
            homo_energy = self.get_homo()  # in eV
            lumo_energy = self.get_lumo()  # in eV

            # Calculate HOMO-LUMO gap
            gap = lumo_energy - homo_energy

            eps = 1e-4
            if abs(gap) < eps:
                warnings.warn(f'Issue with LUMO recognition: The recognized LUMO is equal to the HOMO. This is probably due to degenerate orbitals at the HOMO level. Please check your ')

        except AttributeError as e:
            print(f'Encountered error {e} in calculation {self.dirpath}')
            homo_energy, lumo_energy, gap = np.nan, np.nan, np.nan

        return homo_energy, lumo_energy, gap

    def get_metal_charge(self, method='mulliken'):
        # Make sure charges have been calculated; you may use 'mulliken' or other methods depending on your calculation
        try:
            if method in self.parser.atomcharges:
                mulliken_charges = self.parser.atomcharges[method]

                # Get the index of the metal atom
                atom_numbers = self.parser.atomnos
                metal_index = [idx for idx, Z in enumerate(atom_numbers) if DART_Element(Z).is_transition_metal]
                assert len(metal_index) == 1, f'Multiple transition metal indices found in complex: {metal_index}!'
                metal_index = metal_index[0]

                # Get the charge on the metal atom
                metal_charge = mulliken_charges[metal_index]
            else:
                raise ValueError(f'{method.capitalize()} charges have not been found.')
        except AttributeError as e:
            print(f'Encountered error {e} in calculation {self.dirpath}')
            metal_charge = np.nan

        return metal_charge

    @staticmethod
    def is_finished_calc_dir(calc_path):
        """
        Checks if the directory contains a calculation.
        """
        files = os.listdir(calc_path)
        if any([file.endswith(".com") for file in files]) and \
                any([file.endswith(".log") for file in files]) and \
                any([file.endswith(".chk") for file in files]):
            return True
        else:
            return False

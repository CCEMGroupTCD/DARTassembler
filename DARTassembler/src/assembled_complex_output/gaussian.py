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

    def __init__(self, input_path: Union[str, Path]):
        self.input_path = Path(input_path).resolve()
        self.fchk_file = self._get_fchk_file(self.input_path)
        self.dirpath = self.fchk_file.parent
        self.name = self.dirpath.name

        # This is the parser used for the Gaussian output file. For info how to extract which data, see https://cclib.github.io/data.html
        self.parser = self.load_gaussian_parser()

        self.relaxed_coords = self.parser.atomcoords[-1]
        self.elements = [DART_Element(Z).symbol for Z in self.parser.atomnos]


    def _get_fchk_file(self, input_path: Path) -> Path:
        """
        Returns the fchk file from the input path. If the input path is a directory, it will look for the fchk file in that directory. If the input path is a file, it will return that file.
        :param input_path: Path to the fchk file or directory containing the fchk file.
        :return: Path to the fchk file.
        """
        if not input_path.exists():
            raise FileNotFoundError(f'Path {input_path} not found.')

        if input_path.is_file():
            fchk_file = input_path
        else:
            dirpath = input_path
            fchk_files = list(dirpath.glob('*.fchk'))
            if len(fchk_files) != 1:
                raise ValueError(f'Found multiple .fchk files in directory {dirpath}. Please specify the .fchk file directly.')
            fchk_file = Path(fchk_files[0]).resolve()

        return fchk_file

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
                warnings.warn(f'Issue with LUMO recognition: The recognized LUMO is equal to the HOMO. This is probably due to degenerate orbitals at the HOMO level. Please check your input.')

        except AttributeError as e:
            print(f'Encountered error {e} in calculation {self.dirpath}')
            homo_energy, lumo_energy, gap = np.nan, np.nan, np.nan

        return homo_energy, lumo_energy, gap

    def get_total_energy(self) -> float:
        """
        Returns the total energy in eV.
        :return: total energy (eV)
        """
        energies = self.parser.scfenergies

        # Take the last energy value provided. There are several possibilities how multiple scf energies can end up in the parser, especially when parsing the .log file instead of the .fchk file. For example, it's possible to get 2 values in a single point calculation, if you specify `scf=xqc`. In this case, if Gaussian didn't converge after 33 iterations, it will change the optimizer to xqc and before that, print out the energy that it's at at this moment. Also if Gaussian did a geometry optimization, all the singlepoint energies will be in this attribute and taking the last one is important to get the minimum energy.
        total_energy = energies[-1]  # in eV

        # Usually, the last energy should also be the minimum energy. Print a warning if this is not the case.
        if len(energies) > 1:
            if total_energy > min(energies):
                warnings.warn(f'Last energy in the list of SCF energies is not the minimum energy. This may be due to a non-converged calculation. Please check your input for input file {self.input_path}.')

        return total_energy

    def get_metal_charge(self, method='mulliken'):
        # Make sure charges have been calculated; you may use 'mulliken' or other methods depending on your calculation
        try:
            if method in self.parser.atomcharges:
                charges = self.parser.atomcharges[method]

                # Get the index of the metal atom
                atom_numbers = self.parser.atomnos
                metal_index = [idx for idx, Z in enumerate(atom_numbers) if DART_Element(Z).is_transition_metal]
                assert len(metal_index) == 1, f'Multiple transition metal indices found in complex: {metal_index}!'
                metal_index = metal_index[0]

                # Get the charge on the metal atom
                metal_charge = charges[metal_index]
            else:
                raise ValueError(f'{method.capitalize()} charges have not been found.')
        except AttributeError as e:
            print(f'Encountered error {e} in calculation {self.dirpath}')
            metal_charge = np.nan

        return metal_charge

    def get_all_raw_data(self) -> dict:
        """
        Returns all raw data from the Gaussian output file.
        :returns: Dictionary containing all raw data from the Gaussian output file.
        """
        return vars(self.parser)

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



if __name__ == '__main__':

    # Test the GaussianOutput class
    fchk_path = '/Users/timosommer/PhD/projects/OERdatabase/dev/timo/playground/complex_data/dftsp_output/HOQOKUME_Mn_OH/HOQOKUME_Mn_OH_structure_gfnffrot_gfn2opths_dftspls/HOQOKUME_Mn_OH_structure_gfnffrot_gfn2opths_dftspls.log'  # single point calculation
    # fchk_path = '/Users/timosommer/PhD/projects/RCA/projects/DART/examples/Pd_Ni_Cross_Coupling/dev/output/gaussian_relaxed_complexes/batches/P_N_Donors_Ni_Metal_Centre/complexes/ABADEZIX_PN_Ni/ABADEZIX_PN_Ni_gaussian.log'  # relaxed complex

    gaussian = GaussianOutput(fchk_path)
    homo, lumo, hlgap = gaussian.get_homo_lumo_hlgap()
    metal_charge = gaussian.get_metal_charge()
    data = gaussian.get_all_raw_data()

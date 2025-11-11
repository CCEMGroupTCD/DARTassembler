import warnings
from shutil import rmtree
from typing import Union
import numpy as np
import pandas as pd
from pathlib import Path
import ase
from ase.io import read
import filecmp
import os
from DARTassembler.src.constants.paths import project_path
try:
    import pytest
except ImportError:
    pass

class IntegrationTest(object):
    def __init__(self, new_dir, old_dir, xyz_tol=0):
        """
        Compares two directories and prints the differences. Run compare_all() to get the results.
        """
        self.new_dir = Path(new_dir)
        self.old_dir = Path(old_dir)
        self.only_in_new = []
        self.only_in_old = []
        self.changed = []
        self.small_changes = {} # Dictionary with files that have small changes, where the key is the file name and the value is a string describing the change.
        self.dcmp = filecmp.dircmp(str(self.new_dir), str(self.old_dir))

        # Settings
        self.xyz_tol = xyz_tol

    def compare_all(self) -> Union[bool, list[str]]:
        """
        Compares two directories and prints the differences.
        :return: successful (bool), small_but_relevant_changes (list of str). If successful is True, the directories are the same. If False, small_but_relevant_changes contains the names of files with small but relevant changes.
        """
        print('\nIntegration test: check if the new output is the same as the old output.')
        self._only_in_one(self.dcmp)
        self._check_file_changes(self.dcmp)
        self._check_for_only_small_changes()
        self.only_in_new.sort()
        self.only_in_old.sort()
        self.changed.sort()
        self.dirs_only_in_new = [item for item in self.only_in_new if Path(self.new_dir, item).is_dir()]
        self.files_only_in_new = [item for item in self.only_in_new if Path(self.new_dir, item).is_file()]
        self.dirs_only_in_old = [item for item in self.only_in_old if Path(self.old_dir, item).is_dir()]
        self.files_only_in_old = [item for item in self.only_in_old if Path(self.old_dir, item).is_file()]
        self.successful = not self.only_in_new and not self.only_in_old and not self.changed
        self._print_results()

        # # Raise an error if any concat_....xyz files have more than small changes.
        small_but_relevant_changes = []
        for file in self.changed:
            filename = Path(file).name
            if filename.startswith('concat_') and filename.endswith('.xyz'):
                if file in self.small_changes:
                    message = self.small_changes[file]  # Example message: '-> Same: 27/27. Diff. el: 0/27. Diff. ...'
                    n_same, n_total = message.split('Same: ')[1].split('.')[0].split('/')
                    if n_same != n_total:
                        small_but_relevant_changes.append(filename)
                else:
                    small_but_relevant_changes.append(filename)

        return self.successful, small_but_relevant_changes

    def _check_for_only_small_changes(self):
        for file in self.changed:
            change = ''
            if file.endswith('.csv'):
                change = self._compare_csv_files(Path(self.new_dir, file), Path(self.old_dir, file))
            elif file.endswith('.xyz'):
                change = self._compare_xyz_files(Path(self.new_dir, file), Path(self.old_dir, file))

            if change is None:
                self.small_changes[file] = 'Warning: difference of files expected, but none found'
            elif change != '':
                self.small_changes[file] = change

    def _compare_xyz_files(self, new_file, old_file) -> Union[str]:
        """
        Compares two xyz files and return a string describing the differences.
        """
        xyz_test = XYZIntegrationTest(new_xyz=new_file, old_xzy=old_file, tol=self.xyz_tol)
        return xyz_test.compare_and_return_result_string()

    def _compare_csv_files(self, new_file, old_file) -> Union[str,None]:
        """
        Compares two csv files and prints the differences.
        """
        try:
            df_new = pd.read_csv(new_file)
        except pd.errors.EmptyDataError:
            df_new = pd.DataFrame()
        try:
            df_old = pd.read_csv(old_file)
        except pd.errors.EmptyDataError:
            df_old = pd.DataFrame()

        change = None
        try:
            pd.testing.assert_frame_equal(df_new, df_old)
            change = 'numerical'
        except AssertionError:
            try:
                # Ignore order of rows and columns
                df_new2 = df_new.sort_values(by=list(df_new.columns), axis=0).reset_index(drop=True)
                df_old2 = df_old.sort_values(by=list(df_old.columns), axis=0).reset_index(drop=True)
                pd.testing.assert_frame_equal(df_new2, df_old2, check_like=True)
                change = 'change in order of rows/columns'
            except AssertionError:
                change = 'other'

        return change

    def _only_in_one(self, dcmp):
        self.only_in_new += [str(Path(dcmp.left) / file).replace(str(self.new_dir), '.', 1) for file in dcmp.left_only]
        self.only_in_old += [str(Path(dcmp.right) / file).replace(str(self.old_dir), '.', 1) for file in dcmp.right_only]
        for sub_dcmp in dcmp.subdirs.values():
            self._only_in_one(sub_dcmp)

    def _check_file_changes(self, dcmp):
        self.changed += [str(Path(dcmp.left) / file).replace(str(self.new_dir), '.', 1) for file in dcmp.diff_files]
        for sub_dcmp in dcmp.subdirs.values():
            self._check_file_changes(sub_dcmp)

    def _print_results(self):
        self._print_category('Only in new', self.only_in_new, self.new_dir.name)
        self._print_category('Only in old', self.only_in_old, self.old_dir.name)
        self._print_category('Changed', self.changed)
        self._print_stats()

    def _print_stats(self):
        if self.successful:
            print('Integration test successful: all good!')
            return

        print('\n==========    WARNING: INTEGRATION TEST FOUND ISSUES    ==========')
        if self.dirs_only_in_new:
            print(f'\t# directories missing in old: {len(self.dirs_only_in_new)}')
        if self.files_only_in_new:
            print(f'\t# files missing in old: {len(self.files_only_in_new)}')
        if self.dirs_only_in_old:
            print(f'\t# directories missing in new: {len(self.dirs_only_in_old)}')
        if self.files_only_in_old:
            print(f'\t# files missing in new: {len(self.files_only_in_old)}')
        if self.changed:
            print(f'\t# changed files: {len(self.changed)}')

    def _get_small_change_of_changed_file(self, file) -> str:
        if file in self.small_changes:
            small_change = f'\t\t-> {self.small_changes[file]}'
        else:
            small_change = ''
        return small_change

    def _print_category(self, label, files, dirname=None):
        if not files:
            return
        dirname = ' ' + dirname + os.sep if dirname else ''
        print(f'\n{label}{dirname}: {len(files)} item(s)')
        last_path_parts = files[0].split(os.sep)
        for i, file in enumerate(files):
            file_parts = file.split(os.sep)
            output = ' '
            if i == 0:  # Always print the full path for the first file
                output = '  ' + file.lstrip('.') # Remove leading dot
                output += self._get_small_change_of_changed_file(file)
                print(output)
                last_path_parts = file_parts
                continue
            for j, (last_part, current_part) in enumerate(zip(last_path_parts, file_parts)):
                if last_part == current_part and j != len(file_parts) - 1:
                    output += ' ' * len(last_part)
                else:
                    output += os.sep.join(file_parts[j:])
                    break
                output += os.sep
            output = output.replace(f' {os.sep} ', '   ')
            output += self._get_small_change_of_changed_file(file)
            print(output)
            last_path_parts = file_parts




class XYZIntegrationTest(object):
    def __init__(self, old_xzy, new_xyz, tol=1e-5):
        self.old_xyz = old_xzy
        self.new_xyz = new_xyz
        self.tol = tol
        self.old_mols = self.read_concatenated_xyz_file(old_xzy)
        self.new_mols = self.read_concatenated_xyz_file(new_xyz)

    def read_concatenated_xyz_file(self, filename):
        try:
            mols = read(filename, index=':')
        except ase.io.formats.UnknownFileTypeError:
            # The file is empty
            mols = []
        return mols

    def get_molecules_with_different_interatomic_distances(self):
        mols_with_different_interatomic_distances = []
        for idx, (mol1, mol2) in enumerate(zip(self.old_mols, self.new_mols)):
            dist = np.diagonal(ase.geometry.get_distances(mol1.positions, mol2.positions)[1])
            if np.any(dist > self.tol):
                mols_with_different_interatomic_distances.append(idx)
        return mols_with_different_interatomic_distances

    def compare_individual_xyz_files(self, old_mol, new_mol):
        mol_results = {}
        # Check number of atoms
        same_n_atoms = len(old_mol) == len(new_mol)
        mol_results['same_n_atoms'] = same_n_atoms
        # Check atom types
        mol_results['diff_el'] = np.any(old_mol.get_chemical_symbols() != new_mol.get_chemical_symbols())
        mol_results['diff_el_order'] = np.any(sorted(old_mol.get_chemical_symbols()) != sorted(new_mol.get_chemical_symbols()))
        if same_n_atoms:
            # Check interatomic distances
            dist = np.diagonal(ase.geometry.get_distances(old_mol.positions, new_mol.positions)[1])
            mol_results['sum_interatomic_distances'] = dist.sum()
            mol_results['n_diff_interatomic_distances'] = np.sum(dist > self.tol)
            # Check xyz coordinates
            is_diff_coord = ~np.isclose(old_mol.positions, new_mol.positions, atol=self.tol, rtol=0)
            mol_results['n_diff_xyz_coordinates'] = np.sum(is_diff_coord)
            mol_results['sum_diff_xyz_coordinates'] = np.sum(np.abs(old_mol.positions - new_mol.positions)[is_diff_coord])
        else:
            mol_results['sum_interatomic_distances'] = False
            mol_results['n_diff_interatomic_distances'] = False
            mol_results['n_diff_xyz_coordinates'] = False
            mol_results['sum_diff_xyz_coordinates'] = False

        return mol_results

    def compare_xyz_files(self, print=True):
        assert (len(self.old_mols) == len(self.new_mols))

        df_mol_results = []
        for idx, (old_mol, new_mol) in enumerate(zip(self.old_mols, self.new_mols)):
            mol_results = self.compare_individual_xyz_files(old_mol, new_mol)
            mol_results['Index'] = idx
            df_mol_results.append(mol_results)

        df_mol_results = pd.DataFrame(df_mol_results)

        if print:
            print('=================    Summary:    =================')
            n_changed_atom_numbers = (df_mol_results['n_diff_atom_types'] != 0).sum()
            print(f"Number of molecules with changed atom numbers: {n_changed_atom_numbers}")
            if n_changed_atom_numbers > 0:
                print(f"\t--> Indices of molecules: {df_mol_results[df_mol_results['n_diff_atom_types'] != 0]['Index'].values}")
            n_changed_positions = (df_mol_results['n_diff_xyz_coordinates'] != 0).sum()
            print(f"Number of molecules with changed positions: {n_changed_positions}")
            if n_changed_positions > 0:
                print(f"\t--> Indices of molecules: {df_mol_results[df_mol_results['n_diff_xyz_coordinates'] != 0]['Index'].values}")
            n_changed_interatomic_distances = (df_mol_results['n_diff_interatomic_distances'] != 0).sum()
            print(f"Number of molecules with changed interatomic distances: {n_changed_interatomic_distances}")
            if n_changed_interatomic_distances > 0:
                print(f"\t--> Indices of molecules: {df_mol_results[df_mol_results['n_diff_interatomic_distances'] != 0]['Index'].values}")

            if n_changed_positions > 0 or n_changed_atom_numbers > 0:
                print("!!!Fatal Error!!! --> The code is not working as expected")
            else:
                print("!!!Success!!! --> The code still works --> you can pat yourself on the back")

        return df_mol_results

    def _old_mol_identical_to_one_mol_in_new_mols(self) -> list[bool]:
        """
        Return a Boolean list saying whether every molecule in `self.old_mols`
        occurs (ignoring order) in `self.new_mols`, within `self.tol`.
        :return: list of bools of len `self.old_mols`, where each bool indicates whether the corresponding molecule in `self.old_mols` is the same as any molecule in `self.new_mols`.
        """
        def _fingerprint(mol, tol):
            # 1. composition – order independent
            comp = ''.join(sorted(mol.get_chemical_symbols()))  # e.g. 'CCHHO'
            # 2. geometry signature – sort + round distances to the tol grid
            d = mol.get_all_distances(mic=False).ravel()
            sig = tuple(np.sort(np.rint(d / tol).astype(np.int16)))  # tuple → hashable
            return (comp, sig)

        # Pre-fingerprint all new molecules once
        new_signatures = {_fingerprint(m, self.tol) for m in self.new_mols}

        # Check each old molecule against the set
        is_same = [_fingerprint(m, self.tol) in new_signatures for m in self.old_mols]
        return is_same

    def compare_and_return_result_string(self) -> Union[str, None]:
        """
        Compares the two xyz files and returns a short string describing the differences. If there are no differences, None is returned.
        :return: str or None
        """
        matches = self._old_mol_identical_to_one_mol_in_new_mols()
        n = max(len(self.new_mols), len(self.old_mols))
        result_string = f'Same: {sum(matches)}/{n}. '

        same_n_mols = (len(self.old_mols) == len(self.new_mols))
        if not same_n_mols:
            result_string += f'Diff. n mol: old: {len(self.old_mols)}, new: {len(self.new_mols)}'
            return result_string

        df_mol_results = self.compare_xyz_files(print=False)
        n_same_atoms = df_mol_results['same_n_atoms'].sum()
        n_diff_el = df_mol_results['diff_el'].sum()
        n_diff_el_order = df_mol_results['diff_el_order'].sum()
        sum_interatomic_distances = df_mol_results['sum_interatomic_distances'].mean()
        n_diff_interatomic_distances = df_mol_results['n_diff_interatomic_distances'].mean().round().astype(int)
        n_diff_xyz_coordinates = df_mol_results['n_diff_xyz_coordinates'].mean().round().astype(int)
        sum_diff_xyz_coordinates = df_mol_results['sum_diff_xyz_coordinates'].mean()

        if n_same_atoms != n:
            string = f' for {n - n_same_atoms}/{n} mols' if n > 1 else ''
            result_string += f'Diff. n atoms: {n_same_atoms}/{n}{string}. '
            return result_string

        result_string += f'Diff. el: {n_diff_el}/{n}. '
        result_string += f'Diff. el order: {n_diff_el_order}/{n}. '
        result_string += f'Diff. dist: {n_diff_interatomic_distances} (sum={sum_interatomic_distances:.2g}A). '
        result_string += f'Diff. coords: {n_diff_xyz_coordinates} (sum={sum_diff_xyz_coordinates:.2g}A). '

        if result_string == '':
            result_string = None

        return result_string



if __name__ == '__main__':
    # file1name = str(project_path().extend("src14_Assembly_Unit_Test", "INTEGRATION_TEST_Benchmark_Timo.xyz"))
    # file2name = str(project_path().extend("src14_Assembly_Unit_Test", "INTEGRATION_TEST.xyz"))
    # allowed_differences = 1e-5
    # # df_xyz_diff = compare_xyz_files(file1name, file2name, allowed_differences)
    # df_mol_diff = XYZIntegrationTest(file1name, file2name, tol=allowed_differences).compare_xyz_files()
    path = project_path().extend('src14_Assembly_Unit_Test')
    test = IntegrationTest(new_dir=Path(path, '../../../src14_Assembly_Unit_Test/output'), old_dir=Path(path,
                                                                                                        '../../../src14_Assembly_Unit_Test/output_benchmark'))
    test.compare_all()

    xyz_test = XYZIntegrationTest(Path(path, 'output/INTEGRATION_TEST.xyz'), Path(path, 'output_benchmark/INTEGRATION_TEST.xyz'), tol=1e-2)
    xyz_test.compare_xyz_files()
    print('Done!')


def integration_test(name: str):
    """
    Decorator to wrap an assembler integration test.
    The wrapped function must accept exactly one argument: `outdir` (Path),
    and should perform the run (e.g. Assembler.run_from_yaml(...)).
    It should return the `assembly` object that has `output_directory`.
    """
    def _decorate(run_fn):
        @pytest.mark.usefixtures()
        def _wrapped(change):
            outdir = project_path().extend('tests', 'pytest', name, 'data_output').resolve()
            old_dir = outdir.parent / 'benchmark_data_output'

            cwd = Path.cwd()
            # Assert that the outdir is a subdirectory of the cwd, otherwise we might be in the wrong project
            assert cwd in outdir.resolve().parents, f'Outdir {outdir} is not a subdirectory of the project path {project_path()}. Are you running the test in the correct project?'

            # Fresh start
            if outdir.exists():
                rmtree(outdir)
            try:
                outdir.mkdir(parents=True, exist_ok=True)
                os.chdir(outdir)
                results = run_fn(outdir)   # user-supplied function does the actual run
            except Exception as e:
                pytest.fail(f'Exception while testing {name}: {e}')
            finally:    # Change back to the original working directory
                os.chdir(cwd)

            # Compare with benchmark, set "changed" flag for custom reporting if needed
            if old_dir.exists():
                test = IntegrationTest(new_dir=outdir, old_dir=old_dir)
                successful, small_but_relevant_changes = test.compare_all()
                if small_but_relevant_changes:
                    change.big('Molecules changed majorly: ' + ', '.join(small_but_relevant_changes))
                elif not successful:
                    change.small('Same molecules, just small changes in coordinates and/or other files/things.')
            else:
                pytest.skip(f'Could not find benchmark folder "{old_dir}"!')

            return results
        return _wrapped
    return _decorate

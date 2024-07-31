import warnings
from DARTassembler.src.constants.Paths import project_path
import numpy as np
import pandas as pd
from pathlib import Path
import ase
from ase.io import read
import filecmp
import os

class IntegrationTest(object):
    def __init__(self, new_dir, old_dir):
        """
        Compares two directories and prints the differences. Run compare_all() to get the results.
        """
        self.new_dir = Path(new_dir)
        self.old_dir = Path(old_dir)
        self.only_in_new = []
        self.only_in_old = []
        self.changed = []
        self.dcmp = filecmp.dircmp(str(self.new_dir), str(self.old_dir))

    def compare_all(self):
        """
        Compares two directories and prints the differences.
        """
        print('\nIntegration test: check if the new output is the same as the old output.')
        self._only_in_one(self.dcmp)
        self._check_file_changes(self.dcmp)
        self.only_in_new.sort()
        self.only_in_old.sort()
        self.changed.sort()
        self.dirs_only_in_new = [item for item in self.only_in_new if Path(self.new_dir, item).is_dir()]
        self.files_only_in_new = [item for item in self.only_in_new if Path(self.new_dir, item).is_file()]
        self.dirs_only_in_old = [item for item in self.only_in_old if Path(self.old_dir, item).is_dir()]
        self.files_only_in_old = [item for item in self.only_in_old if Path(self.old_dir, item).is_file()]
        self._print_results()
        # self.compare_all_csv_files()

    def compare_all_csv_files(self):
        for file in self.changed:
            if file.endswith('.csv'):
                self._compare_csv_files(Path(self.new_dir, file), Path(self.old_dir, file))
    def _compare_csv_files(self, new_file, old_file, exact=True):
        """
        Compares two csv files and prints the differences.
        """
        if exact:
            assert filecmp.cmp(new_file, old_file, shallow=False), f'csv files are not the same: \n - {new_file}\n - {old_file}'
        else:
            df_new = pd.read_csv(new_file)
            df_old = pd.read_csv(old_file)
            pd.testing.assert_frame_equal(df_new, df_old, check_like=True)


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
        if not self.only_in_new and not self.only_in_old and not self.changed:
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
                print('  ' + file.lstrip('.'))  # Remove leading dot
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
            print(output)
            last_path_parts = file_parts




class XYZIntegrationTest(object):
    def __init__(self, old_xzy, new_xyz, tol=1e-5):
        self.old_xyz = old_xzy
        self.new_xyz = new_xyz
        self.tol = tol
        self.old_mols = read(old_xzy, index=':')
        self.new_mols = read(new_xyz, index=':')

    def get_molecules_with_different_interatomic_distances(self):
        mols_with_different_interatomic_distances = []
        for idx, (mol1, mol2) in enumerate(zip(self.old_mols, self.new_mols)):
            dist = np.diagonal(ase.geometry.get_distances(mol1.positions, mol2.positions)[1])
            if np.any(dist > self.tol):
                mols_with_different_interatomic_distances.append(idx)
        return mols_with_different_interatomic_distances

    def compare_xyz_files(self):
        if len(self.old_mols) != len(self.new_mols):
            warnings.warn(f"Number of molecules is not the same in the two files. old: {len(self.old_mols)}, new {len(self.new_mols)}. Try to compare anyway.")

        df_mol_results = []
        for idx, (old_mol, new_mol) in enumerate(zip(self.old_mols, self.new_mols)):
            mol_results = {'Index': idx}
            # Check number of atoms
            same_n_atoms = len(old_mol) == len(new_mol)
            mol_results['same_n_atoms'] = same_n_atoms
            # Check atom types
            mol_results['n_diff_atom_types'] = np.sum(old_mol.get_chemical_symbols() != new_mol.get_chemical_symbols())
            if same_n_atoms:
            # Check interatomic distances
                dist = np.diagonal(ase.geometry.get_distances(old_mol.positions, new_mol.positions)[1])
                mol_results['sum_interatomic_distances'] = dist.sum()
                mol_results['n_diff_interatomic_distances'] = np.sum(dist > self.tol)
                # Check xyz coordinates
                mol_results['n_diff_xyz_coordinates'] = np.sum(~np.isclose(old_mol.positions, new_mol.positions, atol=self.tol))
            else:
                mol_results['sum_interatomic_distances'] = False
                mol_results['n_diff_interatomic_distances'] = False
                mol_results['n_diff_xyz_coordinates'] = False
            df_mol_results.append(mol_results)
        df_mol_results = pd.DataFrame(df_mol_results)

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
import warnings

from constants.Paths import project_path
from itertools import zip_longest
import numpy as np
import pandas as pd
import ase
from ase.io import read

def compare_xyz_files(file1name, file2name, allowed_differences=0.1):
    df_test = []
    n_complex = -1
    with open(file1name) as file1, open(file2name) as file2:
        i = 1
        for line1, line2 in zip_longest(file1, file2):
            line1 = str(line1).split()
            line2 = str(line2).split()
            if len(line1) == len(line2):
                if len(line1) == 4:
                    p0 = np.array([float(line1[-3]), float(line1[-2]), float(line1[-1])])
                    p1 = np.array([float(line2[-3]), float(line2[-2]), float(line2[-1])])
                    dist = np.linalg.norm(p0 - p1)
                    df_test.append({
                                        'Line': i,
                                        'Complex': n_complex,
                                        'Benchmark': line1,
                                        'Test': line2,
                                        'same_positions': dist <= allowed_differences,
                                        'same_atom_types': str(line1[-4]) == str(line2[-4]),
                                        'distance': dist,
                                        })
                    if str(line1[-4]) != str(line2[-4]):
                        print(f"!!!Fatal Error!!! --> The atom types have changed [{str(line1[-4])}] != [{str(line2[-4])}] on line [{i}] --> The code is not working as expected")
                        # raise AssertionError


                    elif dist > 0.1:
                        print(f"!!!Fatal Error!!! --> The atom positions have changed [{p0}] != [{p1}] on line [{i}] --> The code is not working as expected")
                        # raise AssertionError
                    else:
                        pass
                elif len(line1) == 1:
                    n_complex += 1
                else:
                    pass
            else:
                print(f"!!!Fatal Error!!! --> The length of the following lines are not equivalent --> The code is not working as expected")
                print(f'Line {i}: File1: {line1}, File2: {line2}')
                raise AssertionError
            i += 1
    df_test = pd.DataFrame(df_test)

    # Check if interatomic distances are the same
    mols1 = read(file1name, index=':')
    mols2 = read(file2name, index=':')
    mols_with_different_interatomic_distances = []
    for idx, (mol1, mol2) in enumerate(zip(mols1, mols2)):
        dist = np.diagonal(ase.geometry.get_distances(mol1.positions, mol2.positions)[1])
        if np.any(dist > 1e-5):
            mols_with_different_interatomic_distances.append(idx)

    print('\n=================    Summary:    =================')
    n_changed_atoms = (~df_test['same_atom_types']).sum()
    n_changed_positions = (~df_test['same_positions']).sum()
    print(f"Number of changed atoms: {n_changed_atoms}")
    print(f"Number of changed positions: {n_changed_positions}")
    if len(mols_with_different_interatomic_distances) > 0:
        print(f"Interatomic distances have changed for {len(mols_with_different_interatomic_distances)} molecules: {mols_with_different_interatomic_distances}!!!")
    else:
        print("All interatomic distances are the same.")
    if n_changed_positions > 0 or n_changed_atoms > 0:
        print("!!!Fatal Error!!! --> The code is not working as expected")
    else:
        print("!!!Success!!! --> The code still works --> you can pat yourself on the back")

    return df_test


class AssemblyIntegrationTest(object):
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

        print('\n=================    Summary:    =================')
        n_changed_atom_numbers = (df_mol_results['n_diff_atom_types'] != 0).sum()
        print(f"Number of molecules with changed atom numbers: {n_changed_atom_numbers}")
        if n_changed_atom_numbers > 0:
            print(f"\t--> Indices of molecules: {df_mol_results[df_mol_results['n_diff_atom_types'] != 0]['Index'].values}")
        n_changed_positions = (df_mol_results['n_diff_xyz_coordinates'] != 0).sum()
        if n_changed_positions > 0:
            print(f"\t--> Indices of molecules: {df_mol_results[df_mol_results['n_diff_xyz_coordinates'] != 0]['Index'].values}")
        print(f"Number of molecules with changed positions: {n_changed_positions}")
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
    file1name = str(project_path().extend("src14_Assembly_Unit_Test", "INTEGRATION_TEST_Benchmark_Timo.xyz"))
    file2name = str(project_path().extend("src14_Assembly_Unit_Test", "INTEGRATION_TEST.xyz"))
    allowed_differences = 1e-5
    # df_xyz_diff = compare_xyz_files(file1name, file2name, allowed_differences)
    df_mol_diff = AssemblyIntegrationTest(file1name, file2name, tol=allowed_differences).compare_xyz_files()

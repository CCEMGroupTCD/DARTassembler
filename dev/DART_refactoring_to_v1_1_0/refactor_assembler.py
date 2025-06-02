import itertools
from copy import deepcopy
from pathlib import Path
from typing import Union
import pandas as pd
from DARTassembler.src.assembly.ligand_geometries import align_donor_atoms, try_all_geometrical_isomer_possibilities
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from ast import literal_eval
import ase
import numpy as np
from DARTassembler.src.ligand_extraction.io_custom import save_to_xyz
import warnings

from dev.Assembler_revision_jan_2025.assembler.utilities import AssembledIsomer
from dev.test.Integration_Test import IntegrationTest

warnings.filterwarnings("ignore", category=UserWarning)
# warnings.simplefilter('error')    # Make warnings raise exceptions


axes = {
    'x': np.array([1, 0, 0]),
    '-x': np.array([-1, 0, 0]),
    'y': np.array([0, 1, 0]),
    '-y': np.array([0, -1, 0]),
    'z': np.array([0, 0, 1]),
    '-z': np.array([0, 0, -1]),
}

def angle_between_vectors(v1, v2, degrees=True) -> float:
    """
    Calculate the angle between two vectors.
    :param v1: Vector 1
    :param v2: Vector 2
    :param degrees: If True, return the angle in degrees. If False, return the angle in radians.
    :return: Angle between the two vectors.
    """
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    if degrees:
        angle = np.degrees(angle)
    return angle




# class AssemblyLigand(object):
#
#     def __init__(
#                     self,
#                     atoms: ase.Atoms,
#                     donor_idc: list[int]
#                     ):
#         self.atoms = atoms
#         self.donor_idc = donor_idc
#         self.denticity = len(donor_idc)
#         self.isomers = None     # Fill with list of ase.Atoms() objects. This is how the assembly should ideally be set up.
#         self.isomer_rssd = None
#
#     def _get_donor_positions(self, donor_idc=None) -> np.ndarray:
#         """
#         Get the positions of the donor atoms.
#         :param donor_idc: Indices of the donor atoms. If None, use all donor atoms.
#         :return: Array of donor atom positions.
#         """
#         if donor_idc is None:
#             donor_idc = self.donor_idc
#         donor_idc = list(donor_idc)
#         return self.atoms.positions[donor_idc]
#
#     def _get_average_donor_vector(self, donor_idc=None) -> np.ndarray:
#         """
#         Get the average vector between the donor atoms. For monodentate ligands, this is the vector from the metal to the donor atom. For bidentate ligands, this is the vector pointing right in the middle between the two donor atoms. For mer-tridentate ligands, the vector points roughly into the direction of the middle atom, but weighted slightly towards the positions of the two outer atoms.
#         :param donor_idc: Indices of the donor atoms. If None, use all donor atoms.
#         :return: Average donor vector as 3D numpy array.
#         """
#         return self._get_donor_positions(donor_idc=donor_idc).mean(axis=0)
#
#     def _get_axis_as_array(self, axis: str) -> np.ndarray:
#         """
#         Get the axis vector as a numpy array from the string representation.
#         :param axis_str: Axis string. Choose from 'x', '-x', 'y', '-y', 'z', '-z'.
#         :return: Axis vector as 3D numpy array.
#         """
#         if isinstance(axis, str):
#             axis = axes[axis]
#         return axis
#
#     def _check_denticity(self, expected_denticity: int):
#         """
#         Check if the denticity of the ligand matches the expected denticity.
#         :param expected_denticity: Expected denticity of the ligand.
#         """
#         if self.denticity != expected_denticity:
#             raise ValueError(f'This function is only for {expected_denticity}-dentate ligands, but got denticity {  self.denticity}.')
#
#     @classmethod
#     def from_ligand_return_geometrical_isomers(cls,
#                                                atoms: ase.Atoms,
#                                                donor_idc: list[int],
#                                                target_vectors: np.array
#                                                ):
#         best_rssd, all_best_vectors, all_best_isomers = try_all_geometrical_isomer_possibilities(atoms=atoms, donor_idc=donor_idc, target_vectors=target_vectors)
#
#         all_isomers_as_ligand_class = [cls(atoms=isomer, donor_idc=donor_idc) for isomer in all_best_isomers]
#
#         return all_isomers_as_ligand_class
#
#     # def rotate_monodentate(self, axis='z') -> ase.Atoms:
#     #     """
#     #     Rotate the monodentate ligand so that its donor atom ends up at the specified target axis. Importantly, the ligand atoms are always rotated with the origin as the center of the rotation so that the placement of the ligand towards the metal center is the same as in the CSD by design.
#     #     :param axis: Target axis for the placement of the donor atom.
#     #     :return: ASE Atoms object of the rotated ligand placed correctly at the target axis.
#     #     """
#     #     self._check_denticity(1)
#     #
#     #     axis = self._get_axis_as_array(axis)
#     #     # Simply rotate the monodentate ligand into position. By using the position of the ligand as extracted from the CSD (in which the metal was centered at the origin) this makes use of the information how the ligand is oriented towards the metal center in the CSD. Therefore, the rotation here is with the center at the origin. This will place e.g. the H of the OH ligand in the correct angle automatically.
#     #     self.atoms.rotate(self._get_average_donor_vector(), axis)   # center=(0, 0, 0) is the default in ASE
#     #
#     #     return self.atoms
#     #
#     # def rotate_bidentate(self, axis1='-x', axis2='-z') -> ase.Atoms:
#     #     """
#     #     Rotate the bidentate ligand so that its two donor atoms end up at the specified target axes.
#     #     :param axis1: Target axis for the placement of the first donor atom.
#     #     :param axis2: Target axis for the placement of the second donor atom.
#     #     :return: ASE Atoms object of the rotated ligand placed correctly at the target axes.
#     #     """
#     #     self._check_denticity(2)
#     #
#     #     axis1, axis2 = self._get_axis_as_array(axis1), self._get_axis_as_array(axis2)
#     #
#     #     self.atoms = align_donor_atoms(self.atoms, self.donor_idc, [axis1, axis2])
#     #
#     #     return self.atoms
#     #
#     # def _get_mer_tridentate_sorted_donor_idc(self) -> list[int]:
#     #     """
#     #     Return the donor indices of a mer-tridentate ligand so that the first and last index are the two outer atoms and the middle index is the middle atom.
#     #     :return: List of donor indices sorted in the correct order.
#     #     """
#     #     self._check_denticity(3)
#     #
#     #     donor_positions = self._get_donor_positions()
#     #     distances = np.linalg.norm(donor_positions[:, None] - donor_positions, axis=-1)
#     #     np.fill_diagonal(distances, -np.inf)
#     #     largest_distance_idc = list(np.unravel_index(np.argmax(distances), distances.shape))
#     #     pseudo_bidentate_idc = [self.donor_idc[i] for i in largest_distance_idc]
#     #     middle_idx = [i for i in self.donor_idc if i not in pseudo_bidentate_idc][0]
#     #     sorted_donor_idc = [pseudo_bidentate_idc[0], middle_idx, pseudo_bidentate_idc[1]]
#     #
#     #     return sorted_donor_idc
#     #
#     # def rotate_mer_tridentate(self, axis1='y', axis2='x', axis3='-y') -> ase.Atoms:
#     #     """
#     #     Rotate the tridentate ligand so that its three donor atoms end up at the specified target axes.
#     #     :param axis1: Target axis for the placement of the first donor atom (one of the outer atoms).
#     #     :param axis2: Target axis for the placement of the second donor atom, the middle atom.
#     #     :param axis3: Target axis for the placement of the third donor atom (the other outer atom).
#     #     :return: ASE Atoms object of the rotated ligand placed correctly at the target axes.
#     #     """
#     #     self._check_denticity(3)
#     #
#     #     axis1, axis2, axis3 = self._get_axis_as_array(axis1), self._get_axis_as_array(axis2), self._get_axis_as_array(axis3)
#     #
#     #     # Make sure these two arrays have the same order or atoms - first one of the outer atoms, then the middle atom, then the other outer atom
#     #     donor_idc = self._get_mer_tridentate_sorted_donor_idc()
#     #     target_axes = [axis1, axis2, axis3]
#     #
#     #     self.atoms = align_donor_atoms(self.atoms, donor_idc, target_axes)
#     #
#     #     return self.atoms
#
#
# class AssembledIsomer(object):
#
#     def __init__(
#                     self,
#                     atoms: ase.Atoms,
#                     ligands: list[AssemblyLigand]
#                     ):
#         """
#         Initialize the complex.
#         """
#         self.atoms = atoms
#         self.ligands = ligands
#
#     @classmethod
#     def from_ligand_combination_return_isomers(cls,
#                                                ligands: list[AssemblyLigand],
#                                                geometry: str,
#                                                metal: str
#                                                ):
#         ligand_isomers = []
#         for ligand, target_vectors in zip(ligands, geometry):
#             # Ensure target vectors is 2D
#             best_isomers = AssemblyLigand.from_ligand_return_geometrical_isomers(
#                 atoms=ligand.atoms,
#                 donor_idc=ligand.donor_idc,
#                 target_vectors=target_vectors
#             )
#             ligand_isomers.append(best_isomers)
#
#         # Get all complex isomers
#         isomers = []
#         ligand_isomer_combinations = list(itertools.combinations(ligand_isomers, len(ligand_isomers)))
#         for ligands in ligand_isomer_combinations:
#             isomer = cls.from_ligands(ligands=ligands, metal=metal)
#             isomers.append(isomer)
#
#         # Do simple 'isomer search' with ff-sp to find the best orientation of the monodentate ligand and maybe others.
#         # ...
#
#         return isomers
#
#     @classmethod
#     def from_ligands(cls,
#                       ligands: list[AssemblyLigand],
#                       metal: str
#                       ):
#         complex = ase.Atoms()
#         complex.append(ase.Atom(symbol=metal, position=[0, 0, 0]))
#         for ligand in ligands:
#             complex.extend(ligand)
#
#         return cls(atoms=complex, ligands=ligands)


    # def rotate_ligands_to_mer_3_2_1(self) -> None:
    #     """
    #     Returns the ligand atoms for the mer-3-2-1 geometry assembled around the origin, where the metal atom will be placed later.
    #     :return: ASE Atoms object of the ligands (without metal) in the mer-3-2-1 geometry.
    #     """
    #     for ligand in self.ligands:
    #         if ligand.denticity == 1:
    #             ligand.rotate_monodentate(axis='z')
    #         elif ligand.denticity == 2:
    #             ligand.rotate_bidentate(axis1='-x', axis2='-z')
    #         elif ligand.denticity == 3:
    #             ligand.rotate_mer_tridentate(axis1='y', axis2='x', axis3='-y')
    #         else:
    #             raise ValueError(f'Unexpected denticity {ligand.denticity}.')
    #
    #     return


if __name__ == '__main__':

    oer_ligand_db_path = '/Users/timosommer/PhD/projects/OERdatabase/data/testbatch/ligand_db/oer_all_ligands.jsonlines'
    OH_ligand_db_path = '/Users/timosommer/PhD/projects/OERdatabase/data/testbatch/ligand_db/oer_OH.jsonlines'
    oer_complexes_csv = '/Users/timosommer/PhD/projects/OERdatabase/data/testbatch/oerdb/oerdb_testbatch_v0.5/Ru_updated/info_table.csv'
    save_concat_xyz = 'data/assembler/test_new_assembler/data_output/concat_complexes.xyz'
    n_max = 10
    n_max_ligands = 300
    global_concat_atoms = []


    # Load ligand databases and append OH ligands to OER ligands to only have one database to search in
    db = LigandDB.load_from_json(path=oer_ligand_db_path, n_max=n_max_ligands)
    OH_db = LigandDB.load_from_json(path=OH_ligand_db_path)
    db.db.update(OH_db.db)

    # Read in OER complex csv and extract the ligand names of the successfully assembled complexes.
    df = pd.read_csv(oer_complexes_csv)
    df = df[df['success']]
    df['ligand names'] = df['ligand names'].apply(literal_eval)
    df = df.set_index('complex name')
    complexes = df.to_dict(orient='index')

    idx = 0
    for complex_name, row in complexes.items():
        # Extract the ligands from the database as old Ligand objects
        tridentate_name, bidentate_name, monodentate_name = row['ligand names']
        try:
            tridentate, bidentate, monodentate = db.db[tridentate_name], db.db[bidentate_name], db.db[monodentate_name]
        except KeyError:
            continue
        # Convert the old Ligand objects to AssemblyLigand objects which are used for the assembly
        # monodentate = AssemblyLigand(atoms=monodentate.mol, donor_idc=monodentate.ligand_to_metal)
        # bidentate = AssemblyLigand(atoms=bidentate.mol, donor_idc=bidentate.ligand_to_metal)
        # tridentate = AssemblyLigand(atoms=tridentate.mol, donor_idc=tridentate.ligand_to_metal)

        # Target vectors for mer-3-2-1 geometry.
        target_vectors = [
            [0, 0, 1],
            [[-1, 0, 0],[0, 0, -1]],
            [[0, 1, 0], [1, 0, 0], [0, -1, 0]],
             ]
        ligand_origins = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        metal_centers = [
            [ase.Atom(symbol='Ru', position=[0, 0, 0])],
            [ase.Atom(symbol='Fe', position=[1, 0, 0]), ase.Atom(symbol='Ru', position=[0, 0, 0])],
            [ase.Atom(symbol='Fe', position=[1, 0, 0])]
        ]
        ligands = [ligand for ligand in [monodentate, bidentate, tridentate]]
        isomers = AssembledIsomer.from_ligands_and_metal_centers(
                                                                    ligands=ligands,
                                                                    target_vectors=target_vectors,
                                                                    ligand_origins=ligand_origins,
                                                                    metal_centers=metal_centers,
                                                                    )

        for c_idx, isomer in enumerate(isomers):
            name = f'{complex_name}_{c_idx}'
            global_concat_atoms.append([isomer.mol, complex_name])
            # ase.visualize.view(atoms) # Uncomment to visualize each assembled complex in ASE

            if idx >= n_max:
                break
            idx += 1

    # Save all assembled complexes as concatenated xyz file so that they can be visualized in VESTA or Mercury
    structures = [atoms for atoms, _ in global_concat_atoms]
    comments = [name for _, name in global_concat_atoms]
    save_to_xyz(outpath=save_concat_xyz, structures=structures, comments=comments)
    print(f'Saved {len(global_concat_atoms)} complexes to {save_concat_xyz}')
    print('Done!')


    #%% ==============    Doublecheck refactoring    ==================
    save_concat_xyz = Path(save_concat_xyz)
    old_dir = save_concat_xyz.parent.parent / Path('benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=save_concat_xyz.parent, old_dir=old_dir)
        test.compare_all()
        print('Test for installation passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')





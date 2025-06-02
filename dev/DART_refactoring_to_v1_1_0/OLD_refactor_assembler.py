import pandas as pd
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from ast import literal_eval
import ase
import numpy as np
from DARTassembler.src.ligand_extraction.io_custom import save_to_xyz

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

class AssemblyLigand(object):

    def __init__(
                    self,
                    atoms: ase.Atoms,
                    donor_idc: list[int]
                    ):
        self.atoms = atoms
        self.donor_idc = donor_idc
        self.denticity = len(donor_idc)

    def _get_donor_positions(self, donor_idc=None) -> np.ndarray:
        """
        Get the positions of the donor atoms.
        :param donor_idc: Indices of the donor atoms. If None, use all donor atoms.
        :return: Array of donor atom positions.
        """
        if donor_idc is None:
            donor_idc = self.donor_idc
        donor_idc = list(donor_idc)
        return self.atoms.positions[donor_idc]

    def _get_average_donor_vector(self, donor_idc=None) -> np.ndarray:
        """
        Get the average vector between the donor atoms. For monodentate ligands, this is the vector from the metal to the donor atom. For bidentate ligands, this is the vector pointing right in the middle between the two donor atoms. For mer-tridentate ligands, the vector points roughly into the direction of the middle atom, but weighted slightly towards the positions of the two outer atoms.
        :param donor_idc: Indices of the donor atoms. If None, use all donor atoms.
        :return: Average donor vector as 3D numpy array.
        """
        return self._get_donor_positions(donor_idc=donor_idc).mean(axis=0)

    def _get_axis_as_arrray(self, axis: str) -> np.ndarray:
        """
        Get the axis vector as a numpy array from the string representation.
        :param axis_str: Axis string. Choose from 'x', '-x', 'y', '-y', 'z', '-z'.
        :return: Axis vector as 3D numpy array.
        """
        if isinstance(axis, str):
            axis = axes[axis]
        return axis

    def rotate_monodentate(self, axis='z') -> ase.Atoms:
        """
        Rotate the monodentate ligand so that its donor atom ends up at the specified target axis.
        :param axis: Target axis for the placement of the donor atom.
        :return: ASE Atoms object of the rotated ligand.
        """
        if self.denticity != 1:
            raise ValueError(f'This function is only for monodentate ligands, but got denticity {self.denticity}.')

        axis = self._get_axis_as_arrray(axis)
        self.atoms.rotate(self._get_average_donor_vector(), axis)

        return self.atoms

    def rotate_bidentate(self, axis1='-x', axis2='-z') -> ase.Atoms:
        """
        Rotate the bidentate ligand so that its two donor atoms end up at the specified target axes.
        :param axis1: Target axis for the placement of the first donor atom.
        :param axis2: Target axis for the placement of the second donor atom.
        :return: ASE Atoms object of the rotated ligand.
        """
        if self.denticity != 2:
            raise ValueError('This function is only for bidentate ligands.')

        axis1, axis2 = self._get_axis_as_arrray(axis1), self._get_axis_as_arrray(axis2)
        # First rotate the ligand so that the average donor vector is between -x and -z.
        target = axis1 + axis2
        self.atoms.rotate(self._get_average_donor_vector(), target)

        # Now rotate the ligand around the right angle so that the two donor atoms are in the xz plane.
        angle1 = angle_between_vectors(self.atoms.positions[self.donor_idc[0]], [0, 0, 1])
        angle2 = angle_between_vectors(self.atoms.positions[self.donor_idc[1]], [0, 0, 1])
        avg_angle = (angle1 + angle2) / 2
        self.atoms.rotate(avg_angle, self._get_average_donor_vector())

        return self.atoms

    def rotate_mer_tridentate(self, axis1='y', axis2='x', axis3='-y') -> ase.Atoms:
        """
        Rotate the tridentate ligand so that its three donor atoms end up at the specified target axes.
        :param axis1: Target axis for the placement of the first donor atom.
        :param axis2: Target axis for the placement of the second donor atom.
        :param axis3: Target axis for the placement of the third donor atom.
        :return: ASE Atoms object of the rotated ligand.
        """
        if self.denticity != 3:
            raise ValueError('This function is only for tridentate ligands.')

        axis1, axis2, axis3 = self._get_axis_as_arrray(axis1), self._get_axis_as_arrray(axis2), self._get_axis_as_arrray(axis3)

        # Get the pseudo bidentate ligand indices, which are the indices of the two donor atoms which are the farthest apart.
        donor_positions = self._get_donor_positions()
        distances = np.linalg.norm(donor_positions[:, None] - donor_positions, axis=-1)
        np.fill_diagonal(distances, -np.inf)
        pseudo_bidentate_idc = np.unravel_index(np.argmax(distances), distances.shape)
        middle_idx = [i for i in range(3) if i not in pseudo_bidentate_idc][0]

        # First rotate the ligand so that the average donor vector hits x.
        target = [1, 0, 0]
        # self.atoms.rotate(self.get_average_donor_vector(donor_idc=[middle_idx]), target)

        # Now rotate the ligand around the right angle so that the two donor atoms are in the xz plane.
        # idx1, idx2 = pseudo_bidentate_idc
        # angle1 = angle_between_vectors(self.mol.positions[self.donor_idc[idx1]], [0, -1, 0])
        # angle2 = angle_between_vectors(self.mol.positions[self.donor_idc[idx2]], [0, 1, 0])
        # avg_angle = (angle1 + angle2) / 2
        # self.mol.rotate(avg_angle, self.get_average_donor_vector(donor_idc=pseudo_bidentate_idc))

        # Todo: Get the angle between the donors at the beggining to directly rotate the ligand to the right position knowing the average ligand deviation from the target.

        return self.atoms


class AssemblyComplex(object):

    def __init__(
                    self,
                    ligands: list[AssemblyLigand],
                    geometry: str,
                    metal: str
                    ):
        """
        Initialize the complex with the ligands, the geometry and the metal.
        :param ligands: List of AssemblyLigand objects.
        :param geometry: Geometry of the complex. Choose from 'mer-3-2-1'.
        :param metal: Chemical symbol of the metal atom.
        """
        self.ligands = ligands
        self.geometry = geometry
        self.metal = metal

    def assemble_complex(self) -> ase.Atoms:
        """
        Assemble the complex from the ligands and the geometry.
        :return: ASE Atoms object of the complex with the metal at the origin and the ligands placed around it according to the geometry.
        """
        # Initialize the complex with the metal atom at the origin
        complex = ase.Atoms()
        complex.append(ase.Atom(symbol=self.metal, position=[0, 0, 0]))

        # Get all the ligand atoms oriented correctly around the origin.
        if self.geometry == 'mer-3-2-1':
            ligand_atoms = self.get_mer_3_2_1_ligands()
        else:
            raise ValueError(f'Unexpected geometry "{self.geometry}". Choose from "mer-3-2-1".')

        # Append the ligand atoms to the complex
        complex.extend(ligand_atoms)

        return complex

    def get_mer_3_2_1_ligands(self) -> ase.Atoms:
        """
        Returns the ligand atoms for the mer-3-2-1 geometry assembled around the origin, where the metal atom will be placed later.
        :return: ASE Atoms object of the ligands (without metal) in the mer-3-2-1 geometry.
        """
        ligand_atoms = ase.Atoms()

        for ligand in self.ligands:
            if ligand.denticity == 1:
                mono_atoms = ligand.rotate_monodentate()
                ligand_atoms.extend(mono_atoms)
            elif ligand.denticity == 2:
                bi_atoms = ligand.rotate_bidentate()
                ligand_atoms.extend(bi_atoms)
            elif ligand.denticity == 3:
                tri_atoms = ligand.rotate_mer_tridentate()
                ligand_atoms.extend(tri_atoms)
            else:
                raise ValueError(f'Unexpected denticity {ligand.denticity}.')

        return ligand_atoms




if __name__ == '__main__':

    oer_ligand_db_path = '/Users/timosommer/PhD/projects/OERdatabase/data/testbatch/ligand_db/oer_all_ligands.jsonlines'
    OH_ligand_db_path = '/Users/timosommer/PhD/projects/OERdatabase/data/testbatch/ligand_db/oer_OH.jsonlines'
    oer_complexes_csv = '/Users/timosommer/PhD/projects/OERdatabase/data/testbatch/oerdb/oerdb_testbatch_v0.5/Ru_updated/info_table.csv'
    save_concat_xyz = 'concat_complexes.xyz'
    n_max = 1

    # Load ligand databases and append OH ligands to OER ligands to only have one database to search in
    db = LigandDB.load_from_json(path=oer_ligand_db_path)
    OH_db = LigandDB.load_from_json(path=OH_ligand_db_path)
    db.db.update(OH_db.db)

    # Read in OER complex csv and extract the ligand names of the successfully assembled complexes.
    df = pd.read_csv(oer_complexes_csv)
    df = df[df['success']]
    df['ligand names'] = df['ligand names'].apply(literal_eval)
    df = df.set_index('complex name')
    complexes = df.to_dict(orient='index')

    comments = []
    structures = []
    for complex_name, row in complexes.items():
        # Extract the ligands from the database as old Ligand objects
        tridentate_name, bidentate_name, monodentate_name = row['ligand names']
        tridentate, bidentate, monodentate = db.db[tridentate_name], db.db[bidentate_name], db.db[monodentate_name]
        # Convert the old Ligand objects to AssemblyLigand objects which are used for the assembly
        monodentate = AssemblyLigand(atoms=monodentate.mol, donor_idc=monodentate.ligand_to_metal)
        bidentate = AssemblyLigand(atoms=bidentate.mol, donor_idc=bidentate.ligand_to_metal)
        tridentate = AssemblyLigand(atoms=tridentate.mol, donor_idc=tridentate.ligand_to_metal)

        complex = AssemblyComplex(
                                    ligands=[monodentate, bidentate, tridentate],
                                    geometry='mer-3-2-1',
                                    metal='Ru'
                                    )
        atoms = complex.assemble_complex()
        structures.append(atoms)
        comments.append(complex_name)
        # ase.visualize.view(atoms) # Uncomment to visualize each assembled complex in ASE

        if len(structures) == n_max:
            break

    # Save all assembled complexes as concatenated xyz file so that they can be visualized in VESTA or Mercury
    save_to_xyz(outpath=save_concat_xyz, structures=structures, comments=comments)





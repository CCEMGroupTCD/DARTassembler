from utilities import angle_between_vectors, axes, align_donor_atoms
import numpy as np
import ase
import numpy as np
import networkx as nx
from DARTassembler.src.ligand_extraction.DataBase import MoleculeDB
from scipy.optimize import differential_evolution


class AssemblyLigand(object):

    def __init__(
            self,
            atoms: ase.Atoms,
            donor_idc: list[int],
            graph: nx.Graph
    ):
        self.atoms = atoms
        self.donor_idc = donor_idc
        self.denticity = len(donor_idc)
        self.m_coords = m_coords = np.ndarray([0, 0, 0])  # implied by the CSD output
        self.graph = graph

    def check_axis_input(self, axis: str):
        """
        Check if the input axis is valid.
        :param axis: Input axis as a string.
        :return: None
        """
        if axis not in axes.keys():
            raise ValueError(f"Invalid axis {axis}. Choose from {list(axes.keys())}.")

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

    def _get_middle_vector(self, vector1: np.ndarray, vector2: np.ndarray) -> np.ndarray:
        """
        Get the vector that bisects the angle between two other vectors.

        :param vector1: First vector as a 3D numpy array.
        :param vector2: Second vector as a 3D numpy array.
        :return: Middle vector as a normalized 3D numpy array.
        """
        middle_vector = vector1 + vector2
        norm = np.linalg.norm(middle_vector)

        # Check for zero magnitude to avoid division by zero
        if norm == 0:
            raise ValueError("The input vectors are directly opposite, resulting in a zero vector.")

        return middle_vector / norm

    def _find_closeness_middle_node(self) -> int:
        pass
        # on second thought this is not a great idea because you could end up with a shorted graph distance between the two outer atoms than the inner atoms

    def rotate_to_minimize_max_x(self, angle_limit: float = 30.0, norm_axis: str = "z", target_axis: str = "x"):
        """
        Rotate an ASE object around the z-axis to minimize the x-coordinate of the atom
        with the largest x-coordinate.
        :param target_axis:  The axis that the ligand is to be aligned to
        :param norm_axis: the axis which the ligand is to be rotated around
        :param angle_limit: Maximum rotation angle in degrees. If the optimal rotation angle is larger than this, the ligand is not rotated.
        :return: Rotated ASE Atoms object and the optimal rotation angle (radians).
        """

        def compute_distances_to_plane(positions, plane_normal, plane_point):
            """
            Compute distances of points to a specified plane in 3D space.

            :param positions: numpy array of shape (n, 3), containing the x, y, z coordinates of the points.
            :param plane_normal: numpy array of shape (3,), the normal vector to the plane [a, b, c].
            :param plane_point: numpy array of shape (3,), a point on the plane (x0, y0, z0).
            :return: numpy array of distances of each point to the plane.
            """
            # Normalize the plane normal vector
            plane_normal = plane_normal / np.linalg.norm(plane_normal)

            # Vector from plane_point to positions
            vectors = positions - plane_point

            # Compute perpendicular distances
            distances = np.abs(np.dot(vectors, plane_normal))

            return distances

        def objective_function(theta):
            atoms_copy = self.atoms.copy()
            atoms_copy.rotate(a=theta[0], v=norm_axis)  # Rotate around z-axis

            # Find the atoms with negative x-coordinates
            negative_indices = np.where(atoms_copy.positions[:, 0] < 0)[0]
            negative_positions = atoms_copy.positions[negative_indices]

            # Compute distances to the plane. The normal of this plane will be the cross product of the norm_axis and target_axis
            plane_norm = np.cross(self._get_axis_as_arrray(norm_axis), self._get_axis_as_arrray(target_axis))
            distances = compute_distances_to_plane(negative_positions, plane_norm, [0, 0, 0])

            # Find the index of the atom with the smallest distance
            return -1 * distances[np.argmin(distances)]  # closest distance to the plane

        # Optimize the rotation angle to minimize the largest x-coordinate
        try:
            result = differential_evolution(func=objective_function, bounds=[(0, 360)])
            optimal_angle = result.x[0]
        except (RuntimeError, ValueError):
            optimal_angle = angle_limit + 1  # If the optimization fails, we don't rotate the ligand

        # Apply the optimal rotation to the original atoms
        if abs(optimal_angle) < angle_limit:
            self.atoms.rotate(optimal_angle, self._get_axis_as_arrray('z'))
            print(f"Angle adjusted [{optimal_angle:.2f}]")
        else:
            pass  # The proposed rotation is too large, so we don't rotate the ligand

    def get_sorted_donor_idc_for_mer_tridentate(self, atoms, donor_idc) -> list[int]:
        """
        Return the indices of the donor atoms of a mer-tridentate ligand so that indices are in the following order: outer atom, central atom, other outer atom
        :param atoms: ASE atoms object of the ligand
        :param donor_idc: List of donor atom indices
        :return: List of donor atom indices in the correct order
        """
        # Define the target vectors so that the central donor atom is the second one
        target_vectors = [[0, 1, 0], [1, 0, 0], [0, -1, 0]]
        all_rssds = []
        for central_idx in donor_idc:
            outer_idc = [i for i in donor_idc if i != central_idx]
            try_idc_order = [outer_idc[0], central_idx, outer_idc[1]]

            # Try to align the donor atoms and get the resulting rssd
            _, rssd = align_donor_atoms(atoms, try_idc_order, target_vectors, return_rssd=True)
            all_rssds.append(rssd)

        # Find the donor indices with the smallest rssd
        central_idx = donor_idc[np.argmin(all_rssds)]
        # Order the donor atoms so that the central atom is the second one
        other_idc = [i for i in donor_idc if i != central_idx]
        donor_idc_ordered = [other_idc[0], central_idx, other_idc[1]]

        return donor_idc_ordered

    def rotate_monodentate(self, axis='z') -> ase.Atoms:
        """
        Rotate the monodentate ligand so that its donor atom ends up at the specified target axis. Importantly, the ligand atoms are always rotated with the origin as the center of the rotation so that the placement of the ligand towards the metal center is the same as in the CSD by design.
        :param axis: Target axis for the placement of the donor atom.
        :return: ASE Atoms object of the rotated ligand placed correctly at the target axis.
        """
        if self.denticity != 1:
            raise ValueError(f'This function is only for monodentate ligands, but got denticity {self.denticity}.')

        axis = self._get_axis_as_arrray(axis)
        # Simply rotate the monodentate ligand into position. By using the position of the ligand as extracted from the CSD (in which the metal was centered at the origin) this makes use of the information how the ligand is oriented towards the metal center in the CSD. Therefore, the rotation here is with the center at the origin. This will place e.g. the H of the OH ligand in the correct angle automatically.
        self.atoms.rotate(self._get_average_donor_vector(), axis)  # center=(0, 0, 0) is the default in ASE

        return self.atoms

    def rotate_bidentate(self, axis1='-x', axis2='-z') -> ase.Atoms:
        """
        Rotate the bidentate ligand so that its two donor atoms end up at the specified target axes.
        :param axis1: Target axis for the placement of the first donor atom.
        :param axis2: Target axis for the placement of the second donor atom.
        :return: ASE Atoms object of the rotated ligand placed correctly at the target axes.
        """
        if self.denticity != 2:
            raise ValueError('This function is only for bidentate ligands.')

        # Step 1: Calculate the normal vector to the donor-atom / metal plane
        d1_coords, d2_coords = self._get_donor_positions()
        normal_vector = np.cross(d1_coords, d2_coords)  # If metal was not at the origin, we would need to subtract the metal position from the donor positions
        normal_vector /= np.linalg.norm(normal_vector)  # Normalize the normal vector to the plane of donor atoms and metal

        # Step 2: Define the y-axis as the target normal
        target_normal = np.array([0, 1, 0])

        # Step 3: Rotate the ligand so that the M-d1-d2 plane is aligned with (essentially inside) the xz-plane
        self.atoms.rotate(a=normal_vector, v=target_normal)  # v is always target and a is always the vector to rotate around

        # Step 4: Rotate the ligand so that the bisector of the d1-d2 vector is aligned 45 degrees between the x and z axes
        d1_coords, d2_coords = self._get_donor_positions()
        middle_vector = self._get_middle_vector(d1_coords, d2_coords)  # This vector is the bisector of the d1-d2 vector
        target_vector = np.array([-1, 0, -1])  # This vector is 45 degrees between the x and z axes
        self.atoms.rotate(a=middle_vector, v=target_vector)

        return self.atoms

    def rotate_mer_tridentate(self, norm_axis: str = "z", target_axis: str = "x") -> ase.Atoms:
        """
        Rotate the tridentate ligand so that its three donor atoms end up at the specified target axes.
        :param norm_axis:   The axis that the normal vector to the donor atoms is to be aligned to
        :param target_axis: The axis the 'bulk' of the molecule is to be aligned to
        :return: ASE Atoms object of the rotated ligand placed correctly at the target axes.
        """
        # Checks ---> orthogonality of the axes and the axis input
        assert norm_axis != target_axis, f"norm: [{norm_axis}] target: [{target_axis}] The normal axis and the target axis must be orthogonal."
        self.check_axis_input(norm_axis)
        self.check_axis_input(target_axis)

        # Step 1: Calculate the normal vector to the donor-atoms plane and orientate the donor atoms such that they are in a cartesian plane
        d1_coords, d2_coords, d3_coords = self._get_donor_positions()
        v1 = d2_coords - d1_coords
        v2 = d3_coords - d1_coords
        normal_vector = np.cross(v1, v2)
        self.atoms.rotate(a=normal_vector, v=norm_axis, center=(0, 0, 0))

        # Step 2: Identify the middle donor atom and its coordinates
        mid_donor_vector = self.atoms.positions[self.get_sorted_donor_idc_for_mer_tridentate(self.atoms.copy(), self.donor_idc)[1]]

        # Step 3: Adjust the middle donor vector such that it is orthogonal to the norm_axis (I noticed issues when this was NOT the case)
        if norm_axis in ['x', '-x']:
            adjusted_vector = (0, mid_donor_vector[1], mid_donor_vector[2])
        elif norm_axis in ['y', '-y']:
            adjusted_vector = (mid_donor_vector[0], 0, mid_donor_vector[2])
        elif norm_axis in ['z', '-z']:
            adjusted_vector = (mid_donor_vector[0], mid_donor_vector[1], 0)
        else:
            raise ValueError(f'Unexpected axis {norm_axis}.')

        # Step 4: calculate the angle between the middle donor vector and the z-axis
        self.atoms.rotate(a=adjusted_vector, v=target_axis)

        negative_indices = np.where(self.atoms.positions[:, 0] < 0)[0]
        if len(negative_indices) > 0:
            self.rotate_to_minimize_max_x(angle_limit=30.0, norm_axis="z", target_axis="x")

        return self.atoms


class AssemblyComplex(object):

    def __init__(
            self,
            ligands: list[MoleculeDB],
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
        self.assembly_ligands = [AssemblyLigand(atoms=ligand.mol, donor_idc=ligand.ligand_to_metal, graph=ligand.get_reindexed_graph()) for ligand in ligands]
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

        # Do simple 'isomer search' to find the best orientation of the monodentate ligand and maybe others. We need to think about in which step this would be the smartest to do.
        # ...

        return complex

    def get_mer_3_2_1_ligands(self) -> ase.Atoms:
        """
        Returns the ligand atoms for the mer-3-2-1 geometry assembled around the origin, where the metal atom will be placed later.
        :return: ASE Atoms object of the ligands (without metal) in the mer-3-2-1 geometry.
        """
        ligand_atoms = ase.Atoms()

        for ligand in self.assembly_ligands:
            if ligand.denticity == 1:
                mono_atoms = ligand.rotate_monodentate()
                ligand_atoms.extend(mono_atoms)
            elif ligand.denticity == 2:
                bi_atoms = ligand.rotate_bidentate()
                ligand_atoms.extend(bi_atoms)  # todo: outcomment this line when working on the bidentate to show the bidentate in the xyz
            elif ligand.denticity == 3:
                tri_atoms = ligand.rotate_mer_tridentate()  # todo: outcomment this line when working on the tridentate to show the tridentate in the xyz
                ligand_atoms.extend(tri_atoms)
            else:
                raise ValueError(f'Unexpected denticity {ligand.denticity}.')

        return ligand_atoms

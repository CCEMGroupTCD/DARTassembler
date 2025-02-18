import ase.io.trajectory
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.lj import LennardJones
from ase.calculators.morse import MorsePotential
from ase.visualize import view
from typing import List
from scipy.spatial.transform import Rotation
from ase.io.trajectory import Trajectory
from ase.calculators.calculator import Calculator, all_changes
import numpy as np

print("Do not implement the following code as it didn not really work as expected and the code is not complete.")

class MyCalculator(Calculator):
    """
    Custom ASE Calculator for a harmonic potential.

    This calculator computes the energy and forces for an atomic system
    using a harmonic potential. It is useful as a starting point for
    developing more complex calculators.
    """
    # Define which properties this calculator can compute.
    implemented_properties = ['energy', 'forces']

    def __init__(self, vectors=None, origins=None, tags_list=None):
        """
        Initialize the custom calculator.

        Parameters:
        vectors : List[np.ndarray]
            List of rotation axes for different fragments.
        origins : List[np.ndarray]
            List of origins for rotation.
        tags_list : List[int]
            List of atom tags identifying different fragments.
        """
        super().__init__()
        self.vectors = vectors  # Store rotation axes
        self.origins = origins  # Store rotation origins
        self.tags_list = tags_list  # Store fragment tags

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """
        Calculate the energy and forces for the provided atomic configuration.
        :param: atoms: The atomic configuration (positions, numbers, etc.) for which the energy and forces will be calculated.
        :param: properties: List of properties to calculate (default is ['energy']).
        :param: system_changes: List of changes that have occurred in the atomic configuration. Defaults to all_changes, meaning that any change will trigger a full recalculation.
        :param: tags_list: List of tags to distinguish ligands in the complex
        :return: None
        """
        # Call the base class calculate() to handle bookkeeping and system_changes.
        if properties is None:
            properties = ['energy']
        super().calculate(atoms, properties, system_changes)

        # ---Step 1--- Get positions (N, 3), distances (N, N) between atoms and the displacement vectors between atoms (N, N, 3)
        positions = atoms.get_positions()
        distance_matrix = atoms.get_all_distances()
        np.fill_diagonal(distance_matrix, 1.0)  # Avoid division by zero
        displacement_vector_matrix = positions[:, None, :] - positions[None, :, :]

        # ---Step 2--- Compute the forces on each atom F = k * (R_ij / D_ij**n)
        # F = Force; R_ij = vector between atom i and j; D_ij = distance between atom i and j; n = power
        # As D_ij increases (atoms become further apart), the force decreases
        k = 0.1  # Spring constant
        n = 3    # Power (n = 3 means we have inverse square law (1/D_ij**2) and the force is normalized 1/D_ij yielding 1/D_ij**3)
        forces = np.sum(k * displacement_vector_matrix / distance_matrix[..., None] ** n, axis=1)   # Yields the force on each atom

        # ---Step 3--- Compute the energy of the system
        energy = -1.0 * np.sum(k * distance_matrix) / 2     # we divide by 2 to avoid double counting

        # ---Step 4--- use the calculated forces to generate torques
        torques = self._calc_torques(atom_positions=positions, atom_forces=forces, axes=self.vectors, axis_origins=self.origins, tags_list=self.tags_list)



        # ---Step 5--- Store the results in the 'results' dictionary
        self.results = {
            'energy': energy,
            'forces': torques
        }

    @staticmethod
    def _calc_torques(atom_positions, atom_forces, axes, axis_origins, tags_list):
        """
        Calculate the torques on each atom in the system using the calculated forces.
        :param: atom_positions: The positions of the atoms in the system.
        :param: atom_forces: The forces acting on each atom in the system.
        :param: axes: The axes of rotation for each fragment.
        :param: axis_origins: The origin of rotation for each fragment.
        :param: tags_list: The list of tags to distinguish ligands in the complex.
        :return: torques: The torques acting on each atom in the system.
        """
        torques = np.zeros_like(atom_positions)  # Initialize torques to zero

        # Get unique fragment tags, excluding tag 0 (metal center)
        unique_tags = [tag for tag in set(tags_list) if tag != 0]

        # Loop over each unique tag in the system and its corresponding axis of rotation and origin
        for tag, axis, origin in zip(set(unique_tags), axes, axis_origins):
            # Get the indices of the atoms with the current tag
            fragment_indices = [i for i, t in enumerate(tags_list) if t == tag]
            # Loop through all the atoms in the fragment
            for i in fragment_indices:
                # torque = (r_i - O) X F ; where r_i is the position of atom i, O is the origin of rotation, F is the force on atom i
                r = atom_positions[i] - origin      # position of atom i relative to the origin
                tau = np.cross(r, atom_forces[i])   # Compute torque
                tau_aligned = np.dot(tau, axis) * axis  # This is teh torque projected onto a specific axis
                torques[i] = tau_aligned  # Store rotational force

        return torques

# Define the Rotation Function (Outside the Calculator)
def apply_rotation(atoms, axes, origins, tags_list, torques, step_size=0.1):
    """ Rotates molecular fragments based on torques. """
    new_positions = atoms.get_positions().copy()

    unique_tags = [tag for tag in set(tags_list) if tag != 0]  # Exclude metal center

    for tag, axis, origin in zip(unique_tags, axes, origins):
        fragment_indices = [i for i, t in enumerate(tags_list) if t == tag]

        net_torque = np.sum(torques[fragment_indices], axis=0)
        net_torque_magnitude = np.linalg.norm(net_torque)

        if net_torque_magnitude > 0:
            I = np.sum(np.dot(new_positions[fragment_indices] - origin, axis) ** 2)

            if I > 0:
                theta = (net_torque_magnitude / I) * step_size  # Compute rotation angle
                rotation = Rotation.from_rotvec(theta * axis)

                for i in fragment_indices:
                    new_positions[i] = origin + rotation.apply(new_positions[i] - origin)

    atoms.set_positions(new_positions)


TMC = ase.io.read('TMC.xyz')

# set tags to distinguish ligands in the complex
TMC[0].tag = 0  # metal center
for i in range(1, 12):  # atoms 1-11 (inclusive) are ligands
    TMC[i].tag = 1  # pyridine ligand 1
for i in range(12, 23):  # atoms 12-22 (inclusive) are ligands
    TMC[i].tag = 2  # pyridine ligand 2
for i in range(23, 34):  # atoms 23-33 (inclusive) are ligands
    TMC[i].tag = 3  # pyridine ligand 3
for i in range(34, 45):  # atoms 34-44 (inclusive) are ligands
    TMC[i].tag = 4  # pyridine ligand 4
for i in range(45, 48):  # atoms 45-47 (inclusive) are ligands
    TMC[i].tag = 5  # H2O ligand 1
for i in range(48, 73):  # atoms 48-72 (inclusive) are ligands
    TMC[i].tag = 6  # big ligand 1

vectors = [np.array([0.0, 1.0, 0.0]),
           np.array([1.0, 0.0, 0.0]),
           np.array([0.0, -1.0, 0.0]),
           np.array([-1.0, 0.0, 0.0]),
           np.array([0.0, 0.0, 1.0]),
           np.array([0.0, 0.0, -1.0])]

origins = [np.array([0.0, 0.0, 0.0]),
           np.array([0.0, 0.0, 0.0]),
           np.array([0.0, 0.0, 0.0]),
           np.array([0.0, 0.0, 0.0]),
           np.array([0.0, 0.0, 0.0]),
           np.array([0.0, 0.0, 0.0])]

TMC.set_calculator(MyCalculator(vectors=vectors, origins=origins, tags_list=TMC.get_tags()))
dyn = BFGS(TMC, trajectory='TMC.traj')
dyn.run(fmax=0.05)

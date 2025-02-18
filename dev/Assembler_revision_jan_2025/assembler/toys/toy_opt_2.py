import ase.io
from ase.visualize import view
from typing import List
from ase import Atoms
import numpy as np
from scipy.optimize import differential_evolution
from scipy.spatial.transform import Rotation as R

def load_atoms():
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

    return TMC, vectors, origins


def rotate(atoms: Atoms, vector: np.array, origin: np.array, idc: List[int], angle: int):
    """
    Rotate the atoms in the Atoms object (only atoms with indices=idc) around the vector by the specified angle.
    :param atoms: Atoms object to rotate.
    :param vector: vector to rotate around.
    :param origin: origin of the rotation.
    :param idc: indices of the atoms to rotate.
    :param angle: the angle to rotate the atoms by in degrees.
    :return: an ase.Atoms object with the rotated atoms.
    """

    # Normalize rotation vector
    vector = np.asarray(vector, dtype=float)
    vector /= np.linalg.norm(vector)

    # Create rotation object
    rotation = R.from_rotvec(np.radians(angle) * vector)

    # Copy the atoms object to avoid modifying the original
    rotated_atoms = atoms.copy()

    # Apply rotation to selected atoms
    for i in idc:
        pos = atoms.positions[i] - origin  # Translate to origin
        rotated_pos = rotation.apply(pos) + origin  # Rotate and translate back
        rotated_atoms.positions[i] = rotated_pos

    return rotated_atoms

def objective(x: np.ndarray, vectors_in, origins_in, TMC_in: Atoms, tags_list_in):
    """
    Objective function to optimize the position of the ligands in the TMC complex.
    :param: tags_list_in: tags to distinguish how the ligands are rotated.
    :param: x: List of parameters to optimize.
    :param: vectors_in: List of vectors that ligands will be rotated around by x[n] degrees.
    :param: origins_in: List of origins that ligands will be rotated around.
    :param: TMC_in: The TMC complex to optimize.
    :return: The objective function value.
    """

    # ---Step 1--- Remove the metal tag from the tag list to avoid rotating the metal center
    unique_tags = [tag for tag in set(tags_list_in) if tag != 0]

    # ---Step 2--- Loop through the unique tags and rotate the ligands based on the parameters
    for angle, axis, origin, tag in zip(x.tolist(), vectors_in, origins_in, unique_tags):

        # Identify the indices of the ligands with the current tag
        fragment_indices = [i for i, t in enumerate(tags_list_in) if t == tag]
        TMC_in = rotate(atoms=TMC_in, vector=axis, origin=origin, idc=fragment_indices, angle=angle).copy()


    # ---Step 3--- Calculate the sum of the distances between atoms in the complex
    distance_matrix = TMC_in.get_all_distances()



    return -1.0*np.sum(distance_matrix)/2




if __name__ == "__main__":
    TMC, vectors, origins = load_atoms()
    view(TMC)

    # Specify bounds for each parameter: (lower, upper) for x and y respectively.
    bounds = [[0, 360] for _ in vectors]

    # Run the global optimizer.
    result = differential_evolution(objective, bounds=bounds, args=(vectors, origins, TMC, TMC.get_tags()))

    unique_tags = [tag for tag in set(TMC.get_tags().tolist()) if tag != 0]
    for angle, axis, origin, tag in zip(result.x.tolist(), vectors, origins, unique_tags):

        fragment_indices = [i for i, t in enumerate(TMC.get_tags().tolist()) if t == tag]
        TMC = rotate(atoms=TMC, vector=axis, origin=origin, idc=fragment_indices, angle=angle).copy()

    view(TMC)
    print(result)

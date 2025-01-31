from pathlib import Path
from typing import Union
import ase
import numpy as np
from scipy.spatial.transform import Rotation as R
from ase.io import write

axes = {
    'x': np.array([1, 0, 0]),
    '-x': np.array([-1, 0, 0]),
    'y': np.array([0, 1, 0]),
    '-y': np.array([0, -1, 0]),
    'z': np.array([0, 0, 1]),
    '-z': np.array([0, 0, -1]),
}

def save_to_xyz(outpath: Union[str,Path], structures: list[ase.Atoms], comments: list[str]):
    """
    Save a list of ASE Atoms objects to an XYZ file.
    :param outpath: Path to the output file.
    :param structures: List of ASE Atoms objects.
    :param comments: List of comments for each structure. Must be the same length as the structures list.
    :return: None
    """
    outpath = Path(outpath)
    if isinstance(comments, str):
        comments = [comments]
    if isinstance(structures, ase.Atoms):
        structures = [structures]
    if len(structures) != len(comments):
        raise ValueError(f"Number of structures ({len(structures)}) and comments ({len(comments)}) do not match.")

    outpath.unlink(missing_ok=True)  # Remove the file if it exists to avoid appending to an existing file.
    outpath.parent.mkdir(parents=True, exist_ok=True)
    for structure, comment in zip(structures, comments):
        with open(outpath, 'a') as f:
            write(f, structure, format='xyz', comment=comment)

    return

def align_donor_atoms(
                        atoms: ase.Atoms,
                        ligand_idc: list[int],
                        target_vectors: list[list[float]],
                        return_rssd: bool = False
                        ):
    """
    Align the donor atoms of a ligand to the target vectors.
    :param atoms: ASE Atoms object of all the atoms of the ligand.
    :param ligand_idc: Indices of the donor atoms of the ligand.
    :param target_vectors: A list of 3D target vectors to which the donor atoms should be aligned.
    :return: ASE Atoms object of the ligand with the donor atoms aligned to the target vectors.
    """
    target_vectors = np.array(target_vectors)
    assert len(target_vectors) == len(ligand_idc), 'The number of target vectors must match the number of donor atoms.'
    assert target_vectors.shape[1] == 3, 'The target vectors must be 3D vectors.'

    donor_vectors = atoms.positions[ligand_idc]
    # Normalize the donor vectors and target vectors to unit vectors to avoid scaling issues
    donor_vectors = donor_vectors / np.linalg.norm(donor_vectors, axis=1)[:, None]
    target_vectors = target_vectors / np.linalg.norm(target_vectors, axis=1)[:, None]
    # Find the correct rotation to align the donor vectors with the target vectors
    rot, rssd = R.align_vectors(a=target_vectors, b=donor_vectors)  # the a and b are unintuitive but correct
    # Apply the rotation to all the atoms of the ligand
    rotated_coords = rot.apply(atoms.positions)
    atoms.set_positions(rotated_coords)

    if return_rssd:
        return atoms, rssd
    else:
        return atoms

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
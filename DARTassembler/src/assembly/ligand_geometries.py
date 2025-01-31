import ase
import numpy as np
import itertools
from scipy.spatial.transform import Rotation as R
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

def xy_angle(angle: float) -> np.ndarray:
    """
    Get the coordinates of a vector in the x-y plane with a given angle from the y-axis.
    :param angle: Angle in degrees
    :return: Vector in the x-y plane
    """
    angle = np.radians(angle)   # Convert angle to radians
    return np.array([np.sin(angle), np.cos(angle), 0])

x = np.array([1, 0, 0])
y = np.array([0, 1, 0])
z = np.array([0, 0, 1])
# The trigonal geometry is defined by three vectors that are 120 degrees apart in the x-y plane.
tri1 = xy_angle(0)
tri2 = xy_angle(120)
tri3 = xy_angle(240)
# The trigonal prismatic geometry is defined by two parallel triangles sandwiching the metal center.
trip1 = tri1 - 0.5*z
trip2 = tri2 - 0.5*z
trip3 = tri3 - 0.5*z
trip4 = tri1 + 0.5*z
trip5 = tri2 + 0.5*z
trip6 = tri3 + 0.5*z
# The tetrahedral geometry is defined by four vectors that are 109.5 degrees apart in 3D space.
tet1 = np.array([1, 1, -1])
tet2 = np.array([-1, 1, 1])
tet3 = np.array([1, -1, 1])
tet4 = np.array([-1, -1, -1])
# The tetradentate 'mer_tridentate_like' geometry is defined by four atoms in the space of a mer-tridentate ligand with an angle of 60 degrees between each atom.
mer1 = xy_angle(0)
mer2 = xy_angle(60)
mer3 = xy_angle(120)
mer4 = xy_angle(180)
# The tetradentate 'monodentate like' geometry is defined by four planar atoms, all in a 45 degree towards the four bottom corners of a cube
tetmono1 = np.array([1, 0, -1])
tetmono2 = np.array([0, 1, -1])
tetmono3 = np.array([-1, 0, -1])
tetmono4 = np.array([0, -1, -1])
# The pentagonal geometry is defined by five vectors that are 72 degrees apart in the x-y plane.
penta1 = xy_angle(0*72)
penta2 = xy_angle(1*72)
penta3 = xy_angle(2*72)
penta4 = xy_angle(3*72)
penta5 = xy_angle(4*72)
# The hexagonal geometry is defined by six vectors that are 60 degrees apart in the x-y plane.
hex1 = xy_angle(0*60)
hex2 = xy_angle(1*60)
hex3 = xy_angle(2*60)
hex4 = xy_angle(3*60)
hex5 = xy_angle(4*60)
hex6 = xy_angle(5*60)

# Specify the name of the geometry, the target vectors and the weight of the geometry. The higher the weight, the less likely the geometry is chosen (default weight is 1.0).
all_geometries = {
    1:
        {
            '1_monodentate': ((x,), 1),
        },
    2:
        {
            '2_cis': ((x, y), 1),
            '2_trans': ((x, -x), 1)
        },
    3:
        {
            '3_meridional': ((y, x, -y), 1),
            '3_facial': ((x, y, z), 1.467),
            '3_trigonal': ((tri1, tri2, tri3), 1.045),
            # '3_trigonal_offset': ((trip1, trip2, trip3), 1.5),    # looks extremely close to a facial geometry so leave it out
        },
    4:
        {
            '4_tetrahedral': ((tet1, tet2, tet3, tet4), 1.06),
            '4_tetragonal': ((x, y, -x, -y), 0.96),
            '4_seesaw': ((x, y, z, -z), 1.45),
            '4_trigonal_pyramidal': ((tri1, tri2, tri3, z), 1.051),
            '4_compressed_meridional': ((mer1, mer2, mer3, mer4), 1.8),
            '4_tetragonal_offset': ((tetmono1, tetmono2, tetmono3, tetmono4), 1.13),
        },
    5:
        {
            '5_square_pyramidal': ((x, y, -x, -y, z), 1),
            # 'trigonal_bipyramidal': ((tri1, tri2, tri3, z, -z), 1),     # Doesn't seem to exist.
            '5_pentagonal': ((penta1, penta2, penta3, penta4, penta5), 0.8)
        },
    6:
        {
            '6_octahedral': ((x, y, z, -x, -y, -z), 1),
            '6_hexagonal': ((hex1, hex2, hex3, hex4, hex5, hex6), 1.30),
            '6_trigonal_prismatic': ((trip1, trip2, trip3, trip4, trip5, trip6), 1),
            '6_pentagonal_pyramidal': ((penta1, penta2, penta3, penta4, penta5, z), 1),
        },
    # For the other denticities, just make up one random geometry so that the code can assemble them in theory, but these should only be assembled on their own, without any other ligands, since they are so crowded.
    7:
        {
            '7_septadentate': ((x, y, z, -x, -y, -z, hex1), 1),                   # made up
        },
    8:
        {
            '8_octadentate': ((x, y, z, -x, -y, -z, hex1, hex2), 1)               # made up
        },
    9:
        {
            '9_nonadentate': ((x, y, z, -x, -y, -z, hex1, hex2, hex3), 1)         # made up
        },
    10:
        {
            '10_decadentate': ((x, y, z, -x, -y, -z, hex1, hex2, hex3, hex4), 1)   # made up
        },
}

def remove_mirrored_isomers_in_planar_ring_like_ligands(best_idc: list[list[int]]) -> list[list[int]]:
    """
    This algorithm affects only planar, ring-like geometries, specifically trigonal, tetragonal, pentagonal and hexagonal ligands. It removes isomers which are simply a rotation of the planar donor atoms around the z-axis. It detects which indices correspond to isomeres mirrored around the z-axis and and only keeps two of them.
    This algorithm only affects trigonal, tetragonal, pentaonal and hexagonal ligands per construction, since otherwise the constructed "mirrored indices" simply dont exist in the list and l1 and l2 will be empty. However, an exception are bidentate ligands, which would also be affected by this algorithm (since a bidentate ligand is "per definition" a planar, ring-like ligand). Therefore, it is important to apply this algorithm only to ligands with denticity >= 3.
    :param best_idc: List of the best indices of the donor atoms for each target vector
    :return: List of the best indices of the donor atoms for each target vector, with mirrored isomers removed.
    """
    denticity = len(best_idc[0])
    # Protect against affecting bidentate ligands
    if denticity <=2:
        return best_idc

    l1, l2 = [], []                                 # Fill with indices and mirrored indices
    for idc in best_idc:
        if idc not in l1 and idc not in l2:
            mirrored_idc = [idc[0]] + idc[1:][::-1]   # Flip the direction of going through the "ring" of donor atoms (e.g. clockwise->anti-clockwise), thereby mirroring the isomer.
            if mirrored_idc in best_idc:
                l1.append(idc)
                l2.append(mirrored_idc)
    if l1 and l2:   # Other geometries will simply pass because l1 and l2 will be empty.
        best_idc = [l1[0], l2[0]]   # Reduce to only two isomers which are mirrored with regard to each other.

    return best_idc

def get_geometrical_isomers_from_trying_out_all_possibilities(
                                    atoms: ase.Atoms,
                                    donor_idc: list[int],
                                    geometry: str
                                    ) -> tuple[list[ase.Atoms], list[list[int]], float]:

    # Some geometries have only one possible isomer, so we can skip the whole process of trying out all isomers. We still want to process the monodentate though, because otherwise it would not correspond to the target vector.
    if geometry in ['7_septadentate', '8_octadentate', '9_nonadentate', '10_decadentate']:
        best_rssd = np.nan
        best_isomers = [atoms]
        best_idc = [donor_idc]
        return best_isomers, best_idc, best_rssd


    # Try out all possible isomers and return all the ones with lowest rssd.
    target_vectors, weight = all_geometries[len(donor_idc)][geometry]
    best_isomers, best_idc, best_rssd = try_all_geometrical_isomer_possibilities(
                                                                                    atoms=atoms,
                                                                                    donor_idc=donor_idc,
                                                                                    target_vectors=target_vectors
                                                                                    )
    best_rssd *= weight

    # Remove symmetrically equivalent isomers for trigonal, tetragonal, pentagonal and hexagonal ligands.
    if geometry in ['3_trigonal', '4_tetragonal', '5_pentagonal', '6_hexagonal']:
        best_idc = remove_mirrored_isomers_in_planar_ring_like_ligands(best_idc)

    # For certain geometries, all generated isomers are equivalent, so we only keep one of them. For tetrahedral, octahedral and trigonal prismatic that is because there will be no other ligands. For square pyramidal, pentagonal pyramidal and trigonal pyramidal, that is because the other isomers are simply rotations of the same isomer around the z-axis (similar to the planar ring-like ligands, just with an additional atom on top).
    if geometry in ['4_tetrahedral', '4_trigonal_pyramidal', '5_square_pyramidal', '6_trigonal_prismatic', '6_octahedral',  '6_pentagonal_pyramidal']:
        best_idc = [best_idc[0]]

    # Apply the reduction of isomers to the isomers list.
    best_isomers = [best_isomers[best_idc.index(idc)] for idc in best_idc]

    return best_isomers, best_idc, best_rssd


def try_all_geometrical_isomer_possibilities(
                                                atoms: ase.Atoms,
                                                donor_idc: list[int],
                                                target_vectors: list[np.ndarray]
                                                ) -> tuple[list[ase.Atoms], list[list[int]], float]:
    """
    Tries out all combinations of mapping the provided donor atoms to the target vectors, for each target vector in the list of target vectors. All combinations are tried out and all isomers with the lowest rssd are returned. Usually, there will be several with the same rssd.
    :param atoms: ase.Atoms() object of the ligand
    :param donor_idc: Indices of the donor atoms in `atoms`
    :param target_vectors: List of 3D vectors of shape (n,3)
    :return: Tuple of:
        - List of ASE Atoms objects of the best isomers
        - List of lists of indices of the best isomers
        - The root sum of squared differences (RSSD) of the best isomers
    """
    target_vectors = np.array(target_vectors)
    if target_vectors.ndim == 1:
        target_vectors = target_vectors[None, :]
    assert target_vectors.shape[-1] == 3, f'Wrong dimension of target_vectors. It\'s {target_vectors.shape[-1]} but must be 3.'

    data = []
    n = len(target_vectors)
    donor_idc_permutations = list(itertools.permutations(donor_idc, n))
    for idc in donor_idc_permutations:
        idc = list(idc)   # Convert tuple to list to allow indexing later on with the indices saved in `data`
        isomer, rssd = align_donor_atoms(atoms, donor_idc=idc, target_vectors=target_vectors, return_rssd=True)
        data.append((rssd, idc, isomer))
    best_rssd, _, _ = min(data, key=lambda x: x[0])

    all_best_rssd = []
    all_best_idc = []
    all_best_isomers = []
    for rssd, vectors, isomer in data:
        if np.isclose(rssd, best_rssd):
            all_best_rssd.append(rssd)
            all_best_idc.append(vectors)
            all_best_isomers.append(isomer)
    best_rssd = np.mean(all_best_rssd)

    return all_best_isomers, all_best_idc, best_rssd

def assign_geometry(atoms: ase.Atoms, donor_idc: list[int]) -> tuple[str, list[ase.Atoms], float, str, float]:
    """
    Assigns the geometry of a ligand.
    :param atoms: ASE Atoms object of the ligand.
    :param donor_idc: Indices of the donor atoms in the ligand.
    :return: Tuple of:
        - The assigned geometry
        - List of ASE Atoms objects of the best isomers
        - The root sum of squared differences (RSSD) of the assigned geometry
        - The second best geometry
        - The weight necessary for a change in geometry
    """
    denticity = len(donor_idc)
    try:
        geometries = all_geometries[denticity]
    except KeyError:
        raise ValueError(f'No geometry defined for denticity {denticity}.')

    all_rssds = []
    for geometry in geometries.keys():
        isomers, best_idc, rssd = get_geometrical_isomers_from_trying_out_all_possibilities(
                                                                                            atoms=atoms,
                                                                                            donor_idc=donor_idc,
                                                                                            geometry=geometry
                                                                                            )
        all_rssds.append((geometry, rssd, isomers, best_idc))
    geometry, rssd, isomers, best_idc = min(all_rssds, key=lambda x: x[1])  # Find geometry with lowest rssd

    # Find second best geometry.
    if len(all_rssds) == 1:     # There is only one geometry, so there is no second best geometry.
        second_best_geometry, second_best_rssd = None, np.nan
    else:
        second_best_geometry, second_best_rssd, _, _ = min([x for x in all_rssds if x[0] != geometry], key=lambda x: x[1])

    # Calculate the weight necessary for a change in geometry. Useful to see how much the geometry would have to change to be the second best geometry.
    weight_necessary_for_change = second_best_rssd / rssd

    return geometry, isomers, rssd, second_best_geometry, weight_necessary_for_change

def align_donor_atoms(
                        atoms: ase.Atoms,
                        donor_idc: list[int],
                        target_vectors: list[list[float]],
                        return_rssd: bool = False
                        ):
    """
    Align the donor atoms of a ligand to the target vectors.
    :param atoms: ASE Atoms object of all the atoms of the ligand.
    :param donor_idc: Indices of the donor atoms of the ligand.
    :param target_vectors: A list of 3D target vectors to which the donor atoms should be aligned.
    :return: ASE Atoms object of the ligand with the donor atoms aligned to the target vectors.
    """
    atoms = atoms.copy()    # Make a copy of the original atoms object so that the original is not changed
    target_vectors = np.array(target_vectors)
    assert len(target_vectors) == len(donor_idc), 'The number of target vectors must match the number of donor atoms.'
    assert target_vectors.shape[1] == 3, 'The target vectors must be 3D vectors.'
    donor_idc = list(donor_idc)   # A tuple wouldn't work for indexing

    donor_vectors = atoms.positions[donor_idc]
    # Normalize the donor vectors and target vectors to unit vectors so that only the direction of the vectors counts, not the magnitude.
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




#
# if __name__ == '__main__':
#
#     n_max = 1500  # None or int. Maximum number of ligands to load
#     denticities = None  # None for all denticities or list of denticities
#     concat_outdir = 'concat_xyz'
#     plot_outdir = 'plots'
#     remove_haptic = False
#     sort_by_rssd = False
#     output_all_isomers = True
#
#
#
#     concat_outdir = Path(concat_outdir)
#     plot_outdir = Path(plot_outdir)
#     db = LigandDB.load_from_json(n_max=n_max)
#     if remove_haptic:
#         db.db = {name: ligand for name, ligand in db.db.items() if not ligand.has_neighboring_coordinating_atoms}
#     db = db.get_db_with_only_certain_denticities(denticities=denticities)
#
#     data = db.save_ligand_geometry_concat_xyz_files(
#                                             outdir=concat_outdir,
#                                             sort_by_rssd=sort_by_rssd,
#                                             output_all_isomers=output_all_isomers
#                                             )
#
#     # Plot distribution of rssd values
#     import seaborn as sns
#     import matplotlib.pyplot as plt
#     for geometry, names_rssd in data.items():
#         plt.figure()
#         _, rssds, _, _, _ = zip(*names_rssd)
#         sns.histplot(rssds, label=geometry, bins=15)
#         plt.title(geometry)
#         plt.xlabel('RSSD')
#         plt.ylabel('Count')
#         outpath = Path(plot_outdir, f'{geometry}.png')
#         plt.savefig(outpath)
#         plt.close()
#
#     # ==============    Doublecheck refactoring    ==================
#     from dev.test.Integration_Test import IntegrationTest
#     new_dir = concat_outdir
#     old_dir = concat_outdir.parent / f'OLD_n={n_max}_concat_xyz'
#     if old_dir.exists():
#         test = IntegrationTest(new_dir=new_dir, old_dir=old_dir)
#         test.compare_all()
#         print('Test for ligand geometries passed!')
#     else:
#         print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')
#
#     print('Done!')

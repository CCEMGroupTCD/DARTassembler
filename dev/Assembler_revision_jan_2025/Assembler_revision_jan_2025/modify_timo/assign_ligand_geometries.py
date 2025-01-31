from collections import defaultdict
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import ase
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
import numpy as np
import itertools
from tqdm import tqdm
from utilities import save_to_xyz
from utilities import align_donor_atoms

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
# The trigonal prismatic geometry is defined by the six following vectors.
trip1 = np.array([1.3, 2/3, -0.5])
trip2 = np.array([-1/6, 0.53269207, -0.5])
trip3 = np.array([-1/6, -1.19935874, -0.5])
trip4 = np.array([1.3, 2/3, 0.5])
trip5 = np.array([-1/6, 0.53269207, 0.5])
trip6 = np.array([-1/6, -1.19935874, 0.5])
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
# The tetradentate 'monodentate_like' geometry is defined by four planar atoms, all in a 45 degree towards the four bottom corners of a cube
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
            'monodentate': ((x,), 1),   # Cian: Is there a better term that also includes haptic ligands such as cp*? If yes, also use similarly for denticities >=7.
                                        # todo: I would propose organizing the ligand geometries by effective ligand coordination number (ELCN) instead of densiticity because then we group
                                        # "assembly-similar" ligands together. we can further annotate ligands with kappa and eta values for denticity and hapticity, respectively.
        },
    2:
        {
            'cis': ((x, y), 1),         # todo good cis bidentate
            'trans': ((x, -x), 1)       # todo good trans bidentate
        },
    3:
        {
            'mer': ((y+0.1*(x+z), x, -y+0.1*(x+z)), 1),             # todo meridional
            'fac': ((x, y, z), 1),                                  # todo facial
            'trigonal': ((tri1, tri2, tri3), 1),                    # todo trigonal
            'trigonal_pyramidal': ((tet1, tet2, tet3), 1.05),       # todo trigonal offset (tbh I feel this could be merged with the trigonal potentially) tridentate
        },
    4:
        {
            'tetrahedral': ((tet1, tet2, tet3, tet4), 1.06),        # todo tetrahedral
            'mer': ((x, y, -x, -y), 0.96),                          # todo tetragonal
            'fac': ((x, y, z, -z), 1.45),                           # todo seesaw
            'tripodal': ((tri1, tri2, tri3, z), 1.051),             # todo trigonal pyramidal
            'tridentate_like': ((mer1, mer2, mer3, mer4), 1.8),     # todo compressed meridional
            'square_antiprismatic': ((tetmono1, tetmono2, tetmono3, tetmono4), 1.13),  # todo tetragonal offset
        },
    5:
        {
            'square_pyramidal': ((x, y, -x, -y, z), 1),                     # todo square pyramidal
            # 'trigonal_bipyramidal': ((tri1, tri2, tri3, z, -z), 1),     # doesn't seem to exist
            'pentagonal': ((penta1, penta2, penta3, penta4, penta5), 0.8)   # todo pentagonal
        },
    6:
        {
            'octahedral': ((x, y, z, -x, -y, -z), 1),                       # todo octahedral
            'hexagonal': ((hex1, hex2, hex3, hex4, hex5, hex6), 0.91),      # todo hexagonal
            'trigonal_prismatic': ((trip1, trip2, trip3, trip4, trip5, trip6), 1), # todo good
            'pentagonal_pyramidal': ((penta1, penta2, penta3, penta4, penta5, z), 1), # todo good
        },
    # For the other denticities, just make up one random geometry so that the code can assemble them in theory,  but these should only be assembled on their own, without any other ligands, since they are so crowded.
    7:
        {
            'septadentate': ((x, y, z, -x, -y, -z, hex1), 1),                   # made up   # todo are there any  completely planar ones?
        },
    8:
        {
            'octadentate': ((x, y, z, -x, -y, -z, hex1, hex2), 1)               # made up
        },
    9:
        {
            'nonadentate': ((x, y, z, -x, -y, -z, hex1, hex2, hex3), 1)         # made up
        },
    10:
        {
            'decadentate': ((x, y, z, -x, -y, -z, hex1, hex2, hex3, hex4), 1)   # made up
        },
}

def try_all_combinations_of_target_vectors(donor_atoms: ase.Atoms, target_vectors: list[np.ndarray]) -> list[np.ndarray]:
    data = []
    ligand_idc = list(range(len(donor_atoms)))
    n = len(target_vectors)
    vector_permutations = itertools.permutations(target_vectors, n)
    for vector_combo in vector_permutations:
        _, rssd = align_donor_atoms(donor_atoms, ligand_idc=ligand_idc, target_vectors=vector_combo, return_rssd=True)
        data.append((rssd, vector_combo))
    best_rssd, best_vectors = min(data, key=lambda x: x[0])

    return best_rssd, best_vectors

def assign_geometry(donor_atoms: ase.Atoms):
    """
    Assigns the geometry of a ligand.
    :param donor_atoms: ASE atoms object of the ligand donor atoms, assuming the metal center was at (0,0,0).
    :return: Geometry of the ligand as a string
    """
    denticity = len(donor_atoms)
    try:
        geometries = all_geometries[denticity]
    except KeyError:
        raise ValueError(f'No geometry defined for denticity {denticity}.')

    if len(geometries) == 1:    # Just return the only geometry since there is no choice
        geometry, (vectors, _) = list(geometries.items())[0]
        rssd = np.nan
        return geometry, rssd, vectors, None, np.nan
    else:
        all_rssds = []
        for geometry, (target_vectors, weight) in geometries.items():
            rssd, vectors = try_all_combinations_of_target_vectors(donor_atoms, target_vectors)
            rssd *= weight
            all_rssds.append((geometry, rssd, vectors))
        # Find geometry with lowest rssd
        geometry, rssd, vectors = min(all_rssds, key=lambda x: x[1])
        # Find second best geometry
        second_best_geometry, second_best_rssd, _ = min([x for x in all_rssds if x[0] != geometry], key=lambda x: x[1])

    geometry = f'{denticity}{geometry}'
    weight_necessary_for_change = second_best_rssd / rssd

    return geometry, rssd, vectors, second_best_geometry, weight_necessary_for_change



if __name__ == '__main__':

    n_max = None
    denticities = None  # None for all denticities or list of denticities
    concat_outdir = 'concat_xyz'
    plot_outdir = 'plots'
    remove_haptic = True

    db = LigandDB.load_from_json(n_max=n_max)
    if remove_haptic:
        db.db = {name: ligand for name, ligand in db.db.items() if not ligand.has_neighboring_coordinating_atoms}
    data = defaultdict(list)
    for name, ligand in tqdm(db.db.items(), total=len(db.db)):
        donor_atoms = ligand.get_effective_donor_atoms()
        if denticities is not None and len(donor_atoms) not in denticities:
            continue
        geometry, rssd, vectors, second_geometry, weight_necessary_for_change = assign_geometry(donor_atoms)
        data[geometry].append((name, rssd, vectors, weight_necessary_for_change, second_geometry))
    # Sort by weight necessary for change
    for geometry, names_rssd in data.items():
        names_rssd.sort(key=lambda x: x[3])

    # Save structures for each geometry in a different concatenated xyz file
    for geometry, names_rssd in data.items():
        names, rssds, _, weights, second_geometries = zip(*names_rssd)
        atoms = [db.db[name].get_ase_molecule_with_metal('Fe') for name in names]
        # Round weight always up and to three decimals
        weights = [np.ceil(weight*1000)/1000 for weight in weights]
        comments = [f'{name} rssd={rssd:.3f} change:{weight:.3f}->{second_geometry}' for name, rssd, vectors, weight, second_geometry in names_rssd]
        outpath = Path(concat_outdir, f'{geometry}.xyz')
        print(f'{geometry}: {len(atoms)} structures')
        save_to_xyz(outpath=outpath, structures=atoms, comments=comments)

    # Plot distribution of rssd values
    for geometry, names_rssd in data.items():
        plt.figure()
        _, rssds, _, _, _ = zip(*names_rssd)
        sns.histplot(rssds, label=geometry, bins=15)
        plt.title(geometry)
        plt.xlabel('RSSD')
        plt.ylabel('Count')
        outpath = Path(plot_outdir, f'{geometry}.png')
        plt.savefig(outpath)
        plt.close()

    print('Done!')

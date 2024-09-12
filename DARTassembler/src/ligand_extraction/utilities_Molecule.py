import collections
import networkx as nx
from rdkit import Chem
import numpy as np
from DARTassembler.src.constants.Periodic_Table import DART_Element
import re

unknown_rdkit_bond_orders = [0, 20, 21]


def get_rdkit_mol_from_smiles(smiles: str, sanitize: bool=False) -> Chem.Mol:
    """
    Get an RDKit molecule from a SMILES string. If sanitize is set to False, the molecule will not be sanitized but an attempt will be made to calculate properties important for other functions.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
    if not sanitize:
        # Let rdkit calculate properties important for other functions, especially HasSubstructMatch(). Usually, those properties would be called by the sanitization process, but here we need to do that manually so that the right properties are set.
        Chem.FastFindRings(mol)
        mol.UpdatePropertyCache(strict=False)

    return mol

def has_smarts_pattern(smarts: str, smiles: str) -> bool:
    """
    Checks whether the molecule matches the given SMARTS pattern.
    @param smarts: SMARTS pattern to match.
    @return: True if the molecule matches the SMARTS pattern, False otherwise.
    """
    mol = get_rdkit_mol_from_smiles(smiles)

    # Check if the molecule matches the SMARTS pattern.
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f'Invalid SMARTS pattern: {smarts}')
    match = mol.HasSubstructMatch(pattern)

    return match

def stoichiometry2atomslist(stoichiometry: str) -> list[str]:
    """
    Convert a stoichiometry string to a list of atoms. For example, "C2H6F" becomes ["C", "C", "H", "H", "H", "H", "H", "H", "F"]. Note that the order of the atoms is not necessarily the same as in the stoichiometry string.
    """
    atoms_list = []
    # Use a regular expression to find all elements and their counts
    for element, count in re.findall(r'([A-Z][a-z]*)(\d*)', stoichiometry):
        # If no count is given, assume it's 1
        count = int(count) if count else 1
        # Extend the list by count times the element
        atoms_list.extend([element] * count)
    return atoms_list

def if_same_stoichiometries(stoichiometry1: str, stoichiometry2: str) -> bool:
    """
    Check if two stoichiometries are the same. Stable against different order of elements in the stoichiometry string and the convention of writing or not writing a count of 1.
    """
    atoms_list1 = sorted(stoichiometry2atomslist(stoichiometry1))
    atoms_list2 = sorted(stoichiometry2atomslist(stoichiometry2))

    return atoms_list1 == atoms_list2


def find_smallest_ring_with_specified_nodes(graph, nodes):
    """
    Returns the indices of nodes in the smallest ring of a networkx graph containing the specified nodes.

    :param graph: networkx graph object
    :param nodes: List of node indices. We're looking for rings containing these nodes.
    :return: The smallest ring represented as a list of node indices, or None if no such ring is found.
    """

    # Get the list of all cycles in the graph
    cycles = nx.cycle_basis(graph)

    # Filter the cycles to include only those containing the specified nodes
    valid_rings = [cycle for cycle in cycles if set(nodes).issubset(set(cycle))]

    # If there are no valid rings, return None
    if not valid_rings:
        return None

    # Sort the valid rings based on their length and return the smallest one
    return sorted(valid_rings, key=len)[0]

def get_max_deviation_from_coplanarity(points: list[tuple]) -> float:
    """
    Check if a list of 3D coordinates are approximately coplanar based on a specified cutoff.

    :param points: List of tuples, each tuple being a 3D point (x, y, z).
    :return: deviation from coplanarity. 0 if perfectly coplanar and > 0 otherwise.
    """
    # We need at least 3 points to define a plane, everything else is per definition coplanar
    if len(points) <= 3:
        return 0

    # Center the points
    points_np = np.array(points)
    centroid = np.mean(points_np, axis=0)
    centered_points = points_np - centroid

    # Compute the SVD
    _, s, _ = np.linalg.svd(centered_points)

    # The smallest singular value is the maximum distance from the plane
    max_dist = s[-1]

    return max_dist

def are_points_coplanar(points, dist=0.1):
    """
    Check if a list of 3D coordinates are approximately coplanar based on a specified cutoff.

    :param points: List of tuples, each tuple being a 3D point (x, y, z).
    :param cutoff: Maximum allowed distance of a point from the plane to be considered coplanar.
    :return: True if points are approximately coplanar, False otherwise.
    """
    max_dist = get_max_deviation_from_coplanarity(points)

    return max_dist <= dist

def graph_to_rdkit_mol(graph: nx.Graph, element_label: str='node_label', bond_label: str='bond_type') -> Chem.Mol:
    """
    DEPRECATED: Use get_rdkit_mol_from_smiles() instead.

    Create an rdkit mol object from a graph. Note that the bond type must be specified in the graph under the attribute called `edge_label`.
    @param graph: input graph of the molecule
    @param element_label: element label for the node dictionary
    @param bond_label: bond type label for the edge dictionary
    @return: rdkit mol object
    """

    # create empty editable mol object
    mol = Chem.RWMol()

    # add atoms to mol and keep track of index
    node_to_idx = {}
    for idx, atom in graph.nodes(data=True):
        a = Chem.Atom(atom[element_label])
        molIdx = mol.AddAtom(a)
        node_to_idx[idx] = molIdx

    # add bonds between adjacent atoms
    for idx1, idx2, bond in graph.edges(data=True):
        bond = bond[bond_label]

        try:
            bond_type = Chem.rdchem.BondType.values[bond]
        except AttributeError:
            ValueError(f'Unknown bond type {bond} in molecule.')

        mol.AddBond(node_to_idx[idx1], node_to_idx[idx2], bond_type)

    # Convert RWMol to Mol object
    mol = mol.GetMol()

    assert mol.GetNumAtoms() == len(graph), f'Number of atoms in rdkit molecule ({mol.GetNumAtoms()}) does not match number of atoms in molecule ({len(graph)})'
    return mol

def rdkit_mol_to_graph(mol: Chem.Mol, element_label: str='node_label', bond_label: str= 'bond_type') -> nx.Graph:
    """
    Create a graph from an rdkit mol object
    @param mol: RDKit mol object
    @param element_label: element label for node dictionary
    @param bond_label: bond type label for edge dictionary
    @return: networkx graph of the molecule
    """
    G = nx.Graph()

    for atom in mol.GetAtoms():
        node = atom.GetIdx()
        label = atom.GetSymbol()
        G.add_node(node, **{element_label: label})

    for bond in mol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()
        label = bond.GetBondTypeAsDouble()
        G.add_edge(u, v, **{bond_label: label})

    return G

def get_all_ligands_by_graph_hashes(all_ligands: list) -> dict:
    """
    Get dictionary of graph hashes with list of ligands with this graph hash
    :param all_ligands: list of all ligands with graph hashes
    :return: dictionary of graph hash: list_of_ligands
    """
    all_hashes = list(set([lig.graph_hash for lig in all_ligands]))
    all_ligands_by_hashes = {h: [] for h in all_hashes}
    for ligand in all_ligands:
        all_ligands_by_hashes[ligand.graph_hash].append(ligand)

    return all_ligands_by_hashes

def calculate_angular_deviation_of_bond_axis_from_ligand_center(atomic_props: dict, metal_position: tuple, donor_indices: list, use_center_of_mass: bool=False) -> float:
    """
    Calculate for monodentate ligands the angular deviation of the bond axis from the ligand center or center of mass.
    :param atomic_props: Dictionary containing atom types and their coordinates.
    :param metal_position: Coordinates of the metal center.
    :param donor_indices: Indices of atoms involved in coordination to the metal center.
    :param use_center_of_mass: Whether to use the center of mass or the center of points.
    :return: Angular deviation from centrosymmetry for monodentates, np.nan for others.
    """
    # Check if monodentate
    if len(donor_indices) != 1:
        return np.nan

    x_centroid_list = np.array(atomic_props['x'])
    y_centroid_list = np.array(atomic_props['y'])
    z_centroid_list = np.array(atomic_props['z'])

    # Calculate centroid
    if use_center_of_mass:
        atoms_list = atomic_props['atoms']
        masses = np.array([DART_Element(atom).atomic_mass for atom in atoms_list])
        x_centroid = np.sum(masses * x_centroid_list) / np.sum(masses)
        y_centroid = np.sum(masses * y_centroid_list) / np.sum(masses)
        z_centroid = np.sum(masses * z_centroid_list) / np.sum(masses)
    else:
        x_centroid = np.mean(x_centroid_list)
        y_centroid = np.mean(y_centroid_list)
        z_centroid = np.mean(z_centroid_list)

    centre_of_points = np.array([x_centroid, y_centroid, z_centroid]) - np.array(metal_position)
    coord_atom_pos = np.array([x_centroid_list[donor_indices[0]],
                               y_centroid_list[donor_indices[0]],
                               z_centroid_list[donor_indices[0]]]) - np.array(metal_position)

    v1 = coord_atom_pos
    v2 = centre_of_points
    v3 = v2 - v1
    if np.all(v1 != v2):
        cosine = np.dot(v1 * (-1), v3) / (np.linalg.norm(v1 * (-1)) * np.linalg.norm(v3))
        cosine = np.clip(cosine, -1, 1)  # Clip the cosine value to avoid ValueError in np.arccos
        angle = np.arccos(cosine)
        angle = np.degrees(angle)
    else:
        angle = 180

    min_ang_dev = min((abs(180 - angle), abs(0 - angle)))

    return min_ang_dev

def group_list_without_hashing(ligand_list: list) -> list:
    """
    Returns a list of list with unique elements grouped together. Works without hashing, just using equity.
    :param ligand_list: list of elements
    :return: list of lists of grouped elements
    """
    groupings = {}
    counter = 0
    for lig1 in ligand_list:

        tmp_groupings = {}
        equal = False

        for i, lig_list in groupings.items():
            lig_representative = lig_list[0]
            equal = lig_representative == lig1

            if equal:

                if i in tmp_groupings:
                    tmp_groupings[i].append(lig1)
                else:
                    tmp_groupings[i] = lig1

                break

        if not equal:
            tmp_groupings[counter] = [lig1]
            counter += 1

        for i in tmp_groupings.keys():

            if i in groupings:
                groupings[i].append(tmp_groupings[i])
            else:
                groupings[i] = tmp_groupings[i]

    groupings = [group for group in groupings.values()]
    return groupings


def original_metal_ligand(ligand):
    """
    We try to find the original metal of a ligand and return None if we couldnt find any
    """

    if hasattr(ligand, "original_metal"):
        return ligand.original_metal_symbol
    elif ligand.global_props is not None:
        try:
            return ligand.global_props['metal_name']
        except KeyError:
            pass
    else:
        return None


def get_standardized_stoichiometry_from_atoms_list(atoms: list) -> str:
    c = collections.Counter(atoms)
    elements = sorted(el for el in c.keys())
    if "C" in elements:
        if 'H' in elements:
            elements = ["C", 'H'] + [el for el in elements if el not in ["C", 'H']]
        else:
            elements = ["C"] + [el for el in elements if el != "C"]
    else:
        if 'H' in elements:
            elements = ["H"] + [el for el in elements if el != "H"]

    # formula = [f"{el}{(c[el]) if c[el] != 1 else ''}" for el in elements] # drop the 1 if an element occurs only once
    formula = [f"{el}{(c[el])}" for el in elements]
    return "".join(formula)


def get_concatenated_xyz_string_from_coordinates(coord_list: list[np.array], element_list: list[list[str]], comment: str = "") -> str:
    """
    Returns a string in xyz format from a list of numpy arrays of coordinates and a list of lists of chemical elements.
    @param coord_list: A list of numpy arrays of coordinates.
    @param element_list: A list of lists of chemical elements.
    @param comment: A comment to be added to the xyz file.
    @return: A string in xyz format.
    """
    if isinstance(coord_list, np.ndarray) and coord_list.ndim == 2: # If only one set of coordinates is given
        coord_list = [coord_list]
    if isinstance(element_list[0], str):
        element_list = [element_list]

    xyz = ''
    for coord, elements in zip(coord_list, element_list):
        xyz += xyz_string_from_coordinates(coord=coord, elements=elements, comment=comment)

    return xyz

def xyz_string_from_coordinates(coord: np.ndarray, elements: list[str], comment: str = "") -> str:
    """
    Returns a string in xyz format from a numpy array of coordinates and a list of chemical elements.
    """
    if coord.shape[1] != 3:
        raise ValueError(f"The coordinates do not have the correct shape: {coord.shape}")
    if len(elements) != coord.shape[0]:
        raise ValueError(
            f"The number of elements does not match the number of coordinates: {len(elements)} != {coord.shape[0]}")

    n_atoms = coord.shape[0]
    xyz = f"{n_atoms}\n"
    xyz += comment + '\n'
    for el, (x, y, z) in zip(elements, coord):
        xyz += f"{el} {x} {y} {z}\n"

    return xyz
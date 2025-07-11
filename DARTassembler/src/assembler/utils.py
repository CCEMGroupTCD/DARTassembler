import hashlib
import random
import ase
import numpy as np
from typing import List, Tuple
from collections import defaultdict
import itertools

def assign_ligands_to_vectors(ligands: List['Ligand'],
                               swap_groups: List[int]) -> List[List['Ligand']]:
    """
    Assigns ligands to vector entries based on swap group IDs.
    Ligands in the same swap group can be swapped among vector entries assigned the same swap group ID.

    :param ligands:       A list of ligand objects.
    :param swap_groups:   A list of integers where the i-th value defines the swap group for the i-th ligand.
                          Ligands with the same group ID can be swapped; those with different IDs cannot.
                          This list must match the length and order of the `vectors` list.
    :return:              A list of ligand groupings for each vector entry.
                          Each entry contains the swappable ligands for that vector,
                          used later for permutation generation.
    :raise:               LoggedValueError if swap_groups are inconsistent with ligand or vector count.
    """
    if swap_groups is None:
        return [ligands]

    n = len(ligands)
    if len(swap_groups) != n:
        raise ValueError("swap_groups must match the number of ligands and vectors.")

    # Group ligand indices by swap group
    group_to_indices = defaultdict(list)
    for idx, grp in enumerate(swap_groups):
        group_to_indices[grp].append(idx)

    # For each vector position, gather allowed ligand indices
    allowed = [group_to_indices[grp] for grp in swap_groups]

    # Build all assignments, filtering out any that reuse the same ligand twice
    results = []
    for combo in itertools.product(*allowed):
        if len(set(combo)) != n:
            continue
        results.append([ligands[i] for i in combo])

    assert results[0] == ligands, "First result must be the original ligand order."

    return results

def are_atoms_equal(atom1: ase.Atom, atom2: ase.Atom) -> bool:
    """
    Compares two ase.Atoms objects to see if they are equal.
    """
    return atom1.symbol == atom2.symbol and np.allclose(atom1.position, atom2.position)

def generate_pronounceable_word(length=5, seed=None, start_with_vowel=None) -> str:
    """
    Generate a pronounceable word by alternating vowels and consonants.

    Parameters:
    length (int): The length of the word to generate. Default is 5.
    seed (int): A seed for the random number generator. This can be used to generate the same word
    repeatedly. Default is None, which means the random number generator is not seeded.
    start_with_vowel (bool): Whether the word should start with a vowel. If None, this is chosen randomly.

    Returns:
    str: The generated word, in uppercase.
    """
    if not isinstance(seed, int):
        seed = int(hashlib.md5(str(seed).encode(encoding='UTF-8', errors='strict')).hexdigest(), 16)

    # Create a new random number generator and seed it if a seed is provided
    rng = random.Random(seed)

    vowels = 'aeiou'
    consonants = 'bcdfghjklmnpqrstvwxyz'
    word = ''

    # Start with a random choice between starting with a vowel or a consonant
    if start_with_vowel is None:
        start_with_vowel = rng.choice([True, False])

    for i in range(length):
        if (i % 2 == 0 and start_with_vowel) or (i % 2 == 1 and not start_with_vowel):
            word += rng.choice(vowels)
        else:
            word += rng.choice(consonants)

    return word.upper()


if __name__ == '__main__':
    n = range(10)
    standard = tuple(f'{generate_pronounceable_word(length=8, start_with_vowel=True)}' for _ in n)
    print(standard)
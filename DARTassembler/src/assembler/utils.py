import hashlib
import random
import ase
import numpy as np
from typing import List, Tuple, Optional
from collections import defaultdict
import itertools


def get_complex_name(seed: str, length: int, suffix: str = '', avoid_names: Optional[set[str]] = None) -> str:
    """
    Generates a unique name for the complex based on the graph hash and the specified length and suffix.
    :param avoid_names: A set of names to avoid. If the generated name is in this set, it will generate a new name with one more character.
    :return: Name (str) of the complex, which is pronounceable and unique.
    """
    if avoid_names is None:
        avoid_names = []
    while True:  # emulate a do-while loop
        # Generate a random, pronounceable name
        name = generate_pronounceable_word(length=length, seed=seed)
        name += suffix  # Add the suffix to the name

        # If the name is already used, redo name generation with one more character. For the next complex, it starts with the original character length again.
        if name in avoid_names:
            length += 1
            continue
        else:
            break  # name is unique, break the loop

    return name

def remove_haptic_dummy_atom(atoms: ase.Atoms, dummy_atom: str) -> ase.Atoms:
    """
    Removes the dummy atom from the generated isomers.
    :param atoms: Ase.Atoms object containing the structure.
    :param dummy_atom: The symbol of the dummy atom to remove (e.g., "Cu").
    :return:
    """
    dummy_idc = [i for i, atom in enumerate(atoms) if atom.symbol == dummy_atom]
    dummy_idc.sort(
        reverse=True)  # This is important so that the larger index is removed first so as not to change the index of the other atoms
    for dummy_idx in dummy_idc:
        atoms.pop(dummy_idx)
    return atoms

def get_list_with_all_possible_swappings(objects: list, swap_groups: List[int]) -> list[list]:
    """
    Returns a list of all possible combinations of the objects in `object_list` based on the provided `swap_groups`. Each group in `swap_groups` indicates which objects can be swapped with each other.
    :param objects: A list of objects to be swapped.
    :param swap_groups: A list of integers where each integer represents a group index. Objects with the same group index can be swapped with each other.
    :return: A list of lists, where each inner list is a unique combination of the objects in `object_list` based on the swap groups.
    """
    if swap_groups is None or len(set(swap_groups)) == 1:
        return [objects]

    n = len(objects)
    if len(swap_groups) != n:
        raise ValueError("`swap_groups` must match the length of `objects`.")

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
        results.append([objects[i] for i in combo])

    assert results[0] == objects, "First result must be the original ligand order."

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
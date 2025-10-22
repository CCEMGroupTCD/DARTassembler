import hashlib
import random
import ase
import numpy as np
from typing import List, Tuple, Optional, Any, Iterable
from collections import defaultdict
import itertools


def get_complex_name(seed: str, length: int, suffix: str = '', avoid_names: Optional[Iterable[str]] = None) -> str:
    """
    Generates a unique name for the complex based on the graph hash and the specified length and suffix.
    :param avoid_names: A set of names to avoid. If the generated name is in this set, it will generate a new name with one more character.
    :return: Name (str) of the complex, which is pronounceable and unique.
    """
    if avoid_names is None:
        avoid_names = set()
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

def join_duplicate_groups_by_union(
        pre_isomers_duplicate_group_names: Iterable[Iterable[str]],
        post_isomers_duplicate_group_names: Iterable[Iterable[str]],
        mode: str,
) -> List[List[str]]:
    """
    Join duplicate groups from the pre- and post-AxialOpt stages.

    Modes
    -----
    - "or"   (default): final duplicates if pre OR post grouped them (transitive closure).
                       → Most aggressive merge (your current behavior).
    - "pre"            : final duplicates follow ONLY the pre groups (post ignored).
                       → Precedent to pre stage (useful to preserve pre-stage distinctions).
    - "post"           : final duplicates follow ONLY the post groups (pre ignored).
                       → Treat AxialOpt as canonical; simplest final equivalence.

    Returns
    -------
    List[List[str]] : Disjoint duplicate groups, each sorted by name.
                      All names present in the chosen stage(s) appear exactly once.
    """
    # --- Normalize inputs: turn each incoming group into a set ---
    pre = [set(g) for g in pre_isomers_duplicate_group_names]
    post = [set(g) for g in post_isomers_duplicate_group_names]

    # Choose which groups to use based on the mode
    m = mode.strip().lower()
    if m == "or":
        source_groups = pre + post
    elif m == "pre":
        source_groups = pre
    elif m == "post":
        source_groups = post
    else:
        raise ValueError(f"_join_duplicate_groups_by_union: unknown mode '{mode}'. Use 'or' | 'pre' | 'post'.")

    # Collect every distinct isomer name seen in the selected groups.
    all_names = sorted({n for g in source_groups for n in g})

    # Edge case: no names (e.g., empty inputs)
    if not all_names:
        return []

    # --- Disjoint Set Union (Union-Find) setup ---
    parent = {n: n for n in all_names}  # representative for each name
    rank = {n: 0 for n in all_names}  # tree rank (for union by rank)

    def find(x: str) -> str:
        """Find the representative of x's set (path compression)."""
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        """Union sets of a and b (union by rank)."""
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    # --- Add edges: union all members within each selected group ---
    # For each group, tie everyone to an anchor so they land in one component.
    for g in source_groups:
        if not g:
            continue
        it = iter(g)
        anchor = next(it)
        for n in it:
            union(anchor, n)

    # --- Collect connected components as final groups ---
    comps: dict[str, List[str]] = {}
    for n in all_names:
        r = find(n)
        comps.setdefault(r, []).append(n)

    joined = [sorted(v) for v in comps.values()]

    # Defensive sanity: ensure disjoint partition over all seen names
    flat = [n for grp in joined for n in grp]
    assert len(flat) == len(all_names) == len(set(flat)), (
        "Joined groups must be a disjoint partition of all isomers seen in the chosen stage(s)."
    )

    return joined

if __name__ == '__main__':
    n = range(10)
    standard = tuple(f'{generate_pronounceable_word(length=8, start_with_vowel=True)}' for _ in n)
    print(standard)
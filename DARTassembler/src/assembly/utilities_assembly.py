import stk
import os
import random

def visualize(input_complex):
    # This method allows to visualize in a blocking way during debug but is not essential at all
    print("initializing visualization")
    stk.MolWriter().write(input_complex, 'input_complex.mol')
    os.system('obabel .mol input_complex.mol .xyz -O  output_complex.xyz')
    os.system("ase gui output_complex.xyz")
    os.system("rm -f input_complex.mol")
    os.system("rm -f output_complex.xyz")
    print("visualization complete")

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
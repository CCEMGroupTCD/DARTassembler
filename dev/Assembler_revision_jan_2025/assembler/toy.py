#########################################################################################
# This is a toy example to demonstrate how to remove duplicates from a list of floats   #
# Timo please ignore this code as it is not relevant to the assembler                   #
#########################################################################################
from typing import List

import ase

isomers = [4.00, 1.00, 2.00, 1.00, 1.00, 3.01, 4.01, 4.9, 11.5, 32.0, 6.0, 14.0, 14.01, 14.0001]
def compare(val1, val2):
    return abs(val1 - val2)


def _reduce_isomers(isomers, threshold) -> List[str]:
    """
    takes a list of 'isomers' and identifies which are super-imposable and removes duplicates
    :return: List[Atoms]
    """
    unique_isomers = []
    for isomer in isomers:
        # Check if this isomer is already in the list within the threshold
        if not any(compare(isomer, u) < threshold for u in unique_isomers):
            unique_isomers.append(isomer)  # Add only if it is not a duplicate

    return unique_isomers


# output = [1.0, 2.0, 3.01, 4.01, 4.9]

a = _reduce_isomers(isomers, threshold=0.1)
a.sort()
print(a)


from ase import Atoms

# Define two separate Atoms objects with different info attributes
atoms1 = Atoms('H2', positions=[[0, 0, 0], [0.7, 0, 0]])
atoms1.info['1'] = 1

atoms2 = Atoms('O2', positions=[[1.5, 0, 0], [2.2, 0, 0]])
atoms2.info['cool'] = {0: ['oxygen'], 1: ['oxygen', 'terminal']}
atoms2.info["2"] = 2

# Merge the two Atoms objects
#merged_atoms = atoms1 + atoms2
atoms1.append(atoms2)

# Check the .info attribute after merging
print(atoms1.info)  # Only retains atoms1's info!

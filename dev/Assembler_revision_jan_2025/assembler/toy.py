#########################################################################################
# This is a toy example to demonstrate how to remove duplicates from a list of floats   #
# Timo please ignore this code as it is not relevant to the assembler                   #
#########################################################################################
from typing import List

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

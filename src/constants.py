import os
import numpy as np      # Does all linear algebra
import sympy as sym         # Symbolic library to make "lambda-like" functions
import matplotlib.pyplot as plt         # Matplotlib imported for plotting
import math
import pandas as pd

'''

some basic chemistry stuff in here
can fully be replaced by the mendeleev package
'''


class element_in_pse:

    def __init__(self, symbol, molecular_mass, name, number):
        self.symbol = symbol
        self.molecular_mass = molecular_mass
        self.name = name
        self.number = number


def create_pse() -> dict:

    pse = dict()

    with open(f"chemical_data/pse.txt", "r") as f:

        while True:
            line1 = f.readline()
            if line1 == "":
                break
            words = line1.split()
            if len(words) == 0:
                continue             #leere Spalten überspringen

            symbol = words[2]
            M = words[3]
            name = words[1]

            pse[int(words[0])] = element_in_pse(symbol=symbol, molecular_mass=M, name=name, number=int(words[0]))

    return pse


PSE = create_pse()
metals_in_pse = [el for a in [[21, 31], [39, 49], [57, 81], [89, 113]] for el in range(a[0], a[1])]


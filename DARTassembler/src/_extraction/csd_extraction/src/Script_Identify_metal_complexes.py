"""
This is the script to identify monometallic TransitionMetalComplexes in the CSD
"""
from ccdc import io
from ccdc.entry import Entry
from mendeleev import element
import json
from tqdm import tqdm
from collections import Counter


class MonoMetallicComplexFilter:

        def __init__(self):

            # probably not the most elegany way, but this serves as an identifier if something is a TM
            self.trans_metals_in_pse = [el for a in [[21, 31], [39, 49],
                                                [57, 81], [89, 113]] for el in range(a[0], a[1])]
            self.metal_symbols = [element(i).symbol for i in self.trans_metals_in_pse]

            # then we need to init some variables
            self.csd_reader = io.EntryReader("CSD")
            self.list_of_identifiers = []

        def has_excatly_one_TM(self, a: Entry) -> bool:
            """
            We check if an entry in the CSD contains exactly one TM
            :param a:
            :return:
            """
            atom_list = [b.atomic_symbol for b in a.molecule.atoms]

            sum_total_metals = 0

            for key, item in Counter(atom_list).items():
                if key in self.metal_symbols:
                    sum_total_metals += int(item)

            return sum_total_metals == 1

        def run_search(self):

            # first we need to get only the complexes with only one metal
            for a in tqdm(self.csd_reader.entries(), desc="Finding Entries with exactly one metal"):
                if self.has_excatly_one_TM(a) is True:
                    self.list_of_identifiers.append(a.identifier)

            # This is kind of the minimum list, we can now filter it down further


if __name__ == "__main__":

    Filter = MonoMetallicComplexFilter()

    Filter.run_search()

    with open("../data/CSD_MM_identifier.json", "w") as file:
        json.dump(Filter.list_of_identifiers, file)

    print(f"Identifier list: {Filter.list_of_identifiers}")



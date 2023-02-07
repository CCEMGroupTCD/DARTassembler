import random
import time
import mendeleev as atom

random.seed(time.time())
# input total ligand charge
# output metal ox_state and multiplicity
# Temp Reference for allowed oxidation states: https://chem.libretexts.org/Bookshelves/Inorganic_Chemistry/Supplemental_Modules_and_Websites_(Inorganic_Chemistry)/Descriptive_Chemistry/Elements_Organized_by_Block/3_d-Block_Elements/1b_Properties_of_Transition_Metals/Oxidation_States_of_Transition_Metals#:~:text=The%20oxidation%20state%20of%20an,electrons)%20other%20atoms%20or%20species.
# used to check some results V.Rough: https://homework.study.com/explanation/calculate-the-spin-multiplicity-for-each-of-the-following-octahedral-complex-ions-and-give-the-ground-state-term-symbol-for-the-complex-ion-a-co-nh-3-6-3-plus-b-fe-cn-6-3-c-cr-nh-3-6-3-plus.html
# Multiplicity is 2S+1
metals_of_interest = {"Sc": [3],
                      "Ti": [4],
                      "V": [2, 3, 4, 5],
                      "Cr": [2, 3, 6],
                      "Mn": [2, 3, 4, 6, 7],
                      "Fe": [2, 3],
                      "Co": [2, 3],
                      "Ni": [2],
                      "Cu": [1, 2],
                      "Zn": [2]}  # Nickel 2+ generally only adopts a square planar configuration.

spins = {0: [1],
         1: [2],
         2: [3],
         3: [4],
         4: [3, 5],
         5: [2, 6],
         6: [1, 5],
         7: [2, 4],
         8: [3],
         9: [2],
         10: [1]}


def ox_spin(lig_charge, Dict, spin):
    tmp_metal_ox_storage = []
    for item in Dict:
        metal = item
        ox_list = Dict[item]
        for ox in ox_list:
            if int(ox) + int(lig_charge) == 0:
                tmp_metal_ox_storage.append([metal, ox])
            else:
                pass
    index = random.randint(0, len(tmp_metal_ox_storage) - 1)
    output_metal = tmp_metal_ox_storage[index][0]
    output_oxidation_state = tmp_metal_ox_storage[index][1]
    d_electron = int(atom.element(output_metal).group_id) - output_oxidation_state
    multiplicity = spin[d_electron][random.randint(0, len(spins[d_electron]) - 1)]
    print("The output metal is:         " + str(output_metal))
    print("The oxidation state is:      " + str(output_oxidation_state))
    print("The period:                  " + str(atom.element(output_metal).group_id))
    print("The d-electron count is:     " + str(d_electron))
    print("The Multiplicity is:         "+str(multiplicity))
    return output_metal, output_oxidation_state, multiplicity


if __name__ == "__main__":
    for i in range(10):
        total_ligand_charge = input("What is the total charge of your ligands")
        ox_spin(total_ligand_charge, metals_of_interest, spins)
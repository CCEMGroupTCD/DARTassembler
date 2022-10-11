from src.Molecule import RCA_Ligand
# some global coordinates

#
#
# list of metal numbers in pse
metals_in_pse = [el for a in [[21, 31], [39, 49], [57, 81], [89, 113]] for el in range(a[0], a[1])]

#
mini_alphabet = ["a", "b", "c", "d", "e", "f"]


def get_monodentate_list():

    monodentate_ligands = []

    #
    #
    # 1.
    coordinates_ = {0: ["O", [0, 0, 1.4361]], 1: ["H", [0.2096, -0.5615, 2.1227]]}
    Hydroxi = RCA_Ligand(coordinates=coordinates_,
                         denticity=1,
                         ligand_to_metal=[1, 0],
                         name="Hydroxi"
                         )
    monodentate_ligands.append(Hydroxi)

    return monodentate_ligands
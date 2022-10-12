from src.Extracted_Molecule import Extracted_Molecule
import numpy as np
from ase import neighborlist
from scipy.sparse.csgraph import connected_components
from mendeleev import element


class Extracted_molecule_w_charge(Extracted_Molecule):
    """
    We will inherit most of what we did earlier and try to enhance the extraction by some functionality
    As we need the same information (almost) we can run the construction similarly
    """

    def __init__(self, coordinates: dict, csd_code: str):
        super().__init__(coordinates, csd_code)
        self.charge_dict = {}
        self.metal_idx = np.where(self.atomic_numbers == self.original_metal)[0][0]

    def get_reduced_coordinates(self):

        shift_vec = None

        metal_symb = element(int(self.original_metal)).symbol

        for index, coord in self.full_coordinates.items():
            if coord[0] == metal_symb:
                shift_vec = np.array(coord[1])
                break

        if shift_vec is None:
            self.status.append("No shift vector could be estimated")
            return self.full_coordinates

        #
        # create new_coords by shifting
        # and throw out the metal element
        shifted_coordinates, new_index = dict(), 0
        for coord in self.full_coordinates.values():
            if coord[0] != metal_symb:
                new_coord = np.array(coord[1]) - shift_vec
                shifted_coordinates[new_index] = [coord[0], list(new_coord), coord[2]]
                new_index += 1

        return shifted_coordinates

    #
    # everything can stay, but we need to update the extraction method
    def extract_ligands(self, **kwargs):

        new_atomic_num_nns = self.w_o_metal.mol.get_atomic_numbers()[self.metal_neighbor_indices_wo_m]

        assert new_atomic_num_nns.all() == self.orig_atomic_num_nns.all()

        cutoff = neighborlist.natural_cutoffs(self.w_o_metal.mol)
        neighborList = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
        neighborList.update(self.w_o_metal.mol)
        graph = neighborList.get_connectivity_matrix(sparse=False)

        connected_comps = connected_components(graph + np.eye(graph.shape[0]))
        neighbor_conn_comp_list = [connected_comps[1][i] for i in self.metal_neighbor_indices_wo_m]

        denticity_dict = {item: neighbor_conn_comp_list.count(item) for item in
                          neighbor_conn_comp_list}  # component:denticity

        used_index_list = []
        ligand_wo_metal_index_list = self.w_o_metal.mol.get_atomic_numbers()
        coordinates_wo_metal = self.get_reduced_coordinates()

        if not self.status:

            for j, (conn_comp_number, conn_comp_denticity) in enumerate(denticity_dict.items()):
                ligand_index_list = [index for index, atomic_number in
                                     enumerate(ligand_wo_metal_index_list) if connected_comps[1][index]
                                     == conn_comp_number]

                used_index_list += ligand_index_list

                partial_charge = sum([float(coordinates_wo_metal[index][2]) for index in ligand_index_list])

                self.charge_dict[f"Lig_{j}"] = partial_charge

        #
        #
        # now we need to add the metal charge: is written in the full coordinates
        self.charge_dict["Metal"] = float(self.full_coordinates[self.metal_idx][2])

        #
        #
        # Now we need to treat the remainders
        other_atom_charges = [float(coordinates_wo_metal[index][2]) for index in coordinates_wo_metal.keys() if index not in used_index_list]

        if other_atom_charges != []:
            self.charge_dict["Other"] = sum(other_atom_charges)

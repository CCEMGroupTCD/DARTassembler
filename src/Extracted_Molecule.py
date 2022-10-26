import numpy as np
from scipy.spatial.transform import Rotation as R

from ase import io, neighborlist
from scipy.sparse.csgraph import connected_components
from mendeleev import element
from pymatgen.core.periodic_table import Element as Pymatgen_Element

from src.utilities import coordinates_to_xyz_str
from src.constants import metals_in_pse, mini_alphabet
from src.Molecule import RCA_Molecule, RCA_Ligand


def original_xyz_indices_to_indices_wo_metal(orig_coordinates: dict) -> list:
    """
    Converts the indices of atoms of the original xyz files to the new indices when the metal atom is deleted from the xyz. That means all atoms after the metal shift one index up.
    :param orig_coordinates: Coordinates of original xyz file
    :return: Dictionary mapping original xyz indices to new wo_metal indices
    """
    orig_indices = [(idx, l[0]) for idx, l in orig_coordinates.items()]
    
    orig_to_wo_metal_indices = {}
    counter = 0
    for orig_idx, el in orig_indices:
        
        is_metal = Pymatgen_Element(el).is_metal
        if is_metal:
            metal_idx = orig_idx
            continue
        
        orig_to_wo_metal_indices[orig_idx] = counter
        counter += 1
    
    # Double checking.
    assert len(orig_to_wo_metal_indices) == len(orig_indices) - 1
    for orig_idx, new_idx in orig_to_wo_metal_indices.items():
        if orig_idx < metal_idx:
            assert orig_idx == new_idx
        elif orig_idx == metal_idx:
            assert not metal_idx in orig_to_wo_metal_indices
        else:
            assert orig_idx == new_idx + 1
    
    return orig_to_wo_metal_indices


def convert_atomic_props_from_original_xyz_indices_to_indices_wo_metal(atomic_props, orig_coords):
    """
    Changes the indices of an atomic property dictionary in order to match the new indices when the metal is deleted in the xyz.
    :param atomic_props:
    :param orig_coords:
    :return:
    """
    orig_to_wo_metal_indices = original_xyz_indices_to_indices_wo_metal(orig_coords)
    
    atomic_props_new_idc = {}
    for orig_idx, prop in atomic_props.items():
        if orig_idx in orig_to_wo_metal_indices:
            new_idx = orig_to_wo_metal_indices[orig_idx]
            atomic_props_new_idc[new_idx] = prop
    
    all_new_elements = [props[0] for props in atomic_props_new_idc.values()]
    assert not any([Pymatgen_Element(el).is_metal for el in all_new_elements]), 'Found metal in ligand? Or Index Error.'
    return atomic_props_new_idc


class Extracted_Molecule:
    """
    allows us to extract and store information on a molecule from a .xyz from the cambridge database
    in an object of the class "Extracted Molecule"
    """

    def __init__(self, coordinates: dict, csd_code: str, atomic_props: dict={}, global_props: dict={}):

        # read basic properties from xyz_file type class
        #
        with open("../tmp/tmp.xyz", "w+") as text_file:
            text_file.write(coordinates_to_xyz_str(coordinates=coordinates))

        self.complete = RCA_Molecule(mol=io.read("../tmp/tmp.xyz"), atomic_props=atomic_props, global_props=global_props)
        self.csd_code = csd_code
        self.atomic_numbers = self.complete.mol.get_atomic_numbers()
        self.full_coordinates = coordinates
        


        #
        # keeps track of errors
        self.status = list()

        #
        # get atomic number of metal atom
        self.original_metal = self.get_metal_atomic_number()

        #
        # init other properties
        self.ligands = list()
        self.modified_coordinates = None
        self.metal_neighbor_indices_wo_m = None
        self.orig_atomic_num_nns = None
        self.indices_without_metal = None
        self.ligand_coordinates = None
        self.metal_neighbor_indices = None
        self.neighborList = None

        #
        # removing metal
        self.w_o_metal = RCA_Molecule(self.get_mol_wo_metal())

    def get_metal_atomic_number(self):
        # löst den Prozess des metal ding findens etwas eleganter als Michael bisher
        metals = set(self.atomic_numbers).intersection(set(metals_in_pse))
        if len(metals) > 1:
            self.status.append("more than 1 metal atom")
        return metals.pop()
        # Atomic mass of metal

    def get_mol_wo_metal(self):
        new_mol = self.complete.mol.copy()

        metal_idx = np.where(self.atomic_numbers == self.original_metal)[0][0]

        # In-between properties
        cutOff = neighborlist.natural_cutoffs(self.complete.mol)
        #  todo:
        #   Default skin is 0.3
        #   have to modify it a little; justification for 0.2?
        self.neighborList = neighborlist.NeighborList(cutOff, skin=0.2, self_interaction=False, bothways=True)

        self.neighborList.update(self.complete.mol)
        self.metal_neighbor_indices = self.neighborList.get_neighbors(metal_idx)[0]
        self.orig_atomic_num_nns = self.atomic_numbers[self.metal_neighbor_indices]

        # Indizes der an metall gebundenen elemente im molecule ohne metal
        self.metal_neighbor_indices_wo_m = sorted([x if x < metal_idx else x - 1 for x in self.metal_neighbor_indices])

        del new_mol[metal_idx]
        return new_mol

    def modify_coordinates(self):
        """
        rotation and translation to origin
        """
        #
        # get shift vector
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
                shifted_coordinates[new_index] = [coord[0], list(new_coord)]
                new_index += 1

        #
        # Rotation by rotation element
        rotation_element_index = self.metal_neighbor_indices_wo_m[0]
        #
        # (a) get rotation matrix
        a_ = np.array(shifted_coordinates[rotation_element_index][1])
        b = np.array([np.linalg.norm(a_), 0, 0])
        a = a_.reshape(-1, 1).T
        C = R.align_vectors(a, b.reshape(-1, 1).T)
        A = C[0].as_matrix()
        #
        # (b) Rotate coorindates
        new_coordinates = dict()
        for index, coord in shifted_coordinates.items():
            v = np.array(coord[1])
            v_new = v.dot(A)
            new_coords = [el if abs(el) > 0.001 else 0 for el in v_new.tolist()]
            new_coordinates[index] = [coord[0], new_coords]

        return new_coordinates

    def extract_ligands(self, denticity_numbers: list):
        new_atomic_num_nns = self.w_o_metal.mol.get_atomic_numbers()[self.metal_neighbor_indices_wo_m]

        if not new_atomic_num_nns.all() == self.orig_atomic_num_nns.all():
            self.status.append("AssertionError")

        #
        # get graph of molecule without metal to identify connected components,
        # which are the ligands

        cutoff = neighborlist.natural_cutoffs(self.w_o_metal.mol)
        neighborList = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
        neighborList.update(self.w_o_metal.mol)
        graph = neighborList.get_connectivity_matrix(sparse=False)

        connected_comps = connected_components(graph + np.eye(graph.shape[0]))

        #
        # to which componenent the neighbors of the metal belong to
        neighbor_conn_comp_list = [connected_comps[1][i] for i in self.metal_neighbor_indices_wo_m]

        self.denticity_dict = {item: neighbor_conn_comp_list.count(item) for item in
                          neighbor_conn_comp_list}  # component:denticity


        cutoff_denticity = 6
        if sum(self.denticity_dict.values()) > cutoff_denticity:
            self.status.append("Too much bindings evaluated")
        assert cutoff_denticity == max(denticity_numbers) or cutoff_denticity == max(denticity_numbers) + 1, 'Cutoff value should be equal to max(denticity_numbers) or max(denticity_numbers)+1, otherwise it is possible that some ligands of a single complex are saved and some are not, which is harmful for charge calculation.'

        #
        # if status: then there are errors
        if not self.status:
            j = 0               # for name

            self.modified_coordinates = self.modify_coordinates()

            # iterieren über alle möglichen liganden
            for conn_comp_number, conn_comp_denticity in self.denticity_dict.items():
                if conn_comp_denticity in denticity_numbers:

                    #
                    # ligand index list: the indices of the ligands in the xyz_file class
                    ligand_index_list = [index for index, atomic_number in
                                         enumerate(self.w_o_metal.mol.get_atomic_numbers()) if connected_comps[1][index]
                                         == conn_comp_number]
                    #
                    # connections to metal
                    ligand_to_metal = [1 if index in self.metal_neighbor_indices_wo_m else 0
                                       for i, index in enumerate(ligand_index_list)]
                    #
                    #
                    ligand_coordinates = {i: self.modified_coordinates[index] for i, index in enumerate(ligand_index_list)}
                    
                    # Get partial charges of ligand
                    ligand_atomic_props = {}
                    for name, props in self.complete.atomic_props.items():
                        # Correct indices to match indices of ´ligand_coordinates´
                        new_idc_props = convert_atomic_props_from_original_xyz_indices_to_indices_wo_metal(props, self.full_coordinates)
                        new_idx_props = {i: new_idc_props[i] for i in ligand_index_list}
                        ligand_atomic_props[name] = {i: new_idc_props[index] for i, index in enumerate(ligand_index_list)}
                        
                        # Double check if all indices and elements are the same.
                        coord_atoms = {idx: l[0] for idx, l in ligand_coordinates.items()}
                        props_atoms = {idx: l[0] for idx, l in ligand_atomic_props[name].items()}
                        assert coord_atoms == props_atoms, f'Atoms for the atomic property {name} not the same as for the coordinates for molecule {self.csd_code}.'
                    
                    
                    ligand = RCA_Ligand(coordinates=ligand_coordinates,
                                        ligand_to_metal=ligand_to_metal,
                                        original_metal=self.original_metal,
                                        denticity=conn_comp_denticity,
                                        name=f'CSD-{self.csd_code}-0{conn_comp_denticity}-{mini_alphabet[j]}',
                                        csd_code=self.csd_code,
                                        atomic_props=ligand_atomic_props
                                        )

                    j += 1

                    self.ligands.append(ligand)
                
                else:
                    raise Warning(f'One of the ligands with denticity {conn_comp_denticity} of complex {self.csd_code} is not saved due to denticity. This is dangerous (especially for charge calculation) because it can lead to the situation that a complex appears in the database but not all of it\'s ligands are there.')





# if __name__ == '__main__':
#     orig_coordinates_list = [
#                                 {0: ['I', [-0.9, 4.9, 8.4]], 1: ['Cu', [0.8, 5.1, 10.1]], 2: ['S', [1.1, 6.8, 11.6]]},
#                                 {0: ['Cu', [0.8, 5.1, 10.1]], 1: ['I', [-0.9, 4.9, 8.4]], 2: ['S', [1.1, 6.8, 11.6]]},
#                                 {0: ['S', [1.1, 6.8, 11.6]], 1: ['I', [-0.9, 4.9, 8.4]], 2: ['Cu', [0.8, 5.1, 10.1]]}
#                             ]
#
#     atomic_props_list = [
#         {0: ['I', [-0.9]], 1: ['Cu', [0.8]], 2: ['S', [1.1]]},
#         {0: ['Cu', [0.8]], 1: ['I', [-0.9]], 2: ['S', [1.1]]},
#         {0: ['S', [1.1]], 1: ['I', [-0.9]], 2: ['Cu', [0.8]]}
#     ]
#
#     for atomic_props, orig_coords in zip(atomic_props_list, orig_coordinates_list):
#         new_atomic_props = convert_atomic_props_from_original_xyz_indices_to_indices_wo_metal(atomic_props, orig_coords)
    
    
    
import pickle
from copy import deepcopy

import numpy as np
from scipy.spatial.transform import Rotation as R

from ase import io, neighborlist
from scipy.sparse.csgraph import connected_components
from mendeleev import element
from pymatgen.core.periodic_table import Element as Pymatgen_Element
from pymatgen.core.composition import Composition

from src.testing import all_same, all_same_number_of_elements, check_partial_charges_add_up, \
    check_if_metal_entries_consistent
from src.utilities import coordinates_to_xyz_str, convert_atomic_props_from_original_xyz_indices_to_indices_wo_metal, \
    get_all_atomic_properties_in_long_array, get_all_atomic_properties_with_modified_coordinates_wo_metal_in_long_array
from src.constants import metals_in_pse, mini_alphabet
from src.Molecule import RCA_Molecule, RCA_Ligand
from warnings import warn
from src.utilities import atomic_props_dict_to_lists




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
        self.original_metal_symbol = element(int(self.original_metal)).symbol

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
            
            self.n_ligands = len(self.ligands)
            self.complete.global_props['n_ligands'] = len(self.ligands)
            if 'metal_node_degree' in self.complete.global_props:
                self.same_n_ligands_and_metal_node_degree = len(self.ligands) == self.complete.global_props['metal_node_degree']
                self.complete.global_props['same_n_ligands_and_metal_node_degree'] = self.same_n_ligands_and_metal_node_degree

    def assert_that_atomic_property_tuples_of_ligands_are_the_same_as_in_the_complex_wo_metal_with_modified_coordinates(self):
        """
        Checks that the tuple of (element, atomic properties, coordinates) of ligands for one atom is exactly the same as in the original complex, apart from the modification of the coordinates which happened in ´self.extract_ligands´.
        :return: None
        """
        idc, atoms, coords = atomic_props_dict_to_lists(self.full_coordinates)
        all_atomic_props_wo_metal = get_all_atomic_properties_with_modified_coordinates_wo_metal_in_long_array(atoms,
                                                                                                               coords,
                                                                                                               self.complete.atomic_props,
                                                                                                               self.modified_coordinates)
        all_atomic_props_wo_metal = list(map(tuple, all_atomic_props_wo_metal))
        
        all_lig_atomic_props = []
        for lig in self.ligands:
            lig_coord_idc, lig_coord_atoms, lig_coord_values = atomic_props_dict_to_lists(lig.coordinates)
            
            lig_atomic_props = get_all_atomic_properties_in_long_array(lig_coord_atoms,
                                                                           lig_coord_values,
                                                                           lig.atomic_props)
            lig_atomic_props = list(map(tuple, lig_atomic_props))
            all_lig_atomic_props.extend(lig_atomic_props)
            
        # Check for each atom that all xyz and atomic properties are exactly the same as in complex.
        assert sorted(all_atomic_props_wo_metal) == sorted(all_lig_atomic_props), 'Some atomic property tuples of ligands seem to not exist in the complex.'
        
        return
    def run_sanity_checks(self, mol_id) -> bool:
        
        assert mol_id == self.csd_code, f'{mol_id}: Molecular IDs don\'t match: {mol_id} vs {self.csd_code}.'

        # Important big check
        self.assert_that_atomic_property_tuples_of_ligands_are_the_same_as_in_the_complex_wo_metal_with_modified_coordinates()

        if 'partial_charge' in self.complete.atomic_props:
            check_partial_charges_add_up(
                                            mol_id=mol_id,
                                            partial_charges=self.complete.atomic_props['partial_charge'],
                                            total_charge=self.complete.global_props['total_charge']
                                            )
        
        check_if_metal_entries_consistent(
                                            global_props=self.complete.global_props,
                                            mol_id=mol_id,
                                            other_metal_entries=[self.original_metal_symbol, self.original_metal]
                                            )
        
        # all_atoms same
        idc, atoms, coords = atomic_props_dict_to_lists(self.full_coordinates)
        assert all_same_number_of_elements(idc, atoms, coords), f'{mol_id}: Different number of element of indices, atoms or coordinate.'
        
        atomic_numbers_correct = [element(int(num)).symbol for num in self.atomic_numbers] == atoms
        assert atomic_numbers_correct, f'{mol_id}: Atomic numbers not consistent.'
        
        for name, prop in self.complete.atomic_props.items():
            prop_idc, prop_atoms, prop_values = atomic_props_dict_to_lists(prop, flatten=True)
            
            assert prop_atoms == atoms, f'{mol_id}: Inconsistent atoms of property {name}.'
            assert prop_idc == idc, f'{mol_id}: Inconsistent indices of property {name}.'
            assert all_same_number_of_elements(prop_idc, prop_atoms, prop_values), f'{mol_id}: Inconsistent number of element of indices, atoms and values of property {name}.'
        
        assert len(atoms) == self.complete.global_props['num_atoms'], f'{mol_id}: Number of atoms not consistent.'
        
        assert all(self.complete.mol.numbers == self.atomic_numbers)
        assert (self.complete.mol.positions == coords).all()
        assert all(self.complete.mol.pbc == False)
        assert Composition(self.complete.mol.symbols.formula._formula) == Composition(self.complete.global_props['stoichiometry'])
        
        assert all_same_number_of_elements(self.denticity_dict, self.ligands)
        
        
        mod_idc, mod_atoms, mod_coords = atomic_props_dict_to_lists(self.modified_coordinates)
        assert all_same_number_of_elements(mod_idc, mod_atoms, mod_coords)
        assert len(mod_atoms) == len(atoms) - 1
        assert len(set(mod_atoms)) == len(set(atoms)) - 1
        
        # Check if indices and elements of atoms neighboring the metal are correct.
        for idx, atomic_num in zip(self.metal_neighbor_indices_wo_m, self.orig_atomic_num_nns):
            el = element(int(atomic_num)).symbol
            assert self.modified_coordinates[idx][0] == el, f'{mol_id}: The elements of the new coordinates and the elements of atomic neighbors of the metal don\'t match.'
        
        all_lig_atoms = []
        all_lig_p_charges = []
        for lig in self.ligands:
            
            all_prop_idc = []
            all_prop_atoms = []
            for name, prop in lig.atomic_props.items():
                
                lig_idc, lig_atoms, lig_values = atomic_props_dict_to_lists(prop)
                all_prop_idc.append(lig_idc)
                all_prop_atoms.append(lig_atoms)
                if name == 'partial_charge':
                    all_lig_p_charges.extend(lig_values)
            assert all_same(*all_prop_idc)
            assert all_same(*all_prop_atoms)
            
            lig_coord_idc, lig_coord_atoms, lig_coord_values = atomic_props_dict_to_lists(lig.coordinates)
            assert all_same(lig_idc, lig_coord_idc)
            assert all_same(lig_atoms, lig_coord_atoms)
            
            coordinated_elements = [el for el, coordinated in zip(lig_atoms, lig.ligand_to_metal) if coordinated]
            assert len(coordinated_elements) == lig.denticity
            
            
            assert all_same_number_of_elements(lig_atoms, lig.ligand_to_metal, lig.graph.nodes, lig.mol.positions)
            assert all_same(lig_atoms, [element(int(num)).symbol for num in lig.mol.numbers])
            
            correct_ligand_name = lig.name == f'CSD-{self.csd_code}-0{lig.denticity}-{lig.name[-1]}'
            assert correct_ligand_name

            all_lig_atoms.extend(lig_atoms)
        assert sorted(mod_atoms) == sorted(all_lig_atoms), f'{mol_id}: Differing atoms between mod_atoms and all_lig_atoms.'

        atoms_without_metal = deepcopy(atoms)
        atoms_without_metal.remove(lig.original_metal_symbol)
        assert sorted(atoms_without_metal) == sorted(all_lig_atoms), f'{mol_id}: Number or element of atoms in ligands and whole complex without metal don\'t match.'
        
        atoms_not_in_ligands = set(atoms).symmetric_difference(set(all_lig_atoms))
        assert len(atoms_not_in_ligands) == 1, f'{mol_id}: More than one atom difference between whole molecule and all ligands.'
        assert Pymatgen_Element(list(atoms_not_in_ligands)[0]).is_metal, f'{mol_id}: Deleted atom is not metal.'
        assert all_same(lig.original_metal_symbol, element(int(lig.original_metal)).symbol, self.complete.global_props['metal'])
        
        if self.status:
            warn(f'Error status of {self.csd_code}: {self.status}')
            
        return
        



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


if __name__ == '__main__':
    
    ligand_db_file = "../data/LigandDatabases/ligand_db_test.pickle"
    mol_id = 'ADAQUU'

    with open(ligand_db_file, 'rb') as file:
        ligand_db = pickle.load(file)
        print('Loaded ligand db from pickle.')
    
    mol = ligand_db.all_Extracted_Molecules[mol_id]
    mol.run_sanity_checks(mol_id)
    
    print('Done!')
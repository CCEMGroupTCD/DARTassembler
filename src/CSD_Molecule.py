import numpy as np
from scipy.spatial.transform import Rotation as R

from ase import io, Atoms, neighborlist
from scipy.sparse.csgraph import connected_components
from mendeleev import element

from get_csd_xyz_dict import xyz_file
from constants import metals_in_pse
from ASE_Molecule import ASE_Molecule, ASE_Ligand


class CSD_Molecule:
    '''
    allows us to extract and store information on a molecule from a .xyz from the cambridge database
    '''

    def __init__(self, xyz: xyz_file):

        # read basic properties from xyz_file type class
        #
        with open("../tmp/tmp.xyz", "w+") as text_file:
            text_file.write(xyz.get_xyz_file_format_string())

        self.complete = ASE_Molecule(io.read("../tmp/tmp.xyz"))
        self.csd_code = xyz.csd_code
        self.atomic_numbers = self.complete.mol.get_atomic_numbers()
        self.full_coordinates = xyz.coordinates

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

        # removing metal
        self.w_o_metal = ASE_Molecule(self.get_mol_wo_metal())

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
        # todo:
        #  Default skin is 0.3
        #  have to modify it a little; justification for 0.2?
        self.neighborList = neighborlist.NeighborList(cutOff, skin=0.2, self_interaction=False, bothways=True)

        self.neighborList.update(self.complete.mol)
        self.metal_neighbor_indices = self.neighborList.get_neighbors(metal_idx)[0]
        self.orig_atomic_num_nns = self.atomic_numbers[self.metal_neighbor_indices]

        # Indizes der an metall gebundenen elemente im molecule ohne metal
        self.metal_neighbor_indices_wo_m = sorted([x if x < metal_idx else x - 1 for x in self.metal_neighbor_indices])

        del new_mol[metal_idx]
        return new_mol

    def modify_coordinates(self):
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
        # conntected_comps[0] : Anzahl Zusammenhangskomponenten
        # connected_comps[1] : gibt die tatsächlichen komponenten, wobei jedem Knoten ein wert von 0 bis #Zshkomp-1
        #   zugeordnet wird je nach dem zu welcher zshkomp er gehört

        #
        # to which componenent the neighbors of the metal belong to
        neighbor_conn_comp_list = [connected_comps[1][i] for i in self.metal_neighbor_indices_wo_m]

        denticity_dict = {item: neighbor_conn_comp_list.count(item) for item in
                          neighbor_conn_comp_list}  # component:denticity

        if sum(denticity_dict.values()) > 6:
            # die methode klappt erstmal nur für ligands die genau aus octahedral TMCs extrahiert wurden
            # todo: Eigentlich können wir auch größere Denticities verwandeln
            self.status.append("Too much bindings evaluated")

        #
        # if status: then there are errors
        if not self.status:
            j = 1               # for name

            self.modified_coordinates = self.modify_coordinates()

            # iterieren über alle möglichen liganden
            for conn_comp_number, conn_comp_denticity in denticity_dict.items():
                if conn_comp_denticity in denticity_numbers:
                    #
                    # building the new ligand
                    properties = dict()
                    #
                    # naming
                    properties['name'] = f'CSD-{self.csd_code}-0{conn_comp_denticity}-0{j}'
                    j += 1
                    #
                    # denticity
                    properties['denticity'] = conn_comp_denticity
                    #
                    # original metal
                    properties['original_metal'] = self.original_metal
                    #
                    # ligand index list: the indices of the ligands in the xyz_file class
                    ligand_index_list = [index for index, atomic_number in
                                         enumerate(self.w_o_metal.mol.get_atomic_numbers()) if connected_comps[1][index]
                                         == conn_comp_number]
                    #
                    # connections to metal
                    ligand_to_metal = [1 if index in self.metal_neighbor_indices_wo_m else 0
                                       for i, index in enumerate(ligand_index_list)]
                    properties['ligand_to_metal'] = ligand_to_metal
                    #
                    #
                    ligand_xyz = xyz_file(atom_number=len(ligand_index_list),
                                          csd_code=self.csd_code,
                                          coordinates={i: self.modified_coordinates[index] for i, index in enumerate(ligand_index_list)}
                                         )

                    ligand = ASE_Ligand(xyz=ligand_xyz, property_dict=properties)

                    self.ligands.append(ligand)

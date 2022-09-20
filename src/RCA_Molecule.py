import numpy as np
from ase import io, Atoms, neighborlist
from ase.visualize import view
import matplotlib.pyplot as plt
import networkx as nx
from sympy import Point3D, Plane
import collections
from networkx import weisfeiler_lehman_graph_hash as graph_hash


# todo: give that a project based name: RandomComplexAssembler
class RCA_Molecule:
    """
    The basis for the molecule objects we handle throughout the process.
    Added by some basic functionality
    """

    def __init__(self, mol: Atoms):
        """
        an ase version of the molecules we have
        """
        self.mol = mol

    def get_graph(self):
        """
        get the graph representing the molecule as a networkx graph
        """
        cutOff = neighborlist.natural_cutoffs(self.mol)
        neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
        neighborList.update(self.mol)

        adjacency_matrix = neighborList.get_connectivity_matrix(sparse=False)
        labels = {i: el for i, el in enumerate(self.mol.get_chemical_symbols())}

        gr = nx.Graph(adjacency_matrix)

        for node in gr.nodes():
            gr.nodes[node]["label"] = labels[node]

        return gr

    def view_graph(self):
        cutOff = neighborlist.natural_cutoffs(self.mol)
        neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
        neighborList.update(self.mol)
        graph = neighborList.get_connectivity_matrix(sparse=False)
        rows, cols = np.where(graph == 1)
        edges = zip(rows.tolist(), cols.tolist())
        gr = nx.Graph()
        gr.add_edges_from(edges)
        labels = {i: el for i, el in enumerate(self.mol.get_chemical_symbols())}
        nx.draw_networkx(gr, node_size=500, with_labels=True, labels=labels)
        plt.show()

    def view_3d(self):
        view(self.mol)


class xyz_file:

    def __init__(self, atom_number: int, coordinates: dict, csd_code: str):
        self.total_atom_number = atom_number
        self.coordinates = coordinates
        # {number in xyz_file:[Name von Atom, Coordinaten als np.array]
        self.csd_code = csd_code

    def get_xyz_file_format_string(self):
        str_ = f"{self.total_atom_number}\n \n"
        for coord in self.coordinates.values():
            str_ += f"{coord[0]} \t {coord[1][0]} \t {coord[1][1]} \t {coord[1][2]} \n"

        return str_

    def xyz_to_ASE(self):
        with open("../tmp/tmp.xyz", "w+") as text_file:
            text_file.write(self.get_xyz_file_format_string())

        return RCA_Molecule(io.read("../tmp/tmp.xyz"))


class RCA_Ligand(RCA_Molecule):
    '''
        inherits its methods from ASE Molecule, while needs some extra information
        '''

    def __init__(self, xyz: xyz_file, **kwargs):
        super().__init__(Atoms([coord[0] for coord in xyz.coordinates.values()],
                               positions=[coord[1] for coord in xyz.coordinates.values()]))

        self.xyz = xyz
        self.csd_code = xyz.csd_code

        if "denticity" in kwargs.keys():
            self.denticity = kwargs['denticity']

        if "ligand_to_metal" in kwargs.keys():
            self.ligand_to_metal = kwargs["ligand_to_metal"]

        if "name" in kwargs.keys():
            self.name = kwargs["name"]

        if 'original_metal' in kwargs.keys():
            self.original_metal = kwargs['original_metal']

        # todo: Braucht noch ein bisschen testing
        self.type = self.get_type()

    def get_assembly_dict(self):
        """
            only to get the attributes required for the assembly with cians script a little faster
            :return: {index: list, type: list, xyz_str: str}
            """
        dict_ = {"index": [i for i, el in enumerate(self.ligand_to_metal) if el == 1]}

        dict_["type"] = [self.xyz.coordinates[i][0] for i in dict_["index"]]
        dict_["str"] = self.xyz.get_xyz_file_format_string()

        return dict_

    def get_type(self):
        if self.denticity == 2 or self.denticity == 1:
            return 'p'
        elif self.denticity == 5:
            return 'np'
        else:
            if self.check_if_planar() is True:
                return 'p'
            else:
                return 'np'

    def check_if_planar(self, eps=1):
        """
            :param eps: durch try'n'error obtained
            eps für (d=4) -> 1
            :return:
            """
        fc = list()

        for index, information in self.xyz.coordinates.items():
            if self.ligand_to_metal[index] == 1:
                fc.append(information[1])

        assert len(fc) == self.denticity

        if self.denticity == 3:
            c1, c2, c3 = Point3D(fc[0]), Point3D(fc[1]), Point3D(fc[2])
            E = Plane(c1, c2, 0)
            if E.distance(c3) < eps:
                return True

        if self.denticity == 4:
            c1, c2, c3, c4 = Point3D(fc[0]), Point3D(fc[1]), Point3D(fc[2]), Point3D(fc[3])
            E = Plane(c1, c2, c3)
            if E.distance(c4) < eps:
                return True

        return False

    def same_sum_formula(self, other):
        sum_formula_1 = [a[0] for a in self.xyz.coordinates.values()]
        sum_formula_2 = [a[0] for a in other.xyz.coordinates.values()]

        return collections.Counter(sum_formula_1) == collections.Counter(sum_formula_2)

    @staticmethod
    def node_check(dict1, dict2):
        return dict1["label"] == dict2["label"]

    def __eq__(self, other):
        if not self.same_sum_formula(other):
            return False

        if nx.is_isomorphic(self.get_graph(), other.get_graph(), node_match=self.node_check):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(graph_hash(self.get_graph()))

    def add_atom(self, symbol: str, coordinates: list[float]):
        if len(coordinates) != 3:
            print("Wrong number of coordinates specified")
            raise ValueError

        self.xyz.coordinates[len(self.xyz.coordinates)] = [symbol, coordinates]

        self.mol = self.xyz.xyz_to_ASE()

    def remove_last_element_in_xyz(self):
        del self.xyz.coordinates[max(self.xyz.coordinates, key=self.xyz.coordinates.get)]
        self.mol = self.xyz.xyz_to_ASE()


from ase import io, Atoms, neighborlist
from ase.visualize import view
import matplotlib.pyplot as plt
import networkx as nx
from sympy import Point3D, Plane
import collections
from networkx import weisfeiler_lehman_graph_hash as graph_hash
from mendeleev import element
import numpy as np

from src.utilities import atomic_props_dict_to_lists


# Package name: RandomComplexAssembler (RCA)
class RCA_Molecule:
    """
    The basis for the molecule objects we handle throughout the process.
    Added by some basic functionality
    """

    def __init__(self, mol: Atoms, atomic_props: dict={}, global_props: dict={}):
        """
        an ase version of the molecules we have
        """
        self.mol = mol
        self.atomic_props = atomic_props
        self.global_props = global_props
        
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
        self.get_graph_hash()
        self.hash = hash(self.graph_hash)
        
        return self.hash
    
    def get_graph_hash(self):
        if self.has_graph_hash():
            return self.graph_hash
        
        self.get_graph()
        self.graph_hash = graph_hash(self.graph, iterations=3, digest_size=16)
        
        return self.graph_hash
    
    def same_sum_formula(self, other):
        sum_formula_1 = [a[0] for a in self.coordinates.values()]
        sum_formula_2 = [a[0] for a in other.coordinates.values()]

        return collections.Counter(sum_formula_1) == collections.Counter(sum_formula_2)
    
    @staticmethod
    def node_check(dict1, dict2):
        return dict1["label"] == dict2["label"]
    
    def get_adjacency_matrix(self):
        cutOff = neighborlist.natural_cutoffs(self.mol)
        neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
        neighborList.update(self.mol)

        return neighborList.get_connectivity_matrix(sparse=False)
    
    def has_graph(self):
        return hasattr(self, 'graph')
    
    def has_graph_hash(self):
        return hasattr(self, 'graph_hash')
    def get_graph(self):
        """
        get the graph representing the molecule as a networkx graph
        returns a labeled, undirected graph
        """
        if self.has_graph():
            return self.graph
        
        adjacency_matrix = self.get_adjacency_matrix()
        labels = {i: el for i, el in enumerate(self.mol.get_chemical_symbols())}

        self.graph = nx.Graph(adjacency_matrix)

        for node in self.graph.nodes():
            self.graph.nodes[node]["label"] = labels[node]
        
        return self.graph
    

    def view_graph(self):
        """
        simple plot of the molecule as a graph, only connectivity, no distances
        """
        gr = self.get_graph()
        labels = {i: el for i, el in enumerate(self.mol.get_chemical_symbols())}
        nx.draw_networkx(gr, node_size=500, with_labels=True, labels=labels)
        plt.show()

    def view_3d(self):
        """
        easy employment of the ase functionality
        """
        view(self.mol)
    
    def run_sanity_tests(self, mol_id):
        # Check partial charge summed up is charge
        
        
        # Check if metal and metal_str_if_exists in global properties are consistent

        pass


class RCA_Ligand(RCA_Molecule):
    """
    Ligands need a little more information, hence we need to expand the RCA_Molecule class
    """

    def __init__(self, coordinates: dict, denticity: int, ligand_to_metal: list, atomic_props: dict={}, global_props: dict={}, **kwargs):
        super().__init__(Atoms([coord[0] for coord in coordinates.values()],
                               positions=[coord[1] for coord in coordinates.values()]),
                         atomic_props=atomic_props,
                         global_props=global_props)

        self.coordinates = coordinates
        self.denticity = denticity
        self.ligand_to_metal = ligand_to_metal

        if "csd_code" in kwargs.keys():
            self.csd_code = kwargs['csd_code']

        if "name" in kwargs.keys():
            self.name = kwargs["name"]

        if 'original_metal' in kwargs.keys():
            self.original_metal = kwargs['original_metal']
            self.original_metal_symbol = element(int(self.original_metal)).symbol

    def planar_check(self, eps3=2, eps4=1):
        """
            :param eps: durch try'n'error obtained
            eps für (d=4) -> 1
            :return:
            """
        fc = list()

        for index, information in self.coordinates.items():
            if self.ligand_to_metal[index] == 1:
                fc.append(information[1])

        assert len(fc) == self.denticity

        if self.denticity == 3:
            c1, c2, c3 = Point3D(fc[0]), Point3D(fc[1]), Point3D(fc[2])
            E = Plane(c1, c2, Point3D([0, 0, 0]))
            if round(E.distance(c3)) < eps3:
                return True

        if self.denticity == 4:
            c1, c2, c3, c4 = Point3D(fc[0]), Point3D(fc[1]), Point3D(fc[2]), Point3D(fc[3])
            E = Plane(c1, c2, c3)
            if round(E.distance(c4)) < eps4:
                return True

        return False

    def get_total_atom_number(self):
        return len(list(self.coordinates.keys()))

    def get_xyz_file_format_string(self):
        """
        returns a string that can be written into an .xyz file
        """
        str_ = f"{self.get_total_atom_number()}\n \n"
        for coord in self.coordinates.values():
            str_ += f"{coord[0]} \t {coord[1][0]} \t {coord[1][1]} \t {coord[1][2]} \n"

        return str_

    def print_to_xyz(self, path: str):
        if not path.endswith(".xyz"):
            print("No valid filename")
            raise NameError
        with open(path, "w+") as file:
            file.write(self.get_xyz_file_format_string())

    def get_assembly_dict(self):
        """
        only to get the attributes required for the assembly with cians script a little faster
        :return: {index: list, type: list, xyz_str: str}
        """
        dict_ = {"index": [i for i, el in enumerate(self.ligand_to_metal) if el == 1]}

        dict_["type"] = [self.coordinates[i][0] for i in dict_["index"]]
        dict_["str"] = self.get_xyz_file_format_string()

        return dict_

    def mol_to_asemol(self):
        self.print_to_xyz(path="../tmp/tmp.xyz")
        return io.read("../tmp/tmp.xyz")

    def add_atom(self, symbol: str, coordinates: list[float]):
        if len(coordinates) != 3:
            print("Wrong number of coordinates specified")
            raise ValueError

        self.coordinates[len(self.coordinates)] = [symbol, coordinates]
        self.mol = self.mol_to_asemol()

    def remove_last_element_in_xyz(self):
        del self.coordinates[max(list(self.coordinates.keys()))]
        self.mol = self.mol_to_asemol()

    # todo: unittests!!
    def NO_check(self):

        return set(self.get_assembly_dict()["type"]).issubset({"N", "O"})

    def betaH_check(self):
        """
        returns True if beta Hydrogen are contained, false if not
        """
        A = self.get_adjacency_matrix()

        B = np.matmul(A, A)
        # The second power of the adjacency matrix, i.e. A^2[i,j] represents the number of paths of length two
        # from i to j. Hence, as we are only interested in hydrogens that have distance two to our functional atoms
        # we can make quick use of that

        for functional_index in self.get_assembly_dict()["index"]:
            for index in self.coordinates.keys():

                if B[functional_index, index] > 0 and self.coordinates[index][0] == "H":
                    # print("Beta Hydrogen detected")
                    # self.view_3d()
                    return True

        return False

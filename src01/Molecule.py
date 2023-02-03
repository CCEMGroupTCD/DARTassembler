# standard Pyhton packages
import hashlib
import warnings
import networkx as nx
import random
import numpy as np
from copy import deepcopy

# some special functions which are required
from pymatgen.core.periodic_table import Element as Pymatgen_Element
from ase.visualize import view
from networkx import weisfeiler_lehman_graph_hash as graph_hash
from scipy.spatial.transform import Rotation as R
from sympy import Point3D, Plane

# collection of molecule objects of other packages
from ase import io, Atoms, neighborlist
from pymatgen.core.structure import Molecule as PyMatMol
from molSimplify.Classes.mol3D import mol3D

# importing own scripts
from src01.utilities_graph import graph_from_graph_dict, graph_to_dict_with_node_labels, view_graph, graphs_are_equal, unify_graph, get_sorted_atoms_and_indices_from_graph, get_reindexed_graph, find_node_in_graph_by_label
from src01.utilities import identify_metal_in_ase_mol, make_None_to_NaN, update_dict_with_warning_inplace
from src01.utilities_Molecule import get_standardized_stoichiometry_from_atoms_list
from src03_Assembly.stk_utils import RCA_Mol_to_stkBB, convert_RCA_to_stk_Molecule


# Package name: RandomComplexAssembler (RCA)
class RCA_Molecule:
    """
    The idea of this class is to build a method which contains an ase molecule for visualization
    but also other convenient features,
    as a graph representation
    and all the global information we have at hand when creating a database
    """

    def __init__(self, mol: Atoms = None,
                 atomic_props: dict = None,
                 global_props: dict = None,
                 graph=None,
                 graph_creating_strategy="default",
                 has_ligands=True,
                 reindex_graph: bool = False,
                 other_props: dict={},
                 **kwargs
                 ):
        """
        :param graph_creating_strategy: If we dont give a graph we might want to specify the graph creating strategy
        :param kwargs: additional parameters which are specified for the graph creation, such as
            skin: float
            cutoff_corrections_for_metals
            strategy: (pymatgen strategy)

        """

        if atomic_props is None:
            atomic_props = {}
        if global_props is None:
            global_props = {}

        self.atomic_props = atomic_props
        self.global_props = global_props
        # Generate mol from atomic_props if possible and no mol given yet
        self.mol = self.get_mol_from_input(mol)

        self.add_additional_molecule_information_to_global_props()

        if has_ligands is True:
            # if we expect ligands, we can set up an empty ligand list
            self.ligands = []

        if graph is None:
            # we pass over all the kwargs, the respective functions will only grab those they need
            # but no kwarg is mandatory as all have default values set
            self.graph = self.make_graph(graph_creating_strategy=graph_creating_strategy,
                                         **kwargs
                                         )
        else:
            if reindex_graph:
                graph = get_reindexed_graph(graph)
            self.graph = graph

        # As the graphs are now not optional anymore we can also make the hashes baseline
        self.graph_hash = self.get_graph_hash()
        self.hash = self.get_hash()

        self.stoichiometry = self.get_standardized_stoichiometry()

        # Set kwargs so that they become properties of the molecule
        self.set_other_props_as_properties(other_props=other_props)

    def get_mol_from_input(self, mol):
        if mol is None:
            coord_list_3D = [[self.atomic_props[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
                             enumerate(self.atomic_props["x"])]
            atom_list = self.atomic_props["atoms"]
            mol = Atoms(atom_list, positions=coord_list_3D)

        return mol

    def set_other_props_as_properties(self, other_props):
        for prop, value in other_props.items():
            differing_prop_already_exists = hasattr(self, prop) and self.__getattribute__(prop) != value
            if differing_prop_already_exists:
                raise Warning(f'Property {prop} is tried to be set but already exists.')
            self.__setattr__(prop, value)

        return

    def add_additional_molecule_information_to_global_props(self):
        n_atoms = len(self.atomic_props['atoms'])
        mol_weight = self.mol.get_masses().sum()
        info = {
                    'n_atoms': n_atoms,
                    'molecular_weight': mol_weight
                    }

        update_dict_with_warning_inplace(
                                            dict_to_update=self.global_props,
                                            dict_with_information=info
                                        )

        return info

    @classmethod
    def make_from_atomic_properties(cls,
                                    atomic_props_mol: dict,
                                    global_props_mol: dict,
                                    graph=None,
                                    graph_creating_strategy: str = "default",
                                    **kwargs
                                    ):
        """
        A more convenient creation method, as in general the atomic properties already imply the information for
        the ase mol (and the pymatgen mol)
        """

        coord_list_3D = [[atomic_props_mol[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
                         enumerate(atomic_props_mol["x"])]
        atom_list = atomic_props_mol["atoms"]

        return cls(mol=Atoms(atom_list, positions=coord_list_3D),
                   atomic_props=atomic_props_mol,
                   global_props=global_props_mol,
                   graph=graph,
                   graph_creating_strategy=graph_creating_strategy,
                   **kwargs
                   )

    # basic view 3D function
    def view_3d(self):
        """
        easy employment of the ase functionality
        """
        view(self.mol)

    # Graph stuff
    def make_graph(self, graph_creating_strategy: str, **kwargs):
        """
        This method also allows to overwrite graphs via the command
        self.graph = self.make_graph(...)
        which will set self.graph to the newly created graph
        """
        from src01.GraphCreation import GraphCreation

        return GraphCreation(
            graph_creating_strategy=graph_creating_strategy,
            molecule=self.mol,
            atomic_props=self.atomic_props,
            **kwargs
        ).G

    def view_graph(self, node_label='node_label', node_size=150):
        """
        simple plot of the molecule as a graph, only connectivity, no distances
        """
        view_graph(self.graph, node_label=node_label, node_size=node_size)

    # Now the hashing
    def has_graph_hash(self):
        return False if self.graph_hash is None else True

    def has_hash(self):
        return False if self.hash is None else True

    def __hash__(self):
        if not self.has_hash():
            self.get_hash()
        return self.hash

    def get_graph_hash(self):

        self.graph_hash = graph_hash(self.graph, node_attr='node_label', iterations=3, digest_size=16)
        self.get_hash()
        return self.graph_hash

    def get_hash(self):
        if not self.has_graph_hash():
            self.get_graph_hash()

        # this hash is not randomized, unlike pythons inbuilt hash()
        self.hash = int(hashlib.md5(self.graph_hash.encode(encoding='UTF-8', errors='strict')).hexdigest(), 16)
        return self.hash

    def get_standardized_stoichiometry(self) -> str:
        """
        Returns a string with the stoichiometry in a standardized way. We use the Hill notation, except that for elements with stoichiometry 1 this 1 is written as well.
        :return: stoichiometry (str)
        """
        formula = get_standardized_stoichiometry_from_atoms_list(self.atomic_props['atoms'])
        return formula

    def __eq__(self, other):
        if not self.stoichiometry == other.stoichiometry:
            return False

        return graphs_are_equal(self.graph, other.graph)

    def __ne__(self, other):
        return not self.__eq__(other)

    # Finally we want to be able to turn our complex into an .xyz file
    def get_xyz_file_format_string(self):
        """
        returns a string that can be written into an .xyz file
        """
        str_ = f"{len(self.atomic_props['x'])}\n \n"
        for i, _ in enumerate(self.atomic_props['x']):
            str_ += f"{self.atomic_props['atoms'][i]}  {self.atomic_props['x'][i]}  {self.atomic_props['y'][i]}  {self.atomic_props['z'][i]} \n"

        return str_

    def print_to_xyz(self, path: str):
        if not path.endswith(".xyz"):
            print("No valid filename")
            raise NameError
        with open(path, "w+") as file:
            file.write(self.get_xyz_file_format_string())

    # helper method for the de-assembly
    def ligand_naming(self, denticity: int, ligand_list) -> (str, str):

        if "CSD_code" in self.global_props.keys():
            lig_key = f'CSD-{self.global_props["CSD_code"]}'
            csd = self.global_props["CSD_code"]
        else:
            lig_key = "NoCSD"
            csd = None
        from src01.constants import mini_alphabet
        j = 0
        while True:
            ligand_name = f'{lig_key}-0{denticity}-{mini_alphabet[j]}'
            if ligand_name not in [lig.name for lig in ligand_list]:
                break
            else:
                j += 1

        return ligand_name, csd

    def shift_metal_to_origin(self):
        """
        Actually, this method is outdated not required.
        However, the idea was to shift the molecule to the origin, so that the metal is right in the origin
        which is the first part - the shift

        and afterwards, rotate the the metal, so that one of the functional atoms is algined with the x-axis
        """
        #
        # get shift vector
        metal_symb = identify_metal_in_ase_mol(self.mol)
        metal_index = self.atomic_props["atoms"].index(metal_symb)

        # do the shifting
        for key in ["x", "y", "z"]:
            self.atomic_props[key] = [value - self.atomic_props[key][metal_index] for value in self.atomic_props[key]]

        """
        # Now we rotate. First we identify the element we want to rotate on
        metal_node = find_node_in_graph_by_label(G=self.graph, label_to_find=metal_symb, expected_hits=1)
        rotation_element_index = list(self.graph.neighbors(metal_node))[0]

        # Next, we obtain the rotation matrix
        rotation_element_vector = np.array([self.atomic_props[key][rotation_element_index] for key in ["x", "y", "z"]])
            # now we can idenftify the desired vector, on which we'd like to rotate the vector
        desired_rotation = np.array([np.linalg.norm(rotation_element_vector), 0, 0])

        # and we can thus find the rotation matrix
        A = R.align_vectors(rotation_element_vector.reshape(-1, 1).T, desired_rotation.reshape(-1, 1).T)[0].as_matrix()

        for index, _ in enumerate(self.atomic_props["x"]):
            location_vec_for_index = np.array([self.atomic_props[key][index] for key in ["x", "y", "z"]])

            # now we rotate
            v_new = location_vec_for_index.dot(A)

            # and we modify the atomic properties
            for i, key in enumerate(["x", "y", "z"]):
                self.atomic_props[key][index] = v_new.tolist()[i] if abs(v_new.tolist()[i]) > 0.001 else 0
        """
        return

    def check_input_inherit_global_properties(self, inherit_global_properties: list) -> list:
        """
        Checks whether `inherit_global_properties` has correct input format.
        :param inherit_global_properties: input
        :return: correctly specified inherit_global_properties
        """
        if inherit_global_properties is None:
            inherit_global_properties = list(self.global_props.keys())
        else:
            unknown_global_property = [prop for prop in inherit_global_properties if not prop in self.global_props]
            if unknown_global_property:
                raise ValueError(
                    f'Unknown values {unknown_global_property}. All properties in inherit_global_properties must be found in self.global_properties.')

        return inherit_global_properties

    def de_assemble(self, Testing: bool = False, inherit_global_properties: list = ['CSD_code']):
        """
        now only graph based, makeslife waaay easier

        All that is based on the following assumption:
        the nodes are denoted as integers and we shall assume that these integers
        correspond to their index in the atomic properties, i.e. the position
        so, for example for node 1, the x-coordinate can be found by
        self.atomic_props["x"][1] or more general self.atomic_props["x"][node]
        Ich glaube das haelt auf jeden Fall fuer alle selber erzeugten Graphen

        den Dummen Testparamater muss ich mitschleppen, um das gegen meine alte, ultra behinderte Grapherstellung testen
        zu koennen
        """
        if not hasattr(self, "ligands"):
            self.ligands = []

        inherit_global_properties = self.check_input_inherit_global_properties(inherit_global_properties)

        # wie gesagt, etwas outdated eigentlich
        self.shift_metal_to_origin()

        atoms, idc = get_sorted_atoms_and_indices_from_graph(self.graph)
        if 'atoms' in self.atomic_props:
            assert atoms == self.atomic_props['atoms'], 'Order of atoms in graph and in atomic_props doesn\'t match.'

        # first we gather some information around the metal in the initial graph
        graph = deepcopy(self.graph)
        metal_in_complex = identify_metal_in_ase_mol(self.mol)
        metal_node = find_node_in_graph_by_label(G=graph, label_to_find=metal_in_complex, expected_hits=1)
        metal_neighbors = list(graph.neighbors(metal_node))  # sind jetzt alle benachbarten nodes vom metal
        metal_neighbor_elements = [graph.nodes[i]['node_label'] for i in metal_neighbors]

        for i, el in zip(idc, atoms):
            graph.nodes[i]['orig_idx'] = i
            graph.nodes[i]['metal_neighbor'] = i in metal_neighbors
            assert graph.nodes[i]['node_label'] == el, 'atom and index don\'t match.'

        # next we create the ripped graph
        ripped_graph = deepcopy(graph)
        ripped_graph.remove_node(metal_node)

        conn_components = [sorted(comp) for comp in
                           nx.connected_components(ripped_graph)]  # sorting of components very important
        conn_components = [comp for comp in
                           sorted(conn_components, key=str)]  # important: sort by string of list, makes order of ligands unique

        for component in conn_components:
            # if this set is empty, the ligand has no connection to the metal
            functional_atom_indices = sorted(list(set(component).intersection(set(metal_neighbors))))
            assert max(component) <= len(
                ripped_graph), 'Indices dont make sense. Most likely this is an implementation error where due to the deletion of the metal atom the indices of original and ripped graph dont match.'

            denticity = len(functional_atom_indices)
            if denticity == 0:
                # because in the assembly 0 is the placeholder for the reactant, whereas -1 means this is just a remainder in the .xyz, not
                # connected to the metal at all
                denticity = -1

            # with that it is insanely easy to determine the atomic properties of the ligand
            assert component == sorted(
                component), 'The list of ligand indices is not sorted, but that is assumed in many parts of this project.'
            ligand_atomic_props = {key_: [el for i, el in enumerate(item_) if i in component] for key_, item_ in
                                   self.atomic_props.items()}
            ligand_atomic_props['original_complex_indices'] = component

            # problem: the functional_atom_indices are the indices in the full original metal, rather than the ligand only
            # so we have to convert them to the index in the ligand_atomic_props

            # dont have to compute the ligand graph new, removes computation time and
            # was a source for errors before
            ligand_graph = deepcopy(graph.subgraph(component))
            atoms_lig, idc_lig = get_sorted_atoms_and_indices_from_graph(ligand_graph)

            local_indices = [component.index(ind) for ind in functional_atom_indices]
            local_elements = [ligand_atomic_props['atoms'][i] for i in local_indices]

            # Doublechecking
            if local_indices != []:   # otherwise error for unconnected ligands
                assert max(local_indices) < len(ligand_graph), 'local_indices make no sense, an index is greater than the number of elements.'
            assert all(el in metal_neighbor_elements for el in local_elements), 'Inconsistent elements of the metal neighbors.'
            assert local_indices == sorted(local_indices), 'local_indices is not sorted but should be.'
            assert atoms_lig == ligand_atomic_props['atoms'], 'elements of graph and atomic_props not consistent'
            assert [atoms[i] for i in component] == ligand_atomic_props['atoms'], 'ligand atoms not consistent with original complex atoms.'

            ligand_name, csd = self.ligand_naming(denticity, self.ligands)

            kwargs = {'skin_': 0.3} if Testing == True else {'csd_code': csd}
            ligand_global_props = {prop: self.global_props[prop] for prop in inherit_global_properties}
            new_lig = RCA_Ligand(denticity=denticity,
                                 ligand_to_metal=local_indices,
                                 atomic_props=ligand_atomic_props,
                                 name=ligand_name,
                                 graph=ligand_graph,
                                 global_props=ligand_global_props,
                                 original_metal=Pymatgen_Element(metal_in_complex).Z,
                                 **kwargs
                                 )

            self.ligands.append(new_lig)

        if self.ligands == []:
            print(f'WARNING: Complex {self.global_props["CSD_code"]} has no ligands extracted.')

    def write_to_mol_dict(self):
        """
        Graph hash not that important for molecules I guess
        """
        return {
            "atomic_props": self.atomic_props,
            "global_props": self.global_props,
            "graph_dict": graph_to_dict_with_node_labels(self.graph),
            # "graph_hash": self.graph_hash
        }

    @classmethod
    def read_from_mol_dict(cls,
                           dict_: dict,
                           graph_creating_strategy: str = "default",
                           **kwargs
                           ):
        """
        generate an RCA mol object from a mol_dict,
        we can either take a graph or a graph dict as input for the graph so far
        :param graph_creating_strategy: which graph strategy to follow in order to create a graph
        :param kwargs: additional requirements for the graph creation, we can dump basically anything
        """
        assert {"atomic_props", "global_props", "graph_dict"}.issubset(set(dict_.keys()))
        if dict_["graph_dict"] is None:
            return cls.make_from_atomic_properties(atomic_props_mol=dict_["atomic_props"],
                                                   global_props_mol=dict_["global_props"],
                                                   graph_creating_strategy=graph_creating_strategy,
                                                   **kwargs
                                                   )
        elif isinstance(dict_["graph_dict"], dict):
            return cls.make_from_atomic_properties(atomic_props_mol=dict_["atomic_props"],
                                                   global_props_mol=dict_["global_props"],
                                                   graph=graph_from_graph_dict(dict_["graph_dict"])
                                                   )
        elif isinstance(dict_["graph_dict"], nx.MultiGraph) or isinstance(dict_["graph_dict"], nx.Graph):
            return cls.make_from_atomic_properties(atomic_props_mol=dict_["atomic_props"],
                                                   global_props_mol=dict_["global_props"],
                                                   graph=unify_graph(G=dict_["graph_dict"])
                                                   )
        else:
            print("Unreadable graph format")
            return cls.make_from_atomic_properties(atomic_props_mol=dict_["atomic_props"],
                                                   global_props_mol=dict_["global_props"],
                                                   graph_creating_strategy=graph_creating_strategy,
                                                   **kwargs
                                                   )

    #
    #
    # converting into other classes of mol objects
    def to_molsimplify_mol(self):
        """
        Here we convert the RCA Molecule into a molsimplify mol3D object
        """

        new_mol = mol3D()
        new_mol.readfromstring(xyzstring=self.get_xyz_file_format_string())

        return new_mol

    def to_stk_mol(self):
        return convert_RCA_to_stk_Molecule(self)

    def to_pymatMol(self):
        return PyMatMol(species=self.atomic_props["atoms"],
                        coords=[[self.atomic_props[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
                                enumerate(self.atomic_props["x"])])


class RCA_Ligand(RCA_Molecule):
    """
    Ligands need a little more information, hence we need to expand the RCA_Molecule class
    """

    def __init__(self,
                 denticity: int,
                 ligand_to_metal: list,
                 atomic_props: dict,
                 name: str,
                 graph=None,
                 global_props: dict = None,
                 other_props={},
                 **kwargs
                 ):
        """

        :param graph: If None, it will be generated custom. Maybe only important for monodentates and reactant
        """
        coord_list_3D = [[atomic_props[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
                         enumerate(atomic_props["x"])]
        atom_list = atomic_props["atoms"]

        super().__init__(mol=Atoms(atom_list, positions=coord_list_3D),
                         atomic_props=atomic_props,
                         global_props=global_props,
                         graph=graph,
                         has_ligands=False,
                         other_props=other_props,
                         **kwargs
                         )

        self.coordinates = {i: [at, coord_list_3D[i]] for i, at in enumerate(atom_list)}

        self.denticity = denticity
        self.name = name

        # the indices and elements where the ligands was bound to the metal
        self.ligand_to_metal = ligand_to_metal
        self.local_elements = self.get_local_elements()

        if "csd_code" in kwargs.keys():
            self.csd_code = kwargs['csd_code']

        if "unique_name" in kwargs:
            self.unique_name = kwargs["unique_name"]

        if 'original_metal' in kwargs.keys():
            self.original_metal = kwargs['original_metal']

        try:
            self.original_metal_symbol = Pymatgen_Element.from_Z(self.original_metal).symbol
        except (ValueError, AttributeError):
            pass

    def get_local_elements(self) -> list:
        """
        Calculates the elements in the first coordination sphere from `self.ligand_to_metal`.
        :return: list of elements in first coordination sphere
        """
        return [self.atomic_props["atoms"][i] for i in self.ligand_to_metal]

    def get_assembly_dict(self):
        """
        only to get the attributes required for the assembly with cians script a little faster
        :return: {index: list, type: list, xyz_str: str}
        """
        dict_ = {}
        dict_["index"] = self.ligand_to_metal
        dict_["type"] = [self.atomic_props["atoms"][i] for i in self.ligand_to_metal]
        dict_["str"] = self.get_xyz_file_format_string()

        return dict_

    def add_atom(self, symbol: str, coordinates: list[float]):
        if len(coordinates) != 3:
            print("Wrong number of coordinates specified")
            raise ValueError

        self.atomic_props["atoms"].append(symbol)
        self.atomic_props["x"].append(coordinates[0])
        self.atomic_props["y"].append(coordinates[1])
        self.atomic_props["z"].append(coordinates[2])

        # now we need to update the self.mol, graph wont be updated
        self.print_to_xyz(path="../tmp/tmp.xyz")
        self.mol = io.read("../tmp/tmp.xyz")

    def functional_atom_check(self, atoms_of_interest: [str, list[str]] = None):
        if atoms_of_interest is None:
            atoms_of_interest = ["N", "O"]
        elif isinstance(atoms_of_interest, str):
            atoms_of_interest = [atoms_of_interest]
        return set(self.get_assembly_dict()["type"]).issubset(set(atoms_of_interest))

    def betaH_check(self):
        """
        returns True if beta Hydrogen are contained, false if not
        """
        A = nx.adjacency_matrix(self.graph).todense()

        B = np.matmul(A, A)
        # The second power of the adjacency matrix, i.e. A^2[i,j] represents the number of paths of length two
        # from i to j. Hence, as we are only interested in hydrogens that have distance two to our functional atoms
        # we can make quick use of that

        # todo: Der functional_index ist nicht der Index im Graph!!!
        #   glaube das waere geloest wenn wir den Graph ordern und zwar kanonisch nach label des nodes

        for functional_index in self.ligand_to_metal:
            for index, atom_symbol in enumerate(self.atomic_props['atoms']):

                if B[functional_index, index] > 0 and atom_symbol == "H":
                    # print("Beta Hydrogen detected")
                    # self.view_3d()
                    return True

        return False

    # todo: unittests!!
    def planar_check(self, eps=2):
        """
        :param eps: durch try'n'error obtained
        eps für (d=4) -> 1
        :return:
        """
        functional_coords = [[self.atomic_props[key][i] for key in ["x", "y", "z"]] for i in self.ligand_to_metal]

        assert len(functional_coords) == self.denticity, f"Error in Planar Check for ligand {self.name}"

        if self.denticity == 3:
            c1, c2, c3 = Point3D(functional_coords[0]), Point3D(functional_coords[1]), Point3D(functional_coords[2])
            E = Plane(c1, c2, Point3D([0, 0, 0]))
            if round(E.distance(c3)) < eps:
                return True

        if self.denticity == 4:
            c1, c2, c3, c4 = Point3D(functional_coords[0]), Point3D(functional_coords[1]), Point3D(
                functional_coords[2]), Point3D(functional_coords[3])
            E = Plane(c1, c2, c3)
            if round(E.distance(c4)) < eps:
                return True

        return False

    def write_to_mol_dict(self):
        return {
            "stoichiometry": self.stoichiometry,
            "atomic_props": self.atomic_props,
            "global_props": self.global_props,
            "graph_dict": graph_to_dict_with_node_labels(self.graph),
            "denticity": self.denticity,
            "ligand_to_metal": self.ligand_to_metal,
            "local_elements": self.local_elements,
            "name": self.name,
            "CSD_code": self.csd_code if hasattr(self, "csd_code") else None,
            "original_metal": self.original_metal if hasattr(self, "original_metal") else None,
            "original_metal_symbol": self.original_metal_symbol if hasattr(self, "original_metal_symbol") else None,
            "graph_hash": self.graph_hash
        }

    @classmethod
    def read_from_mol_dict(cls, dict_: dict, **kwargs):
        """
        Reads the ligand from a provided dictionary.
        """
        necessary_props = ["atomic_props", "global_props", "graph_dict", "denticity", "name", 'ligand_to_metal']
        assert set(necessary_props).issubset(set(dict_.keys())), f'Any of the necessary keys {necessary_props} is not present.'

        other_props = {key: val for key, val in dict_.items() if not key in necessary_props}

        # Optionally add graph if it is present in the dictionary
        kwargs = {'graph_dict': graph_from_graph_dict(dict_['graph_dict'])} if not (dict_['graph_dict'] is None) else {}

        return cls(
            atomic_props=dict_["atomic_props"],
            global_props=dict_["global_props"],
            denticity=dict_["denticity"],
            name=dict_["name"],
            ligand_to_metal=dict_['ligand_to_metal'],
            other_props=other_props,
            **kwargs
        )

    # some stk functionality
    def to_stk_bb(self):
        """
        this is really only designed for ligands as a normal RCA_Molecule doesnt have the required properties
        """
        return RCA_Mol_to_stkBB(self)

class RCA_Complex(RCA_Molecule):

    def __init__(self,
                 mol: Atoms = None,
                 atomic_props: dict = None,
                 global_props: dict = None,
                 pymat_mol=None,
                 graph=None,
                 graph_creating_strategy="default",
                 has_ligands=True,
                 reindex_graph: bool = False,
                 other_props={},
                 **kwargs):

            super().__init__(
                            mol=mol,
                            atomic_props=atomic_props,
                            global_props=global_props,
                            pymat_mol=pymat_mol,
                            graph=graph,
                            graph_creating_strategy=graph_creating_strategy,
                            has_ligands=has_ligands,
                            reindex_graph=reindex_graph,
                            other_props=other_props,
                            **kwargs
                             )
            self.metal = identify_metal_in_ase_mol(self.mol)
            self.metal_oxi_state = make_None_to_NaN(self.global_props['metal_oxi_state'])
            self.charge = make_None_to_NaN(self.global_props['charge'])
            self.metal_idx = self.atomic_props['atoms'].index(self.metal)

            self.add_additional_complex_information_to_global_props()

            return

    def add_additional_complex_information_to_global_props(self):
        info = {}
        if 'partial_charge' in self.atomic_props:
            metal_q = self.atomic_props['partial_charge'][self.atomic_props['atoms'].index(self.metal)]
            info['metal_partial_charge'] = metal_q

        update_dict_with_warning_inplace(
                                            dict_to_update=self.global_props,
                                            dict_with_information=info
                                        )

        return info


    def write_to_mol_dict(self):
        """
        Outputs complex as dictionary.
        """
        return {
            "atomic_props": self.atomic_props,
            "global_props": self.global_props,
            "graph_dict": graph_to_dict_with_node_labels(self.graph),
            'mol_id': self.global_props['CSD_code'],
            'stoichiometry': self.stoichiometry,
            'metal': self.metal,
            'metal_oxi_state': self.metal_oxi_state,
            'total_q': self.charge,
            'ligands': [lig.write_to_mol_dict() for lig in self.ligands],
            'graph_hash': self.graph_hash
        }

    @classmethod
    def read_from_mol_dict(cls, dict_: dict, **kwargs):
        """
        Reads the ligand from a provided dictionary.
        """
        necessary_props = ["atomic_props", "global_props", "graph_dict"]
        assert set(necessary_props).issubset(set(dict_.keys())), f'Any of the necessary keys {necessary_props} is not present.'

        if 'ligands' in dict_:
            dict_['ligands'] = [RCA_Ligand.read_from_mol_dict(lig) for lig in dict_['ligands']]
            # This sounds stupid but this is because otherwise the RCA_molecule class sets up self.ligands = [] which then collides when the actual ligands which are read in here are added because for safety the code checks that it doesn't overwrite anything.
            has_ligands = False
        else:
            has_ligands = True

        if 'total_q' in dict_:
            dict_['charge'] = dict_['total_q']
            del dict_['total_q']

        other_props = {key: val for key, val in dict_.items() if not key in necessary_props}

        # Optionally add graph if it is present in the dictionary
        if 'graph_dict' in dict_ and not (dict_['graph_dict'] is None):
            graph = graph_from_graph_dict(dict_['graph_dict'])
        else:
            graph = None

        return cls(
            atomic_props=dict_["atomic_props"],
            global_props=dict_["global_props"],
            graph=graph,
            has_ligands=has_ligands,
            other_props=other_props,
            **kwargs
        )

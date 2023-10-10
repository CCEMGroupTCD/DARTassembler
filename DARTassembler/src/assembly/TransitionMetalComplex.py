from datetime import date
import stk
from copy import deepcopy
import networkx as nx
import hashlib
import numpy as np
from pathlib import Path
from typing import Union, List, Tuple
import json
import re

from DARTassembler.src.assembly.stk_utils import stkBB_to_networkx_graph
from DARTassembler.src.constants.Periodic_Table import DART_Element as element

from DARTassembler.src.ligand_extraction.utilities_Molecule import original_metal_ligand
from DARTassembler.src.ligand_extraction.Molecule import RCA_Molecule, RCA_Ligand
from DARTassembler.src.ligand_extraction.utilities_graph import graphs_are_equal, \
    get_sorted_atoms_and_indices_from_graph, view_graph
from DARTassembler.src.assembly.utilities_assembly import generate_pronounceable_word

atomic_number_Hg = 80


# todo:
# Auf dem complex fehlen noch ein paar Methoden:

# 3. eine die ein Gaussian erstellt.


class TransitionMetalComplex:

    mol: RCA_Molecule

    def __init__(self,
                 atomic_props: dict,
                 graph: nx.Graph,
                 metal_oxi_state: int,
                 metal_idx: int,
                 charge: int,
                 ligand_props: dict,
                 spin: int = None
                 ):
        """
        :param atomic_props: dict
        :param graph: nx.Graph
        :param metal_oxi_state: int
        :param metal_idx: int
        :param charge: int
        :param ligand_props: dict
        :param spin: int
        """
        self.atomic_props = atomic_props
        self.graph = graph
        self.ligand_props = ligand_props
        self.metal_oxi_state = metal_oxi_state
        self.metal_idx = metal_idx
        self.charge = charge
        self.spin = spin

        self.metal = self.atomic_props["atoms"][self.metal_idx]
        self.total_charge = self.charge  # deprecated, use self.charge instead

        assert sorted(self.graph.nodes) == list(range(len(self.graph.nodes))), f"The graphs indices are not in order: {list(self.graph.nodes)}"
        graph_elements, indices = atoms, _ = get_sorted_atoms_and_indices_from_graph(self.graph)
        assert graph_elements == self.atomic_props["atoms"]
        assert nx.is_connected(self.graph), "The graph is not fully connected!"

        self.mol = RCA_Molecule.make_from_atomic_properties(
                                                            atomic_props_mol=self.atomic_props,
                                                            global_props_mol={},
                                                            graph=self.graph
                                                            )

        self.functional_groups = {key: lig['donor_elements'] for key, lig in ligand_props.items()}
        self.donor_elements = [el for elements in self.functional_groups.values() for el in elements]

        self.mol_id = self.create_random_name()





    @staticmethod
    def get_total_charge(metal_charge_, ligands_):

        charge = metal_charge_

        for ligand in ligands_.values():
            try:
                charge += ligand.charge
            except AttributeError:
                try:
                    charge += ligand.pred_charge
                except AttributeError:
                    return None

        return charge

    @staticmethod
    def assemble_name(metal, ligands):
        """
        this methods encodes our naming scheme
        """
        name = metal
        for i, ligand in ligands.items():
            try:
                name += f"_{ligand.name}"
            except AttributeError:
                # Ligand has no name assigned
                name += f"_{i}dentLig"

        return name

    @staticmethod
    def stk_Constructed_Mol_to_atomic_props(compl: stk.ConstructedMolecule) -> dict:
        atomic_props = {
            "x": [],
            "y": [],
            "z": [],
            "atoms": []
        }

        # to this end we first obtain the indices of interest
        # indices_for_non_Hg = [i for i, atom in enumerate(compl._atoms) if atom.get_atomic_number() != atomic_number_Hg]
        indices = list(range(len(compl._atoms)))

        # Now we can extract the types of atoms
        atomic_props["atoms"] = [element(compl._atoms[i].get_atomic_number()).symbol for i in indices]

        for (x, y, z) in compl.get_atomic_positions(indices):
            atomic_props["x"].append(x)
            atomic_props["y"].append(y)
            atomic_props["z"].append(z)

        return atomic_props

    @staticmethod
    def merge_graph_from_ligands(ligands, metal) -> Tuple[nx.Graph, List, List]:
        """
        Merges the graphs from the ligands into one graph. The metal is added as a node with index 0 and connected to the donor atoms of the ligands.
        :param ligands: dict[RCA_Ligand]
        :param metal: str
        :return: Tuple of the merged graph of the complex, the indices of the ligand atoms and the indices of the ligand donor atoms
        """
        ligand_indices = []
        ligand_donor_indices = []
        old_graphs = [deepcopy(lig.graph) for lig in ligands.values()]

        # relabel the nodes of the old graphs, otherwise merging won't work
        i = 1   # start at 1 because 0 is the metal
        for old_graph in old_graphs:
            node_mapping = {node: i + k for k, node in enumerate(sorted(old_graph.nodes))}
            nx.relabel_nodes(old_graph, mapping=node_mapping, copy=False)
            ligand_indices.append(list(node_mapping.values()))
            i += len(old_graph.nodes)

        # now we create the new graph by merging everything
        graph = nx.Graph()
        graph.add_nodes_from([(0, {"node_label": metal}),])     # add metal node
        for H in old_graphs:
            graph.add_nodes_from(H.nodes(data=True))            # add ligand nodes
            graph.add_edges_from(H.edges())                     # add ligand edges

        # Add bonds to metal
        for old_graph, lig in zip(old_graphs, ligands.values()):
            ligand_donor_indices.append([])
            for i in lig.ligand_to_metal:
                assert lig.atomic_props['atoms'][i] in lig.local_elements, f"Atom {lig.atomic_props['atoms'][i]} is not a donor atom of ligand {lig.name}!"
                donor_idx = sorted(old_graph.nodes)[i]  # The index in ligand_to_metal is the index of the donor in the atomic_properties, so taking the ith node of the sorted graph node gives the donor index.
                graph.add_edge(0, donor_idx)
                ligand_donor_indices[-1].append(donor_idx)

        # Check if everything is valid
        all_donor_elements = [el for lig in ligands.values() for el in lig.local_elements]
        coordinated_elements = [graph.nodes[node]['node_label'] for node in nx.all_neighbors(graph, 0)]
        assert sorted(all_donor_elements) == sorted(coordinated_elements), f"Coordinated elements {coordinated_elements} do not match donor elements {all_donor_elements}!"
        assert nx.is_connected(graph), "The graph is not fully connected!"
        assert all([set(ligand_donor_indices[i]).issubset(set(ligand_indices[i])) for i in range(len(ligand_indices))]), "The ligand donor indices are not subset of the ligand indices!"
        assert sorted(graph.nodes) == list(range(len(graph.nodes))), f"The graphs indices are not in order: {list(graph.nodes)}"
        assert sorted(nx.all_neighbors(graph, 0)) == sorted(i for lig in ligand_donor_indices for i in lig), "The metal is not connected to all donor atoms!"

        return graph, ligand_indices, ligand_donor_indices

    def create_random_name(self, length=8, decimals=6):
        """
        Generate a hash name of the molecule based on its xyz coordinates and elements. If coordinates or elements change, the name will change.
        """
        xyz = self.mol.get_xyz_as_array()
        sorted_indices = np.lexsort((xyz[:, 2], xyz[:, 1], xyz[:, 0]), axis=0)
        xyz = np.round(xyz, decimals=decimals)  # round to 6 decimals to get rid of numerical noise
        xyz = xyz[sorted_indices]
        elements = [el for _, el in sorted(zip(sorted_indices, self.mol.get_elements_list()))] # sort elements according to xyz
        hash_string = str(elements) + str(xyz)  # make hash string

        # Generate a pronounceable word from the hash
        name = generate_pronounceable_word(length=length, seed=hash_string)

        return name

    def to_data_dict(self):
        props = deepcopy(self.__dict__)
        complex_properties = {}

        mol = props["mol"].write_to_mol_dict()
        complex_properties['atomic_props'] = mol['atomic_props']
        complex_properties['global_props'] = mol['global_props']
        complex_properties['graph_dict'] = mol['graph_dict']

        complex_properties['ligand_props'] = props['ligand_props']
        complex_properties['metal_oxi_state'] = props['metal_oxi_state']
        complex_properties['spin'] = props['spin']
        complex_properties['metal'] = props['metal']
        complex_properties['charge'] = props['charge']
        complex_properties['mol_id'] = props['mol_id']

        return complex_properties

    def to_json(self,
                path: Union[Path, str]
                ):
        complex_properties = self.to_data_dict()
        with open(path, "w") as file:
            json.dump(complex_properties, file)

    @classmethod
    def from_stkBB(cls,
                 compl: stk.ConstructedMolecule = None,
                 ligands: dict[RCA_Ligand] = None,
                 metal_charge: int = 2,
                 metal: str = "Fe",
                 metal_idx: int = 0,
                 spin: int = 0
                 ):
        """
        :param compl: stk.ConstructedMolecule
        :param ligands: dict[RCA_Ligand]
        :param metal: str
        :param spin: int
        :param metal_charge: int
        """
        atomic_props = cls.stk_Constructed_Mol_to_atomic_props(compl)
        graph, ligand_indices, ligand_donor_indices = cls.merge_graph_from_ligands(ligands, metal)
        charge = cls.get_total_charge(metal_charge, ligands)

        ligand_props = {
            key: {
                "unique_name": ligand.unique_name,
                'ligand_indices_in_complex': ligand_indices[key],
                'donor_indices': ligand_donor_indices[key],
                'donor_elements': ligand.local_elements,
                'donor_bond_lengths': ligand.stats['coordinating_atom_distances_to_metal'],
                "stoichiometry": ligand.stoichiometry,
                'denticity': ligand.denticity,
                'pred_charge': ligand.pred_charge,
                'pred_charge_is_confident': ligand.pred_charge_is_confident,
                'graph_hash_with_metal': ligand.graph_hash_with_metal,
                'has_good_bond_orders': ligand.has_good_bond_orders,
                'warnings': ligand.warnings,
                'occurrences': ligand.occurrences,
                "n_atoms": ligand.n_atoms,
            }
            for key, ligand in ligands.items()
        }


        complex = cls(
            atomic_props=atomic_props,
            graph=graph,
            ligand_props=ligand_props,
            metal_oxi_state=metal_charge,
            metal_idx=metal_idx,
            charge=charge,
            spin=spin
        )

        return complex

    @classmethod
    def from_json(cls,
                  path: Union[Path, str]
                  ):
        """
        Read TransitionMetalComplex from json file and return it.
        """

        a = cls()
        # vielleicht mit: in der init: if complex is None: pass -> leeres objekt und dann hier fuellen

        with open(path, "r") as file:
            properties = json.load(file)

        if not {'ligand_props', 'total_charge', 'mol', 'name', 'functional_groups',
                'mol_id'}.issubset(set(properties.keys())):
            raise KeyError("Missing keys in input json for TMC generation")

        a.mol = RCA_Molecule.read_from_mol_dict(properties["mol"])
        a.graph = a.mol.graph

        del properties["mol"]

        for k, v in properties.items():
            a.__setattr__(k, v)

        return a

    def get_com_format_string(self,
                              basis_set_dict: dict,
                              cluster_path: str = "/home/michael/molsimp_comfiles/Co_31a_14_OH/Co_31a_14_OH.com"
                              ):

        header_ = """%chk=Co_31a_14_OH_LSb3lyp.chk\n%nprocshared=40\n%mem=100GB\n#p guess=read  gen scrf=(smd, solvent=h2o) pseudo=read scf=xqc ub3lyp pop=(regular, npa)"""

        path_line = f"{cluster_path} auto generated "

        coordinate_part = self.mol.get_xyz_file_format_string().split("\n \n")[1]

        try:
            basis_set_part = "\n****\n".join([f"-{atom_symbol} 0\n{basis_set_dict[atom_symbol]}" for atom_symbol in set(self.mol.atomic_props["atoms"])])
        except:
            basis_set_part = ""
        try:
            metal_instructions_ = basis_set_dict[self.metal].split('\n')[0]
            final_part = f"{self.metal} \n{metal_instructions_}"
        except:
            final_part = ""

        return f"{header_}\n{path_line}\n{coordinate_part}\n{basis_set_part}\n{final_part}"

    def to_com(self,
               path: Union[Path, str],
               basis_set_dict: dict,
               cluster_path: str = "/home/michael/molsimp_comfiles/Co_31a_14_OH/Co_31a_14_OH.com"
               ):

        with open(path, "w") as file:
            file.write(self.get_com_format_string(
                basis_set_dict=basis_set_dict,
                cluster_path=cluster_path
            ))

    def to_gaussian_string(self,
                           filename: str,
                           num_processors: int,
                           memory: int,
                           charge: int,
                           multiplicity: int,
                           metal_basis_set: str,
                           output_directory: str):

        basis_set_dict = {"H": "6-31g(d,p)",

                          "C": "6-31g(d)",
                          "N": "6-31+g(d)",
                          "O": "6-31+g(d)",

                          "P": "6-31+g(d)",
                          "S": "6-31+g(d)",
                          "As": "6-31+g(d)",
                          "Se": "6-31+g(d)",

                          "F": "6-31+g(d)",
                          "Cl": "6-31+g(d)",
                          "Br": "6-31+g(d)",
                          "I":  "6-31+g(d)",

                          "other": "6-31g(d)"}

        header = f"""%chk={filename}.chk
%nprocshared={num_processors}
%mem={memory}GB
#p opt rwb97xd/gen pseudo=read \n
continue calc\n
{charge} {multiplicity}\n"""


        coordinates = self.mol.get_xyz_file_format_string().split("\n \n")[1]

        metal_basis_set = f"""-{self.metal} 0
{metal_basis_set}
F 1 1.0
1.050 1.0
****\n"""
        full_atom_str = ""
        for ligand in self.ligand_props.values():
            for character in ligand["stoichiometry"]:
                if character.isnumeric():
                    pass
                else:
                    full_atom_str = full_atom_str + character

            pass
        full_atom_list = re.split('(?<=.)(?=[A-Z])', full_atom_str)
        reduced_atom_list = list(set(full_atom_list))

        basis_set_string = ""
        for atom in reduced_atom_list:
            print("atoms")
            try:
                basis_set_string = basis_set_string + f"""-{atom} 0
{basis_set_dict[atom]}
****\n"""
            except:
                basis_set_string = basis_set_string + f"""-{atom} 0
{basis_set_dict["other"]}
****\n"""

        pre_link = """Au
lanl2dz\n\n"""
        link = f"""--Link1--
%chk={filename}.chk
#p Geom=AllCheck pseudo=read guess=read rwb97xd/gen pop=nbo7read\n
"""

        final_lines = """\nAu
lanl2dz\n
$nbo aonbo=c $end\n
"""
        gaussian_string = header+coordinates+"\n"+metal_basis_set+basis_set_string+"\n"+pre_link+link+metal_basis_set+basis_set_string+final_lines
        return gaussian_string

    def to_gaussian_string_Frank(self,
                           filename: str,
                           num_processors: int,
                           memory: int,
                           charge: int,
                           multiplicity: int,
                           metal_basis_set: str,
                           output_directory: str):

        basis_set_dict = {"H": "6-31G*",

                          "C": "6-31G*",
                          "N": "6-31G*",

                          "F": "6-31G*",

                          "other": "6-31G*"}

        header = f"""%chk={filename}.chk
%nprocshared={num_processors}
%mem={memory}GB
#p opt freq b3lyp/gen scrf=(smd,solvent=n,n-dimethylformamide) nosymm
pop=(NBO,full,CM5) cphf=conver=7 empiricaldispersion=gd3
int=(acc2e=11,grid=ultrafine) pseudo=cards\n
Title Card Required\n
{charge} {multiplicity}\n"""

        coordinates = self.mol.get_xyz_file_format_string().split("\n \n")[1]

        metal_basis_set = f"""{self.metal} 0
{metal_basis_set}
****\n"""
        full_atom_str = ""
        for ligand in self.ligand_props.values():
            for character in ligand["stoichiometry"]:
                if character.isnumeric():
                    pass
                else:
                    full_atom_str = full_atom_str + character

            pass
        full_atom_list = re.split('(?<=.)(?=[A-Z])', full_atom_str)
        reduced_atom_list = list(set(full_atom_list))

        basis_set_string = ""
        for atom in reduced_atom_list:
            print("atoms")
            try:
                basis_set_string = basis_set_string + f"""{atom} 0
{basis_set_dict[atom]}
****\n"""
            except:
                basis_set_string = basis_set_string + f"""{atom} 0
{basis_set_dict["other"]}
****\n"""

        pre_link = """Au
lanl2dz\n\n"""
        link = f"""--Link1--
%chk={filename}.chk
#p Geom=AllCheck pseudo=read guess=read rwb97xd/gen pop=nbo7read\n
"""

        final_lines = """Cu 0
SDD\n
"""
        gaussian_string = header + coordinates + "\n" + metal_basis_set + basis_set_string + "\n" +  final_lines
        return gaussian_string



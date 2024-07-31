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
import tempfile

from DARTassembler.src.assembly.stk_utils import stkBB_to_networkx_graph
from DARTassembler.src.constants.Periodic_Table import DART_Element as element
from DARTassembler.src.ligand_extraction.io_custom import read_xyz
from DARTassembler.src.ligand_extraction.utilities import angle_between_ab_ac_vectors

from DARTassembler.src.ligand_extraction.utilities_Molecule import original_metal_ligand
from DARTassembler.src.ligand_extraction.Molecule import RCA_Molecule, RCA_Ligand
from DARTassembler.src.ligand_extraction.utilities_graph import graphs_are_equal, \
    get_sorted_atoms_and_indices_from_graph, view_graph, graph_from_graph_dict, get_graph_hash
from DARTassembler.src.assembly.utilities_assembly import generate_pronounceable_word

atomic_number_Hg = 80




class TransitionMetalComplex(object):
    """
    This class represents a mono-metallic transition metal complex assembled by DART. It contains all kind of information about the complex, including the ligands, the metal, the graph, the atomic properties, etc. Importantly, the graph and the atomic properties have consistent indices. A TransitionMetalComplex can be saved as a json file and loaded again.
    """

    mol: RCA_Molecule

    def __init__(self,
                 atomic_props: dict,
                 graph: nx.Graph,
                 metal_oxi_state: int,
                 metal_idx: int,
                 charge: int,
                 ligand_props: dict,
                 ):
        """
        :param atomic_props: dict
        :param graph: nx.Graph
        :param metal_oxi_state: int
        :param metal_idx: int
        :param charge: int
        :param ligand_props: dict
        """
        self.atomic_props = atomic_props
        self.graph = graph
        self.ligand_props = ligand_props
        self.metal_oxi_state = metal_oxi_state
        self.metal_idx = metal_idx
        self.charge = charge
        self.graph_hash = get_graph_hash(self.graph)

        self.metal = self.atomic_props["atoms"][self.metal_idx]
        self.metal_position = [self.atomic_props['x'][self.metal_idx], self.atomic_props['y'][self.metal_idx], self.atomic_props['z'][self.metal_idx]]
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

        self.donor_indices = sorted(self.graph.neighbors(self.metal_idx))
        self.donor_elements = [self.atomic_props['atoms'][idx] for idx in self.donor_indices]
        self.donor_positions = [[self.atomic_props['x'][idx], self.atomic_props['y'][idx], self.atomic_props['z'][idx]] for idx in self.donor_indices]

        self.functional_groups = {key: lig['donor_elements'] for key, lig in ligand_props.items()}
        assert sorted(self.donor_elements) == sorted([el for elements in self.functional_groups.values() for el in elements]), f"The donor indices {self.donor_indices} do not match the donor elements {self.donor_elements}!"

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

    def create_random_name(self, length=8, decimals=6, from_graph=False):
        """
        Generate a hash name of the molecule based on its xyz coordinates and elements. If coordinates or elements change, the name will change.
        """
        if from_graph:
            hash_string = self.graph_hash
        else:
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
        complex_properties['metal'] = props['metal']
        complex_properties['metal_idx'] = props['metal_idx']
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
                 ):
        """
        :param compl: stk.ConstructedMolecule
        :param ligands: dict[RCA_Ligand]
        :param metal: str
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
        )

        return complex

    @classmethod
    def from_dir(cls,
                 path: Union[Path, str],
                 xyz: Union[Path, str, tuple] = None
                 ):
        """
        Read TransitionMetalComplex from a directory containing the `*_data.json` file.
        """
        name = Path(path).name
        json_path = Path(path, f"{name}_data.json")
        return cls.from_json(json_path, xyz=xyz)

    @classmethod
    def from_json(cls,
                    path: Union[Path, str],
                    xyz: Union[Path, str, tuple] = None
                  ):
        """
        Read TransitionMetalComplex from json file and return it.
        :param path: Union[Path, str]
        :param xyz: Union[Path, str] = None. If given, the xyz file will be used for the coordinates instead of the one in the json file.
        """
        try:
            with open(path, "r") as file:
                properties = json.load(file)
        except:
            raise ValueError(f"Could not read json file {path}!")

        try:
            is_old_example_json = 'mol' in properties['complex'].keys()
        except:
            is_old_example_json = False

        if is_old_example_json:
            init_dict = cls.get_init_dict_from_old_json_which_works_for_the_Pd_Ni_example(properties)
        else:
            init_dict = cls.get_init_dict_from_json(properties)

        if xyz is not None:
            try:
                elements, coords = xyz
            except:
                xyz = Path(xyz).resolve()
                if not xyz.exists():
                    raise ValueError(f"The provided xyz file {xyz} does not exist!")

                elements, coords = read_xyz(xyz)
            same_elements = elements == init_dict['atomic_props']['atoms']
            if not same_elements:
                raise ValueError(f"The elements in the xyz file {xyz} do not match the elements in the json file {path}!")

            # Replace coordinates
            init_dict['atomic_props']['x'] = coords[:, 0].tolist()
            init_dict['atomic_props']['y'] = coords[:, 1].tolist()
            init_dict['atomic_props']['z'] = coords[:, 2].tolist()

        return cls(**init_dict)

    @staticmethod
    def get_init_dict_from_json(properties: dict) -> dict:
        comp = properties['complex']

        metal_idx = [idx for idx, el in enumerate(comp['atomic_props']['atoms']) if el == comp['metal']]
        assert len(metal_idx) == 1, f"Expected 1 metal centre, found {len(metal_idx)}."
        metal_idx = metal_idx[0]


        return dict(
                atomic_props=comp['atomic_props'],
                graph=graph_from_graph_dict(comp['graph_dict']),
                metal_oxi_state=comp['metal_oxi_state'],
                metal_idx=metal_idx,
                charge=comp['charge'],
                ligand_props=comp['ligand_props'],
                )

    @staticmethod
    def get_init_dict_from_old_json_which_works_for_the_Pd_Ni_example(properties: dict) -> dict:
        """
        Read TransitionMetalComplex from json file and return it. This is the old version which works for jsons in the old format for the Pd_Ni cross coupling example.
        """
        is_example_json = [str(prop['Metal']) for prop in properties['assembly_input_settings']['Batches']] == ["{'element': 'Pd', 'oxidation_state': 2, 'spin': 1}", "{'element': 'Ni', 'oxidation_state': 2, 'spin': 1}"]
        assert is_example_json, "This method is only for the old example jsons!"

        comp = properties['complex']
        mol = comp['mol']
        charge = comp['total_charge']
        metal_idx = [idx for idx, el in enumerate(mol['atomic_props']['atoms']) if el in ['Pd', 'Ni']][0]
        metal_os = 2  # hard coded because this is only for the old example run where this always was 2 for both Pd and Ni
        ligand_props = comp['ligand_props']
        for idx, lig in ligand_props.items():
            lig['donor_elements'] = comp['functional_groups'][idx]
        return dict(
            atomic_props=mol['atomic_props'],
            graph=graph_from_graph_dict(mol['graph_dict']),
            metal_oxi_state=metal_os,
            metal_idx=metal_idx,
            charge=charge,
            ligand_props=ligand_props,
            )



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

    def get_donor_indices_from_indices_or_elements(self, atoms: List[Union[str,int]]) -> Union[int, List[int]]:
        """
        Convert atom symbols to indices if unique, otherwise raise error.
        :param atoms: List of atom symbols or indices
        :return: List of atom indices or single atom index if only one atom was provided
        """
        if isinstance(atoms, str) or isinstance(atoms, int):
            list_provided = False
            atoms = [atoms]
        else:
            list_provided = True

        atom_indices = []
        for idx, atom in enumerate(atoms):
            if isinstance(atom, int):
                atom_indices.append(atom)
            elif isinstance(atom, str):
                if atom not in self.donor_elements:
                    raise ValueError(f"Element {atom} is not in donor elements {self.donor_elements}!")
                elif self.donor_elements.count(atom) > 1:
                    raise ValueError(f"Element {atom} is not unique in donor elements {self.donor_elements}!")
                else:
                    idx = [idx for idx, el in zip(self.donor_indices, self.donor_elements) if el == atom][0]
                    atom_indices.append(idx)
            else:
                raise ValueError(f"Expected atom index or symbol, got {type(atom)}!")

        if not list_provided:
            atom_indices = atom_indices[0]

        return atom_indices

    def get_bite_angle(self, atoms: Union[int, str, List[Union[str, int]]]) -> float:
        """
        Calculates the bite angle between two coordinating atoms with the metal atom as the vertex.
        """
        if not len(atoms) == 2:
            raise ValueError(f"Expected 2 atom indices, got {len(atoms)}!")

        atom_indices = self.get_donor_indices_from_indices_or_elements(atoms)

        if not all(idx in self.donor_indices for idx in atom_indices):
            raise ValueError(f"Atom indices {atom_indices} are not in the donor indices {self.donor_indices}!")

        positions = [self.metal_position] + [[self.atomic_props['x'][idx], self.atomic_props['y'][idx], self.atomic_props['z'][idx]] for idx in atom_indices]

        # Calculate the angle between the three points, where the returned angle is the one at the first point (the metal)
        angle = angle_between_ab_ac_vectors(*positions, degrees=True)

        return angle

    def get_donor_metal_bond_length(self, atom: Union[str, int]) -> float:
        idx = self.get_donor_indices_from_indices_or_elements(atom)
        pos1 = self.metal_position
        pos2 = [self.atomic_props['x'][idx], self.atomic_props['y'][idx], self.atomic_props['z'][idx]]

        return np.linalg.norm(np.array(pos1) - np.array(pos2))

    def get_xtb_descriptors(self, n_unpaired: int = None) -> dict:
        from dev.src11_machine_learning.utils.utilities_ML import get_xtb_descriptors
        elements, coordinates = self.mol.get_elements_list(), self.mol.get_xyz_as_array()
        desc = get_xtb_descriptors(xyz=(elements, coordinates), charge=self.charge, n_unpaired=n_unpaired)

        return desc


if __name__ == '__main__':

    # old_complex_path = '/Users/timosommer/PhD/projects/RCA/projects/DART/examples/Pd_Ni_Cross_Coupling/dev/output/data_before_restarts/DART_Example_Pd_Ni_Complexes/batches/P_N_Donors_Ni_Metal_Centre/complexes/ABADEZIX_PN_Ni/ABADEZIX_PN_Ni_data.json'
    # new_complex_path = '/Users/timosommer/PhD/projects/RCA/projects/DART/testing/integration_tests/assembly/data_output/batches/Integration_test_1/complexes/AZEPOBOB/AZEPOBOB_data.json'
    # old_complex = TransitionMetalComplex.from_json(old_complex_path)
    # new_complex = TransitionMetalComplex.from_json(new_complex_path)
    # bite_angle2 = old_complex.get_bite_angle(['P', 'N'])

    xtb_comp_path = '/Users/timosommer/PhD/projects/RCA/projects/DART/examples/Pd_Ni_Cross_Coupling/dev/xtb_calculations/relaxations_new_8e-4/complexes/ORIYEBOD_PN_Pd/ORIYEBOD_PN_Pd_data.json'
    xtb_complex = TransitionMetalComplex.from_json(xtb_comp_path)
    xtb = xtb_complex.get_xtb_descriptors()




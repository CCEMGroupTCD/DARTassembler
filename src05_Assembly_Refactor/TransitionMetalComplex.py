from datetime import date
import stk
from mendeleev import element
from copy import deepcopy
import networkx as nx
import hashlib
import pandas as pd
from pathlib import Path
from typing import Union
import json
import re
from constants.constants import project_path

from src01.utilities_Molecule import original_metal_ligand
from src01.Molecule import RCA_Molecule, RCA_Ligand
from src01.utilities_graph import graphs_are_equal

atomic_number_Hg = 80


# todo:
# Auf dem complex fehlen noch ein paar Methoden:

# 3. eine die ein Gaussian erstellt.


class TransitionMetalComplex:
    """
    Here we define how we want the final object to look like
    to bring us in a good position for the post assembly filtering
    """
    mol: RCA_Molecule

    def __init__(self,
                 compl: stk.ConstructedMolecule = None,
                 ligands: dict[RCA_Ligand] = None,
                 metal_charge: int = 2,
                 metal: str = "Fe",
                 spin: int = 0
                 ):
        """
        In general, during the creation process we get
        an stk.ConstructedMolecule with the correct geometries, but with Hg to be remaining
        a dict of ligands, which contain the ligands with old geometries, but the desireable graphs
        """
        if ligands is None or compl is None:
            # we need a way to create an empty object
            return

        # from the stk compl we need to extract the atoms and their 3d positions
        # i.e. we want the atomic_properties for our new Molecule
        atomic_props = self.stk_Constructed_Mol_to_atomic_props(compl)

        # we will ignore global properties, because this doesnt really make sense.
        # however, we will add an attribute which we call ligand_properties, where some information on the ligands is stored
        # I like to store everyhting regarding the ligands in this dict
        self.ligand_props = {
            key: {
                "name": ligand.name,
                "unique name": ligand.unique_name if hasattr(ligand, "unique_name") else None,
                "original metal": original_metal_ligand(ligand),
                "size": len(ligand.atomic_props["x"]),
                "stoichiometry": ligand.stoichiometry,
                "Metal Spin:": spin,
            }
            for key, ligand in ligands.items()
        }

        self.total_charge = self.get_total_charge(metal_charge, ligands)

        self.graph = self.merge_graph_from_ligands(ligands, metal)

        self.mol = RCA_Molecule.make_from_atomic_properties(
            atomic_props_mol=atomic_props,
            global_props_mol={},
            graph=self.graph
        )

        self.date_of_creation = date.today()
        self.name = self.assemble_name(metal=metal, ligands=ligands)

        self.functional_groups = {key_: lig.get_assembly_dict()["type"] for key_, lig in ligands.items()}

        self.working_name = self.create_random_name()

        self.metal = metal

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
        indices_for_non_Hg = [i for i, atom in enumerate(compl._atoms) if atom.get_atomic_number() != atomic_number_Hg]

        # Now we can extract the types of atoms
        atomic_props["atoms"] = [element(compl._atoms[i].get_atomic_number()).symbol for i in indices_for_non_Hg]

        for (x, y, z) in compl.get_atomic_positions(indices_for_non_Hg):
            atomic_props["x"].append(x)
            atomic_props["y"].append(y)
            atomic_props["z"].append(z)

        return atomic_props

    @staticmethod
    def merge_graph_from_ligands(ligands, metal):
        """
        We now have to carefully assemble the graph
        """

        Gs = [deepcopy(lig.graph) for lig in ligands.values()]

        # first thing is that we have to relabel the nodes, so that we have no repititions
        # otherwise the merging wont work
        G = nx.Graph()
        G.add_nodes_from([
            (0, {"node_label": metal}),
        ])

        # relabel all graphs
        i = 1
        for H in Gs:
            node_mapping = {node: i + k for k, node in enumerate(H.nodes)}
            nx.relabel_nodes(H, mapping=node_mapping, copy=False)
            i += len(H.nodes)

        # now we create the new graph by merging everything
        U = nx.Graph()
        U.add_nodes_from(G.nodes(data=True))
        for H in Gs:
            U.add_edges_from(H.edges())
            U.add_nodes_from(H.nodes(data=True))

        # finally we have to bring in the new edges
        for Gr, lig in zip(Gs, ligands.values()):
            for i in lig.ligand_to_metal:
                U.add_edge(0, sorted(Gr.nodes)[i])

        # outputs graph of assembled complex
        return U

    def create_random_name(self):

        # we first generate the reconstructable hash:
        name_hash = int(hashlib.md5(self.name.encode(encoding='UTF-8', errors='strict')).hexdigest(), 16)

        # this is reproducable and 32 digits long
        if int(str(name_hash)[:10]) % 2 == 1:
            gender = "Male"
        else:
            gender = "Female"

        name_row = int(str(name_hash)[10:20]) % 200
        # because we have currently 200 names implemented

        surname_row = int(str(name_hash)[20:30]) % 1000

        # pick the name
        df = pd.read_csv(Path(project_path, "constants", "names", "names.csv"))

        return f"{df.loc[name_row, gender]} {df.loc[surname_row, 'Surname']}"

    def to_json(self,
                path: Union[Path, str]
                ):

        complex_properties = deepcopy(self.__dict__)

        complex_properties["mol"] = complex_properties["mol"].write_to_mol_dict()
        # However, the graph dicts are already contained in the above entry of the dict
        del complex_properties["graph"]

        complex_properties["date_of_creation"] = str(complex_properties["date_of_creation"])

        with open(path, "w") as file:
            json.dump(complex_properties, file)

    @classmethod
    def from_json(cls,
                  path: Union[Path, str]
                  ):
        """
        The inverse method of "to_json"
        """

        a = cls()
        # vielleicht mit: in der init: if complex is None: pass -> leeres objekt und dann hier fuellen

        with open(path, "r") as file:
            properties = json.load(file)

        if not {'ligand_props', 'total_charge', 'mol', 'date_of_creation', 'name', 'functional_groups',
                'working_name'}.issubset(set(properties.keys())):
            raise KeyError("Missing keys in input json for TMC generation")

        a.mol = RCA_Molecule.read_from_mol_dict(properties["mol"])
        a.graph = a.mol.graph

        del properties["mol"]

        for k, v in properties.items():
            a.__setattr__(k, v)

        return a

    def is_equal(self, other):
        """
        other: also an object of type TMC,
        by checking if all attributes are the same
        """

        for key in self.__dict__:

            # comparison for the string valued stuff
            if key in ["total_charge", "name", "working_name", "metal"]:
                try:
                    if self.__dict__[key] == other.__dict__[key]:
                        continue
                    else:
                        return False
                except KeyError:
                    # because then the attribute is missing in the other
                    # and hence they can't be equal
                    return False
            elif key == "date_of_creation":
                try:
                    if str(self.__dict__[key]) == str(other.__dict__[key]):
                        continue
                    else:
                        return False
                except KeyError:
                    # because then the attribute is missing in the other
                    # and hence they can't be equal
                    return False
            elif key in ["ligand_props", "functional_groups"]:
                for k, v in self.__dict__[key].items():
                    try:
                        if other.__dict__[key][k] == v:
                            pass
                        else:
                            return False
                    except KeyError:
                        try:
                            if other.__dict__[key][str(k)] == v:
                                pass
                            else:
                                return False
                        except KeyError:
                            return False
            elif key == "graph":
                try:
                    if graphs_are_equal(other.__dict__[key], self.__dict__[key]):
                        continue
                    else:
                        return False
                except KeyError:
                    return False
            elif key == "mol":
                try:
                    if self.__dict__[key].get_graph_hash() == other.__dict__[key].get_graph_hash():
                        continue
                    else:
                        return False
                except KeyError:
                    return False
            else:
                continue

        return True

    def equal_graphs(self, other):
        try:
            if self.__dict__["mol"].get_graph_hash() == other.__dict__["mol"].get_graph_hash():
                return True
            else:
                return False
        except KeyError:
            print("Warning: One of the TMCs has no .mol attribute")
            return False

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



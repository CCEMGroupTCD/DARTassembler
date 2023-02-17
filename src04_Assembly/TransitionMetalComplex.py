from datetime import date
import stk
from mendeleev import element
from copy import deepcopy
import networkx as nx
import hashlib
import pandas as pd
from src01.utilities_Molecule import original_metal_ligand
from src01.Molecule import RCA_Molecule, RCA_Ligand
atomic_number_Hg = 80

# todo:
# Auf dem complex fehlen noch ein paar Methoden:
# 1. eine .to_json
# 3. eine die ein Gaussian erstellt.


class TransitionMetalComplex:
    """
    Here we define how we want the final object to look like
    to bring us in a good position for the post assembly filtering
    """
    mol: RCA_Molecule

    def __init__(self,
                 compl: stk.ConstructedMolecule,
                 ligands: dict[RCA_Ligand],
                 metal_charge: int,
                 metal: str,
                 spin: int
                 ):
        """
        In general, during the creation process we get
        an stk.ConstructedMolecule with the correct geometries, but with Hg to be remaining
        a dict of ligands, which contain the ligands with old geometries, but the desireable graphs
        """

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
                "charge": ligand.global_props["charge"] if "charge" in ligand.global_props else None,
                "size": len(ligand.atomic_props["x"]),
                "stoichiometry": ligand.stoichiometry,
                "Metal Spin:": spin,
            }
            for key, ligand in ligands.items()
        }

        try:
            charge_list = []
            for ligand in ligands.values():
                charge = ligand.global_props['LCS_pred_charge']
                charge_list.append(int(charge))
            self.total_charge = metal_charge + sum(charge_list)
            #self.total_charge = sum([metal_charge] + [ligand.global_props['LCS_pred_charge'] for ligand in ligands.values()])
        except KeyError:
            # in this case we might not have all of the charges
            self.total_charge = None

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

    def assemble_name(self, metal, ligands):
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
        df = pd.read_csv("../constants/names/names.csv")

        return f"{df.loc[name_row, gender]} {df.loc[surname_row, 'Surname']}"



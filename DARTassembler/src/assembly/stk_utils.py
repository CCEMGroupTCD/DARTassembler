import numpy as np
import stk
from DARTassembler.src.constants.Periodic_Table import DART_Element
import networkx as nx

def create_placeholder_Hg_bb() -> stk.BuildingBlock:
    """
    We use that frequently throughout so we shall return that special building block
    center atom with 6 connections
    """
    return stk.BuildingBlock(smiles='[Hg+2]',
                             functional_groups=(stk.SingleAtom(stk.Hg(0, charge=2)) for i in range(6)),
                             position_matrix=np.ndarray([0, 0, 0])
                             )


def convert_RCA_to_stk_Molecule(mol):
    """
    :param mol: RCA_Molecule
    """

    # create a list of atoms
    atom_list = [stk.Atom(id=i, atomic_number=DART_Element(atom).atomic_number) for i, atom in enumerate(mol.atomic_props["atoms"])]

    # Now we need the bonds from the graph
    # in fact this is good as we have full control that the stk molecules look according to the graphs
    # however, all of them will have default bond order 1
    edges = [e for e in mol.graph.edges]
    nodes = [n for n in mol.graph.nodes]

    bond_list = [stk.Bond(atom1=atom_list[nodes.index(i)], atom2=atom_list[nodes.index(j)], order=1) for (i,j) in edges]

    # Finally the position matrix remains
    A = [[mol.atomic_props["x"][i], mol.atomic_props["y"][i], mol.atomic_props["z"][i]] for i, _ in enumerate(mol.atomic_props["x"])]

    return stk.Molecule(atoms=atom_list, bonds=bond_list, position_matrix=np.array(A))


def RCA_Mol_to_stkBB(mol):
    """
    :param mol: RCA_Ligand
    """

    stk_mol = convert_RCA_to_stk_Molecule(mol)

    func_dict = {type_: getattr(stk, type_) for type_ in mol.get_assembly_dict()["type"]}
    atoms_ = [func_dict[type_](mol.get_assembly_dict()["index"][i]) for i, type_ in enumerate(mol.get_assembly_dict()["type"])]

    functional_groups_ = [stk.GenericFunctionalGroup(atoms=(a,), bonders=(a,), deleters=(),) for a in atoms_]

    return stk.BuildingBlock.init_from_molecule(molecule=stk_mol,
                                                functional_groups=functional_groups_
                                                )


def stkBB_to_networkx_graph(stk_building_block, check_fully_connected=True) -> nx.Graph:
    """
    Convert an stk.BuildingBlock to a networkx graph
    :param stk_building_block: stk.BuildingBlock
    """
    # Initialize an empty graph
    G = nx.Graph()

    # Add atoms as nodes
    for atom in stk_building_block.get_atoms():
        G.add_node(atom.get_id(), element=atom.__class__.__name__)

    # Add bonds as edges
    for bond in stk_building_block.get_bonds():
        atom1 = bond.get_atom1().get_id()
        atom2 = bond.get_atom2().get_id()
        G.add_edge(atom1, atom2)

    if check_fully_connected:
        assert nx.is_connected(G), "The graph is not fully connected!"

    return G







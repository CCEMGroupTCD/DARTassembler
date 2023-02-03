# # # Imply our own script to compute the RACS
from pysmiles import read_smiles
import networkx as nx
import matplotlib.pyplot as plt
from mendeleev import element
inf = "inf"


def smiles_to_graph(smiles_str: str):
    """
    :param smiles_str: smiles str
    :return: adjacecny matrix for smiles str
    """
    mol = read_smiles(smiles_str, explicit_hydrogen=True)
    #print(mol.nodes(data='element'))
    #nx.draw_networkx(mol, node_size=500, with_labels=True, labels={label[0]: label[-1] for label in mol.nodes(data='element')})
    #plt.show()

    A = nx.adjacency_matrix(mol).todense()
    label_dict = {label[0]: label[-1] for label in mol.nodes(data='element')}

    return A, label_dict


def get_graph_distance_from_adj_matrix(i, j, A):
    """
    :param i: Node one
    :param j: Node two
    :param A: adjacency matrix of a graph G
    :return: distance of the nodes in graph G (0 if i=j)
    """
    B = A

    for k in range(len(A)):
        B = B*A
        if B[i, j] > 0:
            return k

    # should never occur as we always assume connected graphs
    return inf


def compute_RAC(Adjacecny_matrix, distance, property_, label_dict):
    """
    :param label_dict:
    :param Adjacecny_matrix:
    :param distance:
    :param property_:
    :return:
    """

    print(property_)
    # bring the property to the right format
    if property_ == "nuclear_charge":
        property_ = "protons"
    elif property_ == "electronegativity":
        property_ = "electronegativity_pauling"

    RAC_value = 0

    for i in range(len(Adjacecny_matrix)):
        for j in range(len(Adjacecny_matrix)):
            if get_graph_distance_from_adj_matrix(i=i, j=j, A=Adjacecny_matrix) == distance:

                if property_ == "identity":
                    RAC_value += 1*1
                else:
                    el_i, el_j = element(label_dict[i]), element(label_dict[j])
                    try:
                        P_i, P_j = getattr(el_i, property_), getattr(el_j, property_)

                        RAC_value += P_i * P_j
                    except TypeError:
                        P_i, P_j = getattr(el_i, property_)(), getattr(el_j, property_)()

                        RAC_value += P_i * P_j

                    except Exception as e:
                        print(f"Something else: {e}")

    return RAC_value


def compute_RACs_from_smiles(smiles_str: str, property_="atomic_number", range_=None) -> dict:

    A, label_dict = smiles_to_graph(smiles_str=smiles_str)

    RAC_dict = {}

    if property_ not in ["atomic_number", "electronegativity", "nuclear_charge", "identity", "atomic_radius"]:
        print("No RACs for the descired property implemented (yet)")
        raise NotImplementedError

    if range_ is None:
        range_ = (0, len(A))

    assert isinstance(range_, tuple), "Wrong format for the range"

    for d in range(range_[0], range_[1]):
        RAC_dict[d] = compute_RAC(Adjacecny_matrix=A,
                                  distance=d,
                                  property_=property_,
                                  label_dict=label_dict
                                  )

    return RAC_dict


def compute_full_RACs(smiles_str_, min_distance=None, max_distance=None):
    """
    :param max_distance:
    :param min_distance:
    :param smiles_str_:
    :return: a Dict {property: {distance: RAC_value} }
    """

    full_RAC_dict = {}
    if max_distance is None:
        range_ = None
    else:
        range_ = (0, max_distance)

    for prop in ["atomic_number", "electronegativity", "nuclear_charge", "identity", "atomic_radius"]:
        full_RAC_dict[prop] = compute_RACs_from_smiles(smiles_str=smiles_str_,
                                                       property_=prop,
                                                       range_=range_
                                                       )

    return full_RAC_dict


if __name__ == "__main__":
    smiles_str = "CH#CH"

    A, _ = smiles_to_graph(smiles_str)

    RAC_dict = compute_full_RACs(smiles_str_=smiles_str)
    print("Done")


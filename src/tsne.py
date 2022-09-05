from collections import defaultdict
import ASE_Molecule
import itertools

import pickle
import os
import matplotlib.pyplot as plt

from molSimplify.Informatics.autocorrelation import *
from molSimplify.Classes import mol3D

from sklearn.manifold import TSNE
from sklearn.cluster import KMeans


def get_ligand_RAC(ligand: ASE_Molecule.ASE_Ligand, depth, properties):
    '''
    :param ligand: object of class ASE_Molecule
    :param properties: properties for which the RACs should be computed
    :param depth: depth for the RACs
    :return: # todo: mal checken wie das aussieht
    '''

    with open("data/tmp2.xyz", "w+") as f:
        f.write(ligand.ligand_xyz_str)

    mol_ = mol3D.mol3D()
    mol_.readfromxyz("data/tmp2.xyz")
    os.remove("data/tmp2.xyz")

    RACs = np.zeros(shape=(len(properties), depth + 1))
    delta = np.zeros(shape=(len(properties), depth))
    for idx, p in enumerate(properties):
        RACs[idx] = atom_only_autocorrelation(mol_, p, depth, ligand.ligand_to_metal, oct=False)
        delta[idx] = atom_only_deltametric(mol_, p, depth, ligand.ligand_to_metal)[1:]
    rac_res = np.concatenate((RACs, delta), axis=None)

    return rac_res


def get_updated_ligand_cluster_dict(ligand_cluster_dict, ligname2rac, n_clusters=100, verbose=False):

    # müssen irgendwie die ligands addressieren können
    index_to_lig =  {}
    lig_racs = list()
    for i, lig in enumerate(list(ligand_cluster_dict.keys())):
        lig_racs.append(ligname2rac[lig.name])
        index_to_lig[i] = lig

    # Dimension reduction via tSNE, to dim=2 (Hyperparameter 1,2)
    n_comps = 2
    X_embedded = TSNE(n_components=n_comps).fit_transform(lig_racs)

    kmc = KMeans(n_clusters=n_clusters, random_state=42)
    kmc.fit(X_embedded)
    if verbose:
        plt.scatter(X_embedded[:, 0], X_embedded[:, 1], alpha=0.6, c=kmc.labels_)
        plt.show()

    # make a dictionary mapping cluster number to indices which are within that cluster
    clus2indices = defaultdict(list)
    for idx, clus in enumerate(kmc.labels_):
        clus2indices[clus].append(idx)

    for cluster_prime in clus2indices.values():
        if len(cluster_prime) == 1:
            continue
        else:
            tuples = list(itertools.combinations(cluster_prime, 2))
            for (idx_1, idx_2) in tuples:
                lig1 = index_to_lig[idx_1]
                lig2 = index_to_lig[idx_2]
                if len(lig1.mol.get_atomic_numbers()) == len(lig2.mol.get_atomic_numbers()):
                    '''
                    If they are in the same cluster, and have the same number of atoms, 
                    we consider them to be equal.
                    '''
                    ligand_cluster_dict[lig2] = ligand_cluster_dict[lig1]

    return ligand_cluster_dict


def get_cluster_from_dict(dict_):

    new_dict = {}
    for lig, cluster_num in dict_.items():
        if cluster_num in new_dict.keys():
            new_dict[cluster_num].append(lig)
        else:
            new_dict[cluster_num] = [lig]

    return new_dict


if __name__ == '__main__':

    dent = 4    # only tetradentates for test reasons
    depth = 5
    properties = ["electronegativity", "nuclear_charge"]

    with open("Old/old_data/ligand_list.pickle", "rb") as handle:
        ligand_list = pickle.load(handle)

    # only certain denticity and no small ligands
    ligand_list_small = [lig for lig in ligand_list if lig.denticity == dent and len(lig.mol.get_atomic_numbers()) >= 10]
    del ligand_list     # to restore some memory

    ligname2rac = {lig.name: get_ligand_RAC(ligand=lig, depth=depth, properties=properties) for lig in ligand_list_small}

    n_clus = 10      # für test, später größer; man könnte das als funktion der Länge von ligand_list_small festhalten
    ligand_cluster_dict = {lig: i for i, lig in enumerate(ligand_list_small)}

    while True:

        ligand_cluster = get_cluster_from_dict(ligand_cluster_dict)         #getting the actual cluster from the dict
        n_old_ligand_cluster = len(ligand_cluster)

        ligand_cluster_dict = get_updated_ligand_cluster_dict(ligand_cluster_dict, ligname2rac, n_clusters=n_clus, verbose=True)

        ligand_cluster = get_cluster_from_dict(ligand_cluster_dict)

        if n_old_ligand_cluster - len(ligand_cluster) == 0:
            # if clustering doesnt anything: end
            break
        # only after the clustering doesn't change anything, we break

    final_cluster = get_cluster_from_dict(ligand_cluster_dict)



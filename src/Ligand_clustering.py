import pickle
import ASE_Molecule
import os
from collections import defaultdict
import itertools
from tqdm import tqdm

import matplotlib.pyplot as plt

from molSimplify.Informatics.autocorrelation import *
from molSimplify.Classes import mol3D

from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def get_ligand_RAC(ligand: ASE_Molecule.ASE_Ligand, depth, properties):
    """
    :param ligand: object of class ASE_Molecule
    :param properties: properties for which the RACs should be computed
    :param depth: depth for the RACs
    :return:
    """

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


def dim_red(lig_racs):
    """
    specify the method of dimension reduction
    Dimension reduction via tSNE, to dim=2 (Hyperparameter 1,2)
    """
    n_comps = 2
    return TSNE(n_components=n_comps).fit_transform(lig_racs)


def get_updated_ligand_cluster_dict(ligand_cluster_dict, ligname2rac, n_clusters=100, verbose=False):

    # müssen irgendwie die ligands addressieren können
    index_to_lig = {}
    lig_racs = list()
    for i, lig in enumerate(list(ligand_cluster_dict.keys())):
        lig_racs.append(ligname2rac[lig.name])
        index_to_lig[i] = lig

    X_embedded = dim_red(lig_racs)

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


def optimal_number_clusters(ligand_cluster_dict, ligname2rac):

    silhouette_scores = list()

    # smallest number of clusters:
    min_cluster_n = 10

    for n in tqdm(range(min_cluster_n, len(ligand_cluster_dict)+1)):
        try:
            lig_racs = list()
            for i, lig in enumerate(list(ligand_cluster_dict.keys())):
                lig_racs.append(ligname2rac[lig.name])

            X_embedded = dim_red(lig_racs)

            kmc = KMeans(n_clusters=n)
            kmc.fit(X_embedded)
            cluster_labels = kmc.labels_

            # silhouette score
            silhouette_scores.append(silhouette_score(X_embedded, cluster_labels))
        except Exception:
            pass

    return silhouette_scores.index(max(silhouette_scores)) + 10


def cluster_ligands(dent_number: int, target_path: str, depth=5, verbose=False):

    properties = ["electronegativity", "nuclear_charge"]

    with open(f"{target_path}/ligand_list.pickle", "rb") as handle:
        ligand_list = pickle.load(handle)

    # only certain denticity and no small ligands
    ligand_list_small = [lig for lig in ligand_list if lig.denticity == dent_number and len(lig.mol.get_atomic_numbers()) >= 10]
    del ligand_list     # to restore some memory

    ligname2rac = {lig.name: get_ligand_RAC(ligand=lig, depth=depth, properties=properties) for lig in ligand_list_small}
    ligand_cluster_dict = {lig: i for i, lig in enumerate(ligand_list_small)}

    n_clus = optimal_number_clusters(ligand_cluster_dict, ligname2rac)
    if verbose:
        print(f"number of clusters: {n_clus} with a total use of {len(ligand_cluster_dict)} ligands")

    while True:

        ligand_cluster = get_cluster_from_dict(ligand_cluster_dict)         # getting the actual cluster from the dict
        n_old_ligand_cluster = len(ligand_cluster)

        ligand_cluster_dict = get_updated_ligand_cluster_dict(ligand_cluster_dict, ligname2rac, n_clusters=n_clus, verbose=False)

        ligand_cluster = get_cluster_from_dict(ligand_cluster_dict)

        if n_old_ligand_cluster - len(ligand_cluster) == 0:
            # if clustering doesnt anything: end
            break
        # only after the clustering doesn't change anything, we break

    # Now we have the final dict
    # cluster_number : [list of ASE_ligands] w.r.t to a certain denticity
    final_cluster = get_cluster_from_dict(ligand_cluster_dict)

    with open(f"{target_path }/clustered_ligands_w_dent_{dent_number}.pickle", "wb") as handle:
        pickle.dump(final_cluster, handle)

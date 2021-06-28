"""
For more informations on the contents of this module:
- help(clustering.GenotypeMatrix)
- help(clustering.cluster_mutations)

--------

Module that exposes the clustering algorithm presented at
https://github.com/AlgoLab/celluloid

Simple example workflow:

from plastic import clustering as cl

to_cluster = cl.GenotypeMatrix.from_files('to_cluster.txt', mutations_file = 'mutations.txt')
# Reduce the size of the input down to 50 to speed up some complex computation
# (for instance SASC tree inference)
clustered = cl.cluster_mutations(to_cluster, k = 50)

# Get the clustered mutations as comma separated lists of simple mutations
muts = clustered.mutations()

# Save the matrix and use it for some intensive computation
cl.GenotypeMatrix.to_files('clustered.txt', mutations_file = 'clustered_mutations.txt')
from phylo import sasc as sc
tree = sc.infer_tree(clustered, ...)['inferred_tree']
"""

from ._core.genotypematrix import GenotypeMatrix
import numpy as np
from kmodes.kmodes import KModes
from collections import defaultdict


def cluster_mutations(
        genotype_matrix,
        k,
        n_inits=10,
        max_iter=100,
        verbose=False,
        **kwargs):
    """
    Clusters the mutations in a genotype matrix by applying kmodes

    Parameters:
        genotype_matrix(GenotypeMatrix):
            A matrix representing the results of single-cell sequencing.
        k(int):
            The number of clustered mutations in the output matrix.
            Note that empty clusters will be discarded after clustering.
        n_inits(int):
            The number of initiliazations in the clustering process.
        max_iter(int):
            The maximum number of iterations in the clustering process.
        verbose (bool)
        **kwargs:
            Additional arguments passed to KModes process.

    Returns:
        GenotypeMatrix:
            The result of the clustering process. Each column in the matrix
            will be the centroid of a non-empty cluster, and will be labeled with
            a comma-separated list of the labels of the mutations within the cluster.
            Cell labels are left unaltered.
    """

    if type(k) != int or k < 1:
        raise ValueError(f'the number of clusters must be a positive integer, but {k} is not.')
    if type(max_iter) != int or max_iter < 1:
        raise ValueError(f'the number of iterations must be a positive integer, but {max_iter} is not.')
    if type(n_inits) != int or n_inits < 1:
        raise ValueError(f'the number of initializations must be a positive integer, but {n_inits} is not.')

    return _celluloid(genotype_matrix, k, n_inits, max_iter,verbose,**kwargs)

def _conflict_dissim(a, b, **_):
    v = np.vectorize(lambda ai, bi: ai != 2 and bi != 2 and ai != bi)
    return np.sum(v(a, b), axis=1)


def _celluloid(
        genotype_matrix,
        k,
        n_inits,
        max_iter,
        verbose,
        **kwargs
):
    """
    Clusters the mutations in a genotype matrix by applying kmodes

    Parameters:
        genotype_matrix(GenotypeMatrix):
            A matrix representing the results of single-cell sequencing.
        k(int):
            The number of clustered mutations in the output matrix.
            Note that empty clusters will be discarded after clustering.
        n_inits(int):
            The number of initiliazations in the clustering process.
        max_iter(int):
            The maximum number of iterations in the clustering process.
        verbose (bool)
        **kwargs:
            Additional arguments passed to KModes process.

    Returns:
        GenotypeMatrix:
            The result of the clustering process. Each column in the matrix
            will be the centroid of a non-empty cluster, and will be labeled with
            a comma-separated list of the labels of the mutations within the cluster.
            Cell labels are left unaltered.
    """

    mutations_as_points = np.array(genotype_matrix.matrix(), dtype='int').transpose()
    mutation_labels = genotype_matrix.mutation_labels

    km = KModes(
        n_clusters=k,
        cat_dissim=_conflict_dissim,
        init='huang',
        n_init=n_inits,
        max_iter=max_iter,
        verbose=(1 if verbose else 0),
        **kwargs
    )

    clusters = km.fit_predict(mutations_as_points)

    # Each cluster will be labeled with the labels of its components.
    clusters_of_mutations = km.labels_
    clustered_mutation_labels = defaultdict(list)
    for mutation_label, mutation_cluster in zip(mutation_labels, clusters_of_mutations):
        clustered_mutation_labels[mutation_cluster].append(mutation_label)

    nonempty_clusters = clustered_mutation_labels.keys()

    # build the output matrix and the mutation labels as strings
    cluster_centroids = km.cluster_centroids_
    clustered_mutation_labels_strings = [','.join(clustered_mutation_labels[cluster_id]) for cluster_id in
                                         sorted(nonempty_clusters)]
    out_matrix = [cluster_centroids[cluster_id] for cluster_id in sorted(nonempty_clusters)]

    # the matrix needs to be transposed back to its original orientation
    out_matrix = np.array(out_matrix).transpose()

    return GenotypeMatrix(out_matrix, cell_labels=genotype_matrix.cell_labels,
                          mutation_labels=clustered_mutation_labels_strings)

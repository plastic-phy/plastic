"""
For more info about this module's functionalities:
- help(phylogeny.SASCPhylogeny)
- help(plastic.GenotypeMatrix)
- help(phylogeny.infer_tree)

---------

Module that exposes the phylogeny inference algorithm presented at
https://github.com/sciccolella/sasc

Simple example workflow:
from phylo import sasc as sc
from multiprocessing import cpu_count

# load the matrix and the mutation labels
gmat = sc.GenotypeMatrix.from_files(
    matrix_file = 'genotype.txt', 
    mutations_file = 'mutation_labels.txt', 
    matrix_parser = sc.SASCParser() # sc.SASCParser, sc.SPHYRParser and sc.SCITEParser are available
                                    # but anything with a parse() method that puts out a matrix is fine
)

# infer a tree and get the expected genotype matrix for it, then store these to files.
out = sc.infer_tree(gmat, alphas = 0.3, beta = 0.1, k = 1, cores = cpu_count())
out['inferred_tree'].to_file('tree.gv')
out['inferred_tree'].with_visualization_features() \
    .draw_to_file('tree.png', show_support = True, show_color = False)
out['expected_genotype_matrix'].to_files('expected_genotype.txt')
"""


from libc.stdlib cimport malloc, free
from libc.stdio cimport printf
from plastic._core.genotypematrix import GenotypeMatrix
from plastic._core.phylogenytree import SASCPhylogeny
import numpy as np
import pandas as pd
import networkx as nx
cimport sascpycapi as sca
from libc.limits cimport INT_MAX


def sasc(
        labeled_genotype_matrix,
        alphas,
        beta,
        k,
        max_deletions=INT_MAX,
        repetitions=5,
        start_temp=10**4,
	    cores=1,
        el_a_variance=0,
        el_b_variance=0,
        el_g_variance=0,
        monoclonal=False,
        gammas=1,
        get_cells=False,
        cooling_rate=0.01 #this cannot be written as 10**-2, otherwise Cython builds it from 10L and assigns 0 to it
):
    """
    Infers a tree from a matrix representing the mutations for a sample of cells using
    simulated annealing to build a phylogeny using a k-Dollo model.

    Parameters:
        labeled_genotype_matrix(GenotypeMatrix):
            The genotype matrix from which to infer the phylogeny tree.
        alphas(float or list(float)):
            The probability of a false negative for each mutation (that is, for each column in the matrix).
            this can be a single float that will be used for all mutations, or a list of floats
            with an element for each mutation.
        gammas(float or list(float)), by default 1:
            The prior loss probability for each mutation. If a single float is specified,
            then it will be used for all mutations, and if nothing is specified the used value will be 1.
            Otherwise, a list of floats must be provided with a gamma value for each mutation.
        beta(float):
            The probability for a false positive for each mutation. A single value must be provided.
        k(int):
            The k used in the k-Dollo model for tree inference. This is the maximum amount of times
            each mutation can be lost after being acquired in the phylogeny tree.
        max_deletions(int), by default MAX_INT:
            The maximum amount of deletions globally allowed in the tree.
        repetitions(int), by default 5:
            The number of times the inference process will be repeated.
            The inferred tree with the best score will be used for the output.
        monoclonal(bool), by default False:
            If this is true, then SASC will be forced to infer a monoclonal tree,
            with exactly one child for the germline.
        start_temp(float), by default 10000:
            The starting temperature for simulated annealing.
        cooling_rate(float), by default 0.01:
            The cooling rate for simulated annealing.
        cores(int), by default 1:
            The number of cores used by SASC.
        el_a_variance(float), by default 0:
            The variance of the false negative rates for error learning.
            If a single alpha was specified, then all the false negative rates will be adjusted
            simultaneously. Otherwise, each alpha will change independently from the others.
        el_b_variance(float), by default 0:
            The variance of the false positive rates for error learning.
            If a single gamma was specified, then all the prior loss probabilities will be adjusted
            simultaneously. Otherwise, each gamma will change independently from the others.
        el_g_variance(float), by default 0:
            The variance of the prior loss probabilities for error learning.
        get_cells(bool), by default False:
            If this is true, the tree in the output will also include the nodes
            for the cells in the sample.

    Returns:
        dict:
            A dictionary with the following key-value pairs:
            * inferred_tree (PhylogenyTree):
                The tree that was inferred by SASC, with its confidence score.
            * expected_genotype_matrix (GenotypeMatrix):
                A matrix where the missing information about the mutations has been inferred
                starting from the tree.
            * inferred_alphas (list(float)):
                The false negative rates inferred through error learning.
            * inferred_beta (float):
                The false positive rate inferred through error learning.
            * inferred_gammas (list(float)):
                The prior loss probabilities inferred through error learning.
    """


    cdef sca.sasc_in_t args_struct
    cdef sca.sasc_in_t* arguments = &args_struct

    genotype_matrix = labeled_genotype_matrix.matrix()
    arguments.N = genotype_matrix.shape[0]
    arguments.M = genotype_matrix.shape[1]
    cdef int N = arguments.N
    cdef int M = arguments.M

    #input validation
    if beta < 0 or beta > 1:
        raise ValueError(f'beta must be within 0 and 1, but {beta} is not')
    if el_a_variance < 0 or el_b_variance < 0 or el_g_variance < 0:
        raise ValueError(f'error learning parameters must be positive, but one in {[el_a_variance, el_g_variance, el_b_variance]} is not')
    if int(k) != k or k < 0:
        raise ValueError(f'k must be a non-negative integer, but {k} is not')
    if int(max_deletions) != max_deletions  or max_deletions < 0:
        raise ValueError(f'the maximum number of deletions must be a non-negative integer, but {max_deletions} is not')
    if int(repetitions) != repetitions or repetitions < 1:
        raise ValueError(f'the number of repetitions must be a positive integer, but {repetitions} is not')
    if start_temp < 0:
        raise ValueError(f'the starting temperature must be positive, but {start_temp} is not')
    if cooling_rate <= 0:
        raise ValueError(f'the cooling rate must be positive, but {cooling_rate} is not')
    if int(cores) != cores or cores < 1:
        raise ValueError(f'the number of cores must be a positive integer, but {cores} is not')

    single_alpha = not isinstance(alphas, list)
    single_gamma = not isinstance(gammas, list)
    if single_alpha:
        alphas = [alphas] * M
    if single_gamma:
        gammas = [gammas] * M
    if len(alphas) != M:
        raise ValueError(f'multiple alphas are specified in {alphas}, but they are more or less than the number of mutations.')
    if len(gammas) != M:
        raise ValueError(f'multiple gammas are specified in {gammas}, but they are more or less than the number of mutations.')
    if any([alpha < 0 or alpha > 1 for alpha in alphas]):
        raise ValueError(f'all alpha values must be within 0 and 1, but at least one in {alphas} is not.')
    if any([gamma < 0 or gamma > 1 for gamma in gammas]):
        raise ValueError(f'all gamma values must be within 0 and 1, but at least one in {gammas} is not.')

    # Some inputs need to be processed and marshalled "manually" before being fed to the c function.
    # Let Cython handle the rest of them by itself.
    
    # Automatic marshalling.
    arguments.beta = beta #;print(beta)
    arguments.el_a_variance = el_a_variance
    arguments.el_b_variance = el_b_variance
    arguments.el_g_variance = el_g_variance
    if k == 0 or max_deletions == 0:
        k = 0
        max_deletions = 0
    arguments.k = k
    arguments.max_deletions = max_deletions
    arguments.repetitions = repetitions
    arguments.start_temp = start_temp
    arguments.cooling_rate = cooling_rate
    arguments.cores = cores

    # Marshalling of the matrix
    arguments.genotype_matrix = <int**>malloc(N*sizeof(int*))
    if arguments.genotype_matrix == NULL:
        raise MemoryError()
    for i in range(N):
        arguments.genotype_matrix[i] = <int*>malloc(M*sizeof(int))
        if arguments.genotype_matrix[i] == NULL:
            raise MemoryError()
        for j in range(M):
            arguments.genotype_matrix[i][j] = genotype_matrix[i][j]

    # To work around the string size limitation for SASC, we pass their indexes in the lists in string form to SASC
    # instead of the labels.
    arguments.mutation_labels = <char**>malloc(M * sizeof(char*))
    automatic_labels = [bytes(str(i), 'ascii') for i in range(M)]
    for i in range(M):
        arguments.mutation_labels[i] = automatic_labels[i]

    # Arrays with error parameters must be allocated and filled
    arguments.alphas = <double*>malloc(M * sizeof(double))
    arguments.gammas = <double*>malloc(M * sizeof(double))
    if (
        arguments.gammas == NULL
        or arguments.alphas == NULL
    ):
        raise MemoryError()
    
    for i in range(M):
        arguments.alphas[i] = alphas[i]
        arguments.gammas[i] = gammas[i]

    # Python bools must be converted into C integers
    arguments.single_alpha = 1 if single_alpha else 0
    arguments.single_gamma = 1 if single_gamma else 0
    arguments.monoclonal = 1 if monoclonal else 0

    # The output structure is created beforehand to avoid memory allocations inside the c call.
    cdef sca.sasc_out_t out_struct
    cdef sca.sasc_out_t* c_out = &out_struct

    c_out.gtp_matrix = <int**> malloc(N * sizeof(int*))
    if c_out.gtp_matrix == NULL:
        raise MemoryError()
    for i in range(N):
        c_out.gtp_matrix[i] = <int*> malloc(M * sizeof(int))
        if c_out.gtp_matrix[i] == NULL:
            raise MemoryError()
    
    c_out.ids_of_leaves = <int*> malloc(N * sizeof(int))
    c_out.el_alphas = <double*> malloc(M * sizeof(double))
    c_out.el_gammas = <double*> malloc(M * sizeof(double))
    if (
        c_out.el_alphas == NULL
        or c_out.el_gammas == NULL
        or c_out.ids_of_leaves == NULL
    ):
        raise MemoryError()

    cdef int comp_result
    cdef sca.node_t* root

    try:
        # THE C A L L
        comp_result = sca.compute(arguments, c_out)
    
        # And now for the unmarshalling. We'll output a tuple with the tree as a networkx graph, the matrix as
        # a numpy array, and the rest of the values as simple ints/doubles/strings.
    
        # Unmarshalling the tree.
        best_tree = nx.DiGraph()
        _unmarshal_tree(c_out.best_tree, best_tree)

        # The labels for the tree are now indexes to the original mutation labels, aside from special labels that are
        # added artificially to the tree during the execution of SASC and only need to be converted
        # to a string.
        # We assume that every label that cannot be converted to an integer needs to be left as-is.
        def converts_to_index(label):
            try:
                return int(label) >= 0
            except ValueError:
                return False
    
        for [node, attributes] in best_tree.nodes(data = True):
            if 'label' in attributes:
                if converts_to_index(attributes['label']):
                    attributes['label'] = labeled_genotype_matrix.mutation_labels[ int(attributes['label']) ]
                else:
                    attributes['label'] = str(attributes['label'], 'ascii')
    
        best_tree.graph['label'] = f"Confidence score: {c_out.calculated_likelihood}"
        best_tree.graph['labelloc'] = 't'

        if get_cells:
            for i, cell in enumerate(labeled_genotype_matrix.cell_labels):
                cell_id = f'cell: {cell}'
                best_tree.add_node(cell_id, shape = 'box')
                best_tree.add_edge(str(c_out.ids_of_leaves[i]), cell_id)
    
        # Unmarshalling the expected genotype matrix
        expected_matrix = np.ndarray([N, M])
        for i in range(N):
            for j in range(M):
                expected_matrix[i][j] = c_out.gtp_matrix[i][j]
    
        # Unmarshalling the error learning parameters
        el_alphas = [None] * M
        el_gammas = [None] * M
        for i in range(M):
            el_gammas[i] = c_out.el_gammas[i]
            el_alphas[i] = c_out.el_alphas[i]
        el_beta = c_out.el_beta
    
        # Building the output
        best_tree = SASCPhylogeny(best_tree)
        expected_matrix = GenotypeMatrix(
            expected_matrix,
            labeled_genotype_matrix.cell_labels,
            labeled_genotype_matrix.mutation_labels
        )
    
        out = {
            'inferred_tree': best_tree,
            'expected_genotype_matrix': expected_matrix,
            'inferred_alphas': el_alphas,
            'inferred_beta': el_beta,
            'inferred_gammas': el_gammas
        }

        return out
    
    finally:
        # We do the memory cleanup after everything is done to avoid memory leaks.
        free(c_out.el_alphas)
        free(c_out.el_gammas)
        for i in range(N):
            free(c_out.gtp_matrix[i])
        free(c_out.gtp_matrix)
        free(c_out.ids_of_leaves)
        free(arguments.mutation_labels)
        sca.destroy_tree(c_out.best_tree)


cdef _unmarshal_tree(sca.node_t* node, G):
    if (node == NULL):
        return

    # Recursion if the node has children and/or siblings.
    _unmarshal_tree(node.next_sibling, G)
    _unmarshal_tree(node.first_child, G)

    # If the node is a deletion, color it in red.
    if (node.loss == 1):
        G.add_node(str(node.id), label = node.label, fillcolor = 'indianred1', style = 'filled')
    else:
        G.add_node(str(node.id), label = node.label)

    # Add the arc from the parent of the node to the node (if this is not the root node).
    if(node.parent != NULL):
        G.add_edge(str(node.parent.id), str(node.id))


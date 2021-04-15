import numpy as np
import pandas as pd
import networkx as nx
cimport sascpycapi as sca
from libc.limits cimport INT_MAX
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf
from phylo.core.labeledmutationmatrix import LabeledMutationMatrix
from phylo.core.phylogenytree import PhylogenyTree


def compute(
        labeled_mutation_matrix,
        alphas,
        beta,
        k,
        max_deletions = INT_MAX,
        repetitions = 5,
        start_temp = 10**4,
        cooling_rate = 0.01, #this cannot be written as 10**-2, otherwise Cython builds it from 10L and assigns 0 to it.
	cores = 1,
        el_a_variance = 0,
        el_b_variance = 0,
        el_g_variance = 0,
        monoclonal = False,
        gammas = 1,
        get_cells = False
):
    """
    Infers a tree from a matrix representing the mutations for a sample of cells using
    simulated annealing to build the tree. More informations at !!!insert link to SASC here!!!

    Parameters:
    - labeled_mutation_matrix(LabeledMutationMatrix): The matrix from which to infer the phylogeny tree.
    - alphas(float or list(float)): The probability of a false negative for each mutation (that is, for each column in the matrix).
      this can be a single float, in which case the value will be used for all mutations, or a collection of
      floats with an element for each mutation.
    - gammas(float or list(float)), by default 1: The prior loss probability for each mutation. If a single float is specified, then
      it will be used for all mutations, and if nothing is specified the used value will be 1. Otherwise, a 
      collection of floats must be provided with a gamma value for each mutation.
    - beta(float): The probability for a false positive for each mutation. A single value must be provided.
    - k: The k used in the k-Dollo model for tree inference. This is the number of times a mutation can
      undergo a deletion in the tree.
    - max_deletions(int), by default infinity: The maximum amount of deletions globally allowed in the tree. 
    - repetitions(int), by default 5: The number of times the inference process will be repeated. 
      The tree with the best score will be used for the output.
    - monoclonal(bool), by default False: If this is true, then SASC will be forced to infer a monoclonal tree, 
      with exactly one child for the germline.
    - start_temp(float), by default 10**4: The starting temperature for simulated annealing.
    - cooling_rate(float), by default 10**2: The cooling rate for simulated annealing.
    - cores, by default 1: The number of cores used by SASC. 
    - el_a_variance(float), by default 0: The variance of the false negative rates for error learning.
    - el_b_variance(float), by default 0: The variance of the false positive rates for error learning.
    - el_g_variance(float), by default 0: The variance of the prior loss probabilities for error learning.
    - get_cells(bool), by default False: If this is true, the tree in the output will also include the nodes
      for the cells in the sample.

    Returns:
    - inferred_tree (PhylogenyTree): the tree that was inferred by SASC, with its confidence score.
    - expected_genotype_matrix (LabeledMutationMatrix): the input matrix in which the missing information
      on the genotypes of the cells has been inferred.
    - inferred_alphas (list(float)): the inferred false negative rates after error learning.
    - inferred_beta (float): the inferred false positive rate after error learning.
    - inferred_gammas (list(float)): the inferred prior loss probabilities after error learning.
    """


    cdef sca.sasc_in_t args_struct
    cdef sca.sasc_in_t* arguments = &args_struct

    mutations_matrix = labeled_mutation_matrix.matrix()
    arguments.N = len(mutations_matrix)
    arguments.M = len(mutations_matrix[0])    
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
        raise ValueError(f'all alpha values must be within 0 and 1, but at least one in {alphas} is not')
    if any([gamma < 0 or gamma > 1 for gamma in gammas]):
        raise ValueError(f'all gamma values must be within 0 and 1, but at least one in {gammas} is not')

    
    # Some inputs need to be processed and marshalled "manually" before being fed to the c function.
    # Let Cython handle the rest of them by itself.

    
    # Automatic marshalling.
    arguments.beta = beta #;print(beta)
    arguments.el_a_variance = el_a_variance #;print(el_a_variance) 
    arguments.el_b_variance = el_b_variance #;print(el_b_variance)
    arguments.el_g_variance = el_g_variance #;print(el_g_variance)
    if k == 0 or max_deletions == 0:
        k = 0
        max_deletions = 0
    arguments.k = k #;print(k)
    arguments.max_deletions = max_deletions #;print(max_deletions)
    arguments.repetitions = repetitions #;print(repetitions)
    arguments.start_temp = start_temp #;print(start_temp)
    arguments.cooling_rate = cooling_rate #;print(cooling_rate)
    arguments.cores = cores #;print(cores)

    
    # Marshalling of the matrix
    
    arguments.mutations_matrix = <int**>malloc(N*sizeof(int*))
    if arguments.mutations_matrix == NULL:
        raise Exception('out of memory')
    for i in range(N):
        arguments.mutations_matrix[i] = <int*>malloc(M*sizeof(int))
        if arguments.mutations_matrix[i] == NULL:
            raise Exception('out of memory')
        for j in range(M):
            arguments.mutations_matrix[i][j] = mutations_matrix[i][j]

    
    # To work around the string size limitation for SASC, we pass their indexes in the lists in string form to SASC instead of the
    # labels.
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
        raise Exception('out of memory')
    
    for i in range(M):
        arguments.alphas[i] = alphas[i]
        arguments.gammas[i] = gammas[i]

        
    # Python bools must be converted into C integers
    arguments.single_alpha = 1 if single_alpha else 0
    arguments.single_gamma = 1 if single_gamma else 0
    arguments.monoclonal = 1 if monoclonal else 0

    cdef sca.sasc_out_t out_struct
    cdef sca.sasc_out_t* c_out = &out_struct

    c_out.gtp_matrix = <int**> malloc(N * sizeof(int*))
    if c_out.gtp_matrix == NULL:
        raise Exception('out of memory')
    for i in range(N):
        c_out.gtp_matrix[i] = <int*> malloc(M * sizeof(int))
        if c_out.gtp_matrix[i] == NULL:
            raise Exception('out of memory')
    
    c_out.ids_of_leaves = <int*> malloc(N * sizeof(int))
    c_out.el_alphas = <double*> malloc(M * sizeof(double))
    c_out.el_gammas = <double*> malloc(M * sizeof(double))
    if (
        c_out.el_alphas == NULL
        or c_out.el_gammas == NULL
        or c_out.ids_of_leaves == NULL
    ):
        raise Exception('out of memory')

    cdef int comp_result
    cdef sca.node_t* root

    try:
        
        # THE C A L L
        with nogil:
            comp_result = sca.compute(arguments, c_out)
    
        # And now for the unmarshalling. We'll output a tuple with the tree as a networkx graph, the matrix as
        # a numpy array, and the rest of the values as simple ints/doubles/strings.
    
        # Unmarshalling of the tree.
        best_tree = nx.DiGraph()
        unmarshal_tree(c_out.best_tree, best_tree)

        # The labels for the tree are now indexes to the original mutation labels, aside from special labels that are
        # added artificially to the tree during the execution of SASC and only need to be converted
        # to a string.

        # We assume that every label that cannot be converted to an integer needs to be left as-is
        def converts_to_index(label):
            try:
                return int(label) >= 0
            except ValueError:
                return False
    
        for [node, attributes] in best_tree.nodes(data = True):
            if 'label' in attributes:
                if converts_to_index(attributes['label']):
                    attributes['label'] = labeled_mutation_matrix.mutation_labels[ int(attributes['label']) ]
                else:
                    attributes['label'] = str(attributes['label'], 'ascii')
    
        best_tree.graph['label'] = f"Confidence score: {c_out.calculated_likelihood}"
        best_tree.graph['labelloc'] = 't'

        if get_cells:
            for i, cell in enumerate(labeled_mutation_matrix.cell_labels):
                cell_id = f'cell: {cell}'
                best_tree.add_node(cell_id, shape = 'box')
                best_tree.add_edge(str(c_out.ids_of_leaves[i]), cell_id)

        sca.destroy_tree(c_out.best_tree)
    
        # Unmarshalling of the expected genotype matrix
        expected_matrix = np.ndarray([N, M])
        for i in range(N):
            for j in range(M):
                expected_matrix[i][j] = c_out.gtp_matrix[i][j]
    
        # Unmarshalling of the error learning parameters
        el_alphas = [None] * M
        el_gammas = [None] * M
        for i in range(M):
            el_gammas[i] = c_out.el_gammas[i]
            el_alphas[i] = c_out.el_alphas[i]
    
        el_beta = c_out.el_beta
    
        # Building the output
        best_tree = PhylogenyTree(best_tree)
        expected_matrix = LabeledMutationMatrix(expected_matrix, labeled_mutation_matrix.cell_labels, labeled_mutation_matrix.mutation_labels)
    
        out = {
            'inferred_tree': best_tree,
            'expected_genotype_matrix': expected_matrix,
            'inferred_alphas': el_alphas,
            'inferred_beta': el_beta,
            'inferred_gammas': el_gammas
        }
    
        # Cleanup and return
        return out
    
    finally:
        free(c_out.el_alphas)
        free(c_out.el_gammas)
        for i in range(N):
            free(c_out.gtp_matrix[i])
        free(c_out.gtp_matrix)
        free(c_out.ids_of_leaves)
        free(arguments.mutation_labels)


cdef unmarshal_tree(sca.node_t* node, G):

    if (node == NULL):
        return;

    unmarshal_tree(node.next_sibling, G)
    unmarshal_tree(node.first_child, G)

    # If the node is a deletion, color it in red.
    if (node.loss == 1):
        G.add_node(str(node.id), label = node.label, color = 'indianred1', style = 'filled')
    else:
        G.add_node(str(node.id), label = node.label)

    # Add the arc from the parent of the node to the node (if this is not the root node).
    if(node.parent != NULL):
        G.add_edge(str(node.parent.id), str(node.id))


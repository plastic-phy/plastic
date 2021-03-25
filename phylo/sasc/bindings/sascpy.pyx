import numpy as np
import pandas as pd
import networkx as nx
cimport sascpycapi as sca
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf
from phylo.core.labeledmutationmatrix import LabeledMutationMatrix
from phylo.core.phylogenytree import PhylogenyTree

def compute(
        labeled_mutation_matrix,
        alphas,
        beta,
        gammas,
        el_a_variance,
        el_b_variance,
        el_g_variance,
        k,
        max_deletions,
        repetitions,
        force_monoclonal,
        start_temp,
        cooling_rate,
        cores,
        get_leaves = False
):
    
    cdef sca.sasc_in_t* arguments = <sca.sasc_in_t*>malloc(sizeof(sca.sasc_in_t))
    # Some inputs need to be processed and marshalled "manually" before being fed to the c function.
    # Let Cython handle the rest of them by itself.
    
    # Marshalling of the matrix, assuming that the rows represent cells and the columns represent
    # mutations.

    mutations_matrix = labeled_mutation_matrix.matrix()
    
    arguments.N = len(mutations_matrix)
    arguments.M = len(mutations_matrix[0])
    
    cdef int N = arguments.N
    cdef int M = arguments.M
    
    arguments.mutations_matrix = <int**>malloc(N*sizeof(int*))
    for i in range(N):
        arguments.mutations_matrix[i] = <int*>malloc(M*sizeof(int))
        for j in range(M):
            arguments.mutations_matrix[i][j] = mutations_matrix[i][j]
    
    # Marshalling of the cell and mutation labels. We convert the label strings into bytes assuming ASCII format.
    
    arguments.cell_labels = <char**>malloc(N * sizeof(char*))
    arguments.mutation_labels = <char**>malloc(M * sizeof(char*))
    
    cell_labels_bytes = [bytes(lb, 'ascii') for lb in labeled_mutation_matrix.cell_labels] 
    mutation_labels_bytes = [bytes(lb, 'ascii') for lb in labeled_mutation_matrix.mutation_labels]
    
    for i in range(N):
        arguments.cell_labels[i] = cell_labels_bytes[i]
    for i in range(M):
        arguments.mutation_labels[i] = mutation_labels_bytes[i]
    
    
    
    # Arrays with error parameters must be allocated and filled
    arguments.alphas = <double*>malloc(M * sizeof(double))
    arguments.gammas = <double*>malloc(M * sizeof(double))

    single_alpha = isinstance(alphas, float)
    single_gamma = isinstance(gammas, float)

    # Python bools must be converted into C integers
    arguments.single_alpha = 1 if single_alpha else 0
    arguments.single_gamma = 1 if single_gamma else 0
    arguments.force_monoclonal = 1 if force_monoclonal else 0

    if single_alpha:
        alphas = [alphas] * M
    if single_gamma:
        gammas = [gammas] * M

    if len(alphas) != M:
        raise ValueError(f'multiple alphas are specified in {alphas}, but they are more or less than the number of mutations.')
    if len(gammas) != M:
        raise ValueError(f'multiple gammas are specified in {gammas}, but they are more or less than the number of mutations.')
    for i in range(M):
        arguments.alphas[i] = alphas[i]
        arguments.gammas[i] = gammas[i]
    
    # Automatic marshalling.
    arguments.beta = beta
    arguments.el_a_variance = el_a_variance
    arguments.el_b_variance = el_b_variance
    arguments.el_g_variance = el_g_variance
    arguments.k = k
    arguments.max_deletions = max_deletions
    arguments.repetitions = repetitions
    arguments.start_temp = start_temp
    arguments.cooling_rate = cooling_rate
    arguments.cores = cores
    
    cdef sca.sasc_out_t out_struct
    cdef sca.sasc_out_t* c_out = &out_struct

    c_out.gtp_matrix = <int**> malloc(N * sizeof(int*))
    for i in range(N):
        c_out.gtp_matrix[i] = <int*> malloc(M * sizeof(int))

    c_out.ids_of_leaves = <int*> malloc(N * sizeof(int));
    c_out.el_alphas = <double*> malloc(M * sizeof(double));
    c_out.el_gammas = <double*> malloc(M * sizeof(double));

    cdef int comp_result
    # THE C A L L
    with nogil:
        comp_result = sca.compute(arguments, c_out)
    print(comp_result)
    
    # And now for the unmarshalling. We'll output a tuple with the tree as a networkx graph, the matrix as
    # a numpy array, and the rest of the values as simple ints/doubles/strings.
    
    # Unmarshalling of the tree.
    cdef sca.node_t* root = c_out.best_tree
    best_tree = nx.DiGraph()
    
    unmarshal_tree(root, best_tree)

    best_tree.graph['label'] = f"Confidence score: {c_out.calculated_likelihood}"
    best_tree.graph['labelloc'] = 't'

    if get_leaves:
        for i, cell in enumerate(cell_labels):
            cell_id = f'cell: {cell}'
            best_tree.add_node(cell_id, shape = 'box')
            best_tree.add_edge(str(c_out.ids_of_leaves[i]), cell_id)

    sca.destroy_tree(c_out.best_tree)
    free(c_out.ids_of_leaves)
    
    # Unmarshalling of the expected genotype matrix
    expected_matrix = np.ndarray([N, M])
    for i in range(N):
        for j in range(M):
            expected_matrix[i][j] = c_out.gtp_matrix[i][j]
        free(c_out.gtp_matrix[i])
    free(c_out.gtp_matrix)
    
    # Unmarshalling of the error learning parameters
    el_alphas = [None] * M
    el_gammas = [None] * M
    for i in range(M):
        el_gammas[i] = c_out.el_gammas[i]
        el_alphas[i] = c_out.el_alphas[i]

    free(c_out.el_alphas)
    free(c_out.el_gammas)
    
    el_beta = c_out.el_beta
    
    # Building the output
    best_tree = PhylogenyTree(best_tree)
    expected_matrix = LabeledMutationsMatrix(expected_matrix, mutations_matrix.cell_labels, mutations_matrix.mutation_labels)
    
    out = best_tree, expected_matrix, el_alphas, el_beta, el_gammas
    
    # Cleanup and return
    return out


cdef unmarshal_tree(sca.node_t* node, G):

    if (node == NULL):
        return;

    unmarshal_tree(node.next_sibling, G)
    unmarshal_tree(node.first_child, G)

    # If the node is a deletion, color it in red.
    if (node.loss == 1):
        G.add_node(str(node.id), color = 'indianred1', style = 'filled', label = str(node.label, 'ascii'))
    else:
        G.add_node(str(node.id), label = str(node.label, 'ascii'))

    # Add the arc from the parent of the node to the node (if this is not the root node).
    if(node.parent != NULL):
        G.add_edge(str(node.parent.id), str(node.id))


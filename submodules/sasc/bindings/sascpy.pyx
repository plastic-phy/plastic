import numpy as np
import pandas as pd
import networkx as nx
cimport sascpycapi as sca
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

def compute(
        mutations_matrix,
        mutation_labels,
        cell_labels,
        alphas,
        beta,
        gammas,
        single_alpha,
        single_gamma,
        el_a_variance,
        el_b_variance,
        el_g_variance,
        k,
        max_deletions,
        repetitions,
        force_monoclonal,
        start_temp,
        cooling_rate,
        cores
):
    
    cdef sca.sasc_in_t* arguments = <sca.sasc_in_t*>malloc(sizeof(sca.sasc_in_t))
    # Some inputs need to be processed and marshalled "manually" before being fed to the c function.
    # Let Cython handle the rest of them by itself.
    
    # Marshalling of the matrix, assuming row-major list representation (that is, each column is a contiguous
    # sublist of the original one.
    
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
    
    cell_labels_bytes = [bytes(lb, 'ascii') for lb in cell_labels] 
    mutation_labels_bytes = [bytes(lb, 'ascii') for lb in mutation_labels]
    
    for ix in range(N):
        arguments.cell_labels[ix] = cell_labels_bytes[ix]
    for ix in range(M):
        arguments.mutation_labels[ix] = mutation_labels_bytes[ix]
    
    # Python bools must be converted into C integers
    arguments.single_alpha = 1 if single_alpha else 0
    arguments.single_gamma = 1 if single_gamma else 0
    arguments.force_monoclonal = 1 if force_monoclonal else 0
    
    # Arrays with error parameters must be allocated and filled
    arguments.alphas = <double*>malloc(M * sizeof(double))
    arguments.gammas = <double*>malloc(M * sizeof(double))
    
    for ix in range(M):
        arguments.alphas[ix] = alphas[ix]
        arguments.gammas[ix] = gammas[ix]
    
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
    
    cdef sca.sasc_out_t* c_out
    
    # THE C A L L.
    c_out = sca.compute(arguments)
    
    # And now for the unmarshalling. We'll output a tuple with the tree as a networkx graph, the matrix as
    # a numpy array, and the rest of the values as simple ints/doubles
    
    
    
    # Unmarshalling of the tree.
    cdef sca.node_t* root = c_out.best_tree
    best_tree = nx.DiGraph()
    best_tree.add_node(root.id, label = str(root.label), loss = root.loss, mutation_index = root.mut_index)
    unmarshal_tree(root, best_tree)

    sca.destroy_tree(c_out.best_tree)
    
    #Unmarshalling of the expected matrix
    expected_matrix = np.ndarray([N, M])
    for i in range(N):
        for j in range(M):
            expected_matrix[i][j] = c_out.gtp_matrix[i][j]
        free(c_out.gtp_matrix[i])
    free(c_out.gtp_matrix)
    print(expected_matrix)
    
    # Unmarshalling of the error learning parameters
    el_alphas = [None] * M
    el_gammas = [None] * N
    for ix in range(M):
        el_gammas[ix] = c_out.el_gammas[ix]
        el_alphas[ix] = c_out.el_alphas[ix]
    
    el_beta = c_out.el_beta
    
    # Build the output 
    out = best_tree, c_out.calculated_likelihood, expected_matrix, el_alphas, el_beta, el_gammas
    
    # Cleanup and return
    free(c_out)
    return out


cdef unmarshal_tree(sca.node_t* root, G):
    
    cdef sca.node_t* curr_child = root.first_child
    while curr_child != NULL:
        G.add_node(curr_child.id, label = str(curr_child.label), loss = curr_child.loss, mutation_index = curr_child.mut_index)
        G.add_edge(root.id, curr_child.id)
        unmarshal_tree(curr_child, G)
        curr_child = curr_child.next_sibling;


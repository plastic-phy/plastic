import numpy as np
import pandas as pd
import networkx as nx
cimport sasc-compute-api as sca
from libc.stdlib cimport malloc, free

def compute(
        mutations_matrix,
        mutation_labels,
        cell_labels,
        alphas,
        betas,
        gammas,
        single_alpha,
        single_gamma,
        el_a_variances,
        el_b_variance,
        el_g_variances,
        k,
        max_deletions,
        repetitions,
        force_monoclonal,
        start_temp,
        cooling_rates,
        cores
):
    
    cdef sca.sasc_in_t* arguments = <sca.sasc_in_t*>malloc(sizeof(sca.sasc_in_t))
    # Some inputs need to be processed and marshalled "manually" before being fed to the c function.
    # Let Cython handle the rest of them by itself.

    # Marshalling of the matrix, assuming row-major list representation (that is, each column is a contiguous
    # sublist of the original one.

    arguments.N = len(mutations_matrix)
    arguments.M = len(mutations_matrix[0])
    arguments.mutations_matrix = <int**>malloc(N*sizeof(int*))
    for i in range(N):
        arguments.mutations_matrix[i] = <int*>malloc(M*sizeof(int))
        for j in range(M):
            arguments.mutations_matrix[i][j] = mutations_matrix[i][j]
    
    # Marshalling of the cell and mutation labels. We convert the label strings into bytes assuming ASCII format.
    
    cdef char* arglabcells[N]
    cdef char* arglabmuts[M]

    for ix in range(N):
        arglabcells[ix] = bytes(cell_labels[ix], 'ascii')
    for ix in range(M):
        arglabmuts[ix] = bytes(mutation_labels[ix], 'ascii')

    arguments.cell_labels = <char**>arglabcells
    arguments.mutation_labels = <char**>arglabmuts
    
    # Python bools must be converted into C integers
    arguments.single_alpha = 1 if single_alpha else 0
    arguments.single_gamma = 1 if single_gamma else 0

    # Automatic marshalling.
    arguments.alphas = alphas
    arguments.betas = betas
    arguments.gammas = gammas
    arguments.el_a_variances = el_a_variances
    arguments.el_b_variances = el_b_variance
    arguments.el_g_variances = el_g_variances
    arguments.k = k
    arguments.max_deletions = max_deletions
    arguments.repetitions = repetitions
    arguments.force_monoclonal = force_monoclonal
    arguments.start_temp = start_temp
    arguments.cooling_rates = cooling_rates
    arguments.cores = cores
    
    # THE C A L L. NoGIL tells Cython it isn't doing anything on Python objects.
    with nogil:
        cdef sca.sasc_out_t* c_out = sca.compute(arguments)

    # And now for the unmarshalling. We'll output a tuple with the tree as a networkx graph, the matrix as
    # a numpy array, and the rest of the values as simple ints/doubles

    # Unmarshalling of the tree. 
        #TODO
    sca.destroy_tree(out.best_tree)

    #Unmarshalling of the expected matrix
    expected_matrix = np.ndarray([N, M])
    for i in range(N):
        for j in range(M):
            gtp_matrix[i][j] = c_out.gtp_matrix[i][j]
        free(c_out.gtp_matrix[i])
    free(c_out.gtp_matrix)

    # Unmarshalling of the error learning parameters
    el_alphas = [None] * M
    el_gammas = [None] * N
    for ix in range(M):
        el_gammas[ix] = c_out.el_gammas[ix]
        el_alphas[ix] = c_out.el_alphas[ix]

    el_beta = out.el_beta
    
    # Build the output 
    out = best_tree, out.calculated_likelihood, expected_matrix, el_alphas, el_beta, el_gammas 

    # Cleanup and return
    free(c_out)
    return out


    

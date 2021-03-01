#ifndef SASC_COMPUTE_H
#define SASC_COMPUTE_H

/* Ideally this is all we need to expose to python, and we can let the python side do I/O, input validity checks
   and so forth. */

typedef struct sasc_output {
    node_t* best_tree,
    double calculated_likelihood,
    int** gtp_matrix,
    double el_alphas,
    double el_beta,
    double el_gammas
} sasc_out_t;

// This struct might correspond to a couple different ones on the python side.
typedef struct sasc_input {
/*-------SASC PARAMETERS-------*/
    int K,
    double* alphas,     // If a single alpha was specified for the input, this will be filled with that value 
    double beta,
    double* gammas,     // Similar behaviour for single gammas
    int max_deletions,
    int repetitions,
    int force_monoclonal,
    double start_temp,
    double cooling_rate,
    int cores,
/*-------ERROR LEARNING PARAMETERS-------*/
    double el_a_variance,
    double el_b_variance,
    double el_g_variance,
/*-------MUTATIONS MATRIX-------*/
    int** mutations_matrix,
    int N,
    int M,
    buffer255_t** mutation_labels,
    buffer255_t** cell_labels  
} sasc_in_t;

/* Tis might simplify allocation during marshalling (less mallocs+frees), but it might be better to leave it as a pointer.
   Truth be told, I don't know much C at all, so I'd like feedback on this. What I do know is that we probably need to check that 
   we don't go beyond bounds while filling this during marshalling, as C doesn't stop us from doing that.*/
typedef struct fxdsizebuf {
    char b[255]
} buffer255_t;

sasc_out_t* compute(sasc_in_t* arguments);

#endif

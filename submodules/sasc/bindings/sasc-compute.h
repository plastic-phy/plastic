#ifndef SASC_COMPUTE_H
#define SASC_COMPUTE_H

/* Ideally this is all we need to expose to python, and we can let the python side do I/O, input validity checks
   and so forth. */

typedef struct Node node_t;

typedef struct sasc_output {
    node_t* best_tree;
    double calculated_likelihood;
    int** gtp_matrix;
    double* el_alphas;
    double el_beta;
    double* el_gammas;
} sasc_out_t;

// This struct might correspond to a couple different ones on the python side.
typedef struct sasc_input {
/*-------SASC PARAMETERS-------*/
    int k;
    int max_deletions;
    int repetitions;
    int force_monoclonal;
    double start_temp;
    double cooling_rate;
    int cores;
/*-------MUTATIONS MATRIX-------*/
    int** mutations_matrix;
    int N;
    int M;
    char **mutation_labels;
    char **cell_labels;
/*-------ERROR PARAMETERS-------*/
    double* alphas; int single_alpha;    // If a single alpha was specified for the input, this will be filled with that value 
    double beta;
    double* gammas; int single_gamma;     // Similar behaviour for single gammas
/*-------ERROR LEARNING PARAMETERS-------*/
    double el_a_variance;
    double el_b_variance;
    double el_g_variance;
} sasc_in_t;


sasc_out_t* compute(sasc_in_t* arguments); 

#endif

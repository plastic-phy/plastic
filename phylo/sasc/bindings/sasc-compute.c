/*
    MIT License

    Copyright (c) 2017-2019 Simone Ciccolella

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mt19937ar.h"
#include <time.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "tree.h"
#include "utils.h"
#include "sastep.h"
#include "vector.h"
#include "sasc-compute.h"

#ifdef NDEBUG
#include <assert.h>
#else
#define assert(ignore)((void) 0)
#endif

/* 
    Business logic extracted from the SASC main; parameters are expressed
    as a data structure rather than in files, and the output is a data 
    structure that contains the best philogeny tree found through SA, its likelihood,
    the expected genotype matrix and the error rates acquired through error learning.
 */

int MAX_ID_TREE = 0;

int compute (sasc_in_t* arguments, sasc_out_t* out)
{
    int N = arguments->N;
    int M = arguments->M;
    int K = arguments->k;

    int SINGLE_ALPHA = arguments->single_alpha;
    int SINGLE_GAMMA = arguments->single_gamma;
    
    int MAX_LOSSES = arguments->max_deletions;
    int MONOCLONAL = arguments->monoclonal;

    printf("Starting SASC.\n");

    // Set random seed
    srand((unsigned)time(NULL));
    
    // Read cell and mutation labels from input
    char MUT_NAMES[M][255];
    
    for (int i = 0; i < M; i++) {
      strcpy(MUT_NAMES[i], arguments->mutation_labels[i]);
    }

    // Get error parameters
    double MULTI_ALPHAS[M];
    double MULTI_GAMMAS[M];
    
    for (int i = 0; i < M; i++) {
        MULTI_ALPHAS[i] = arguments->alphas[i];
        MULTI_GAMMAS[i] = arguments->gammas[i];
    }

    double BETA = arguments->beta;

    // Itialize MT19937
    unsigned long init[10], length=10;
    FILE *f;
    f = fopen("/dev/urandom", "r");

    for (int i = 0; i < length; i++){
        unsigned long randval;
        fread(&randval, sizeof(randval), 1, f);
        init[i] = randval;
    }
    fclose(f);

    init_by_array(init, length);

    // Initialize the input matrix
    int** INPUT_MATRIX = malloc(N * (sizeof *INPUT_MATRIX));
    if (INPUT_MATRIX == NULL) {
        printf("out of memory");
	exit(1);
    }
    for (int i = 0; i < N; i++) {
        INPUT_MATRIX[i] = malloc(M * (sizeof *INPUT_MATRIX[i]));
	if (INPUT_MATRIX[i] == NULL) {
	    printf("out of memory");
	    exit(1);
	}
	for (int j = 0; j < M; j++) {
	    INPUT_MATRIX[i][j] = arguments->mutations_matrix[i][j];
	}
    }

    double START_TEMP = arguments->start_temp;
    double COOLING_RATE = arguments->cooling_rate;
    double MIN_TEMP = 0.001;

    int REPETITIONS = arguments->repetitions;
    node_t *best_tree = NULL;
    double best_loglike = -DBL_MAX;
    int best_sigma[N]; for (int i = 0; i < N; i++) { best_sigma[i] = 0; }
    vector best_tree_vec;
    vector_init(&best_tree_vec);
    vector best_losses_vec;
    vector_init(&best_losses_vec);

    double a_mu[M];
    double a_xs[M];
    double g_mu[M];
    double g_xs[M];

    for (int i = 0; i < M; i++) {
        a_mu[i] = MULTI_ALPHAS[i];
        a_xs[i] = MULTI_ALPHAS[i];
        g_mu[i] = MULTI_GAMMAS[i];
        g_xs[i] = MULTI_GAMMAS[i];
    }
    elpar_t* el_params = set_el_params(SINGLE_ALPHA, M,
				       MULTI_ALPHAS, a_mu, arguments->el_a_variance, a_xs,
                                       &BETA, BETA, arguments->el_b_variance,
                                       MULTI_GAMMAS, g_mu, arguments->el_g_variance, g_xs,
                                       SINGLE_GAMMA);

    // Generate Cj
    int Cj[M];
    for (int i = 0; i < M; i++) {
        Cj[i] = 0;
    }

    for (int r = 0; r < REPETITIONS; r++) {
        printf("Iteration: %d\n", r+1);
        
        // Generate a RANDOM BTREE
        vector ml_tree_vec;
        vector ml_losses_vec;
        vector_init(&ml_tree_vec);
        vector_init(&ml_losses_vec);

        node_t *root = node_new("germline", -1, 0);
        vector_add(&ml_tree_vec, root);

        int rantree[M];
        for (int i = 0; i < M; i++) {
            rantree[i] = i;
        }
        shuffle(rantree, M);

        if (MONOCLONAL == 0) {
            int app_node = 0;
            for (int i = 0; i < M; i++) {
                node_t *cnode1 = node_new(MUT_NAMES[rantree[i]], rantree[i], vector_total(&ml_tree_vec));
                vector_add(&ml_tree_vec, cnode1);
                node_append(vector_get(&ml_tree_vec, app_node), cnode1);
                i += 1;

                if (i < M){
                    node_t *cnode2 = node_new(MUT_NAMES[rantree[i]], rantree[i], vector_total(&ml_tree_vec));
                    vector_add(&ml_tree_vec, cnode2);
                    node_append(vector_get(&ml_tree_vec, app_node), cnode2);
                }
                app_node += 1;
            }


        } else {
            node_t *first_clone = node_new(MUT_NAMES[rantree[0]], rantree[0], vector_total(&ml_tree_vec));
            vector_add(&ml_tree_vec, first_clone);
            node_append(vector_get(&ml_tree_vec, 0), first_clone);

            int app_node = 1;
            for (int i = 1; i < M; i++) {
                node_t *cnode1 = node_new(MUT_NAMES[rantree[i]], rantree[i], vector_total(&ml_tree_vec));
                vector_add(&ml_tree_vec, cnode1);
                node_append(vector_get(&ml_tree_vec, app_node), cnode1);
                i += 1;

                if (i < M){
                    node_t *cnode2 = node_new(MUT_NAMES[rantree[i]], rantree[i], vector_total(&ml_tree_vec));
                    vector_add(&ml_tree_vec, cnode2);
                    node_append(vector_get(&ml_tree_vec, app_node), cnode2);
                }
                app_node += 1;
            }
        } 

        // Generate SIGMA
        int SIGMA[N];
        for (int i = 0; i < N; i++) {
            SIGMA[i] = 0;
        }

        // get log-likelihood
        double lh = greedy_tree_loglikelihood(root, ml_tree_vec, SIGMA, INPUT_MATRIX, N, M, MULTI_ALPHAS, BETA, MULTI_GAMMAS, Cj, arguments->cores);
        // for (int i = 0; i < N; i++) { printf("%d ", SIGMA[i]); } printf("\n");
        // printf("Start log-like: %lf\n", lh);
        node_t *ml_tree = anneal(root, ml_tree_vec, N, M, K, MULTI_ALPHAS, BETA, INPUT_MATRIX, START_TEMP, COOLING_RATE,
                MIN_TEMP, MAX_LOSSES, el_params, MULTI_GAMMAS, Cj, MONOCLONAL, arguments->cores);
        
        vector_free(&ml_tree_vec);
        vector_init(&ml_tree_vec);
        vector_free(&ml_losses_vec);
        vector_init(&ml_losses_vec);
        
        ml_tree = treecpy(ml_tree, &ml_tree_vec, &ml_losses_vec, N);

//        for (int i = 0; i < M; i++) { printf("%d ", Cj[i]); } printf("\n");

        double current_lh = greedy_tree_loglikelihood(ml_tree, ml_tree_vec, SIGMA, INPUT_MATRIX, N, M, MULTI_ALPHAS, BETA, MULTI_GAMMAS, Cj, arguments->cores);
//        printf("Maximum log-likelihood: %lf\n", current_lh);

        if (current_lh > best_loglike) {
            if (best_tree != NULL)
                destroy_tree(best_tree);

            vector_free(&best_tree_vec);
            vector_free(&best_losses_vec);

            vector_init(&best_tree_vec);
            vector_init(&best_losses_vec);

            best_tree = treecpy(ml_tree, &best_tree_vec, &best_losses_vec, N);
            best_loglike = current_lh;
        }

        vector_free(&ml_tree_vec);
        destroy_tree(root);
        destroy_tree(ml_tree);
    }

    double best_calculated_likelihood = greedy_tree_loglikelihood(best_tree, best_tree_vec, best_sigma, INPUT_MATRIX, N, M, MULTI_ALPHAS, BETA, MULTI_GAMMAS, Cj, arguments->cores);

    
    int gtpo[M];
    

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) { gtpo[j] = 0; }

	node_t *node = vector_get(&best_tree_vec, best_sigma[i]);
	assert(node != NULL);
	get_genotype_profile(vector_get(&best_tree_vec, best_sigma[i]), gtpo);

	for (int j = 0; j < M; j++) {
	  out->gtp_matrix[i][j] = gtpo[j]; 
	}
    }
    
    out->best_tree = best_tree;
    out->calculated_likelihood = best_calculated_likelihood;

    for (int i = 0; i < N; i++) {
        node_t *leaf = vector_get(&best_tree_vec, best_sigma[i]);
	out->ids_of_leaves[i] = leaf->id;
    }

    for(int i = 0; i < M; i++) {
        out->el_alphas[i] = a_xs[i];
        out->el_gammas[i] = g_xs[i];
    }
    out->el_beta = el_params->b_x;
    
    for (int i = 0; i < N; i++) { free(INPUT_MATRIX[i]); }
    free(INPUT_MATRIX);

    return 0;
}

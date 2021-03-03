cdef extern from "tree.h":
    
    ctypedef struct Node node_t:
        int id
        char label[255]
        int loss
        int mut_index
        struct Node *first_child
        struct Node *next_sibling
        struct Node *previous_sibling
        struct Node *parent

    void destroy_tree(node_t* node)
        
cdef extern from "sasc-compute.h":
    
    ctypedef struct sasc_input sasc_in_t:
        int k
        int max_deletions
        int repetitions
        int force_monoclonal
        double start_temp
        double cooling_rate
        int cores
        int** mutations_matrix
        int N
        int M
        char **mutation_labels
        char **cell_labels
	double* alphas
	int single_alpha
        double beta
        double* gammas
	int single_gamma
	double el_a_variance
        double el_b_variance
        double el_g_variance
    
    ctypedef struct sasc_output sasc_out_t:
        node_t* best_tree
        double calculated_likelihood
        int** gtp_matrix
        double el_alphas
        double el_beta
        double el_gammas
    
    sasc_out_t* compute(sasc_int_t* arguments)

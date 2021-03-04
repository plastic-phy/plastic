cdef extern from "tree.h":
    
    struct Node:
        int id
        char label[255]
        int loss
        int mut_index
        Node *first_child
        Node *next_sibling
        Node *previous_sibling
        Node *parent
    
    ctypedef Node node_t

    void destroy_tree(node_t* node)
        
cdef extern from "sasc-compute.h":
    
    struct sasc_input:
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

    ctypedef sasc_input sasc_in_t
    
    struct sasc_output:
        node_t* best_tree
        double calculated_likelihood
        int** gtp_matrix
        double* el_alphas
        double el_beta
        double* el_gammas

    ctypedef sasc_output sasc_out_t
    
    sasc_out_t* compute(sasc_in_t* arguments)

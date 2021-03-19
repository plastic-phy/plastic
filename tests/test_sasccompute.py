from sasc import compute
from core.labeledmutationmatrix import LabeledMutationMatrix

def test_dummy_input():

    inp = LabeledMutationMatrix.from_files('matrepr_files/clustered_example')
    mat = inp.matrix()
    M = len(mat[0])
    
    out = compute(
        mutations_matrix = mat,
        mutation_labels = inp.mutation_labels,
        cell_labels = inp.cell_labels,
        alphas = [0.1] * M,
        beta = 0.00001,
        gammas = [0.2] * M,
        single_alpha = True,
        single_gamma = True,
        el_a_variance = 0.01,
        el_b_variance = 0.01,
        el_g_variance = 0.01,
        k = 1,
        max_deletions = 2,
        repetitions = 5,
        force_monoclonal = False,
        start_temp = 10**4,
        cooling_rate =  10**-1,
        cores = 8
    )

    print([edge for edge in out[0].edges])
    print(out[1])
    print(out[2])
    print(out[3])
    print(out[4])
    print(out[5])


if __name__ == '__main__':
    test_dummy_input()

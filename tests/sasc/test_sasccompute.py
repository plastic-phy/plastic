from phylo.sasc import compute
from phylo.core.labeledmutationmatrix import LabeledMutationMatrix

def test_dummy_input(get_cells = False):

    inp = LabeledMutationMatrix.from_files('matrepr_files/clustered_example')
    
    out = compute(
        labeled_mutation_matrix = inp,
        alphas = 0.1,
        beta = 0.00001,
        gammas = 0.2,
        el_a_variance = 0.01,
        el_b_variance = 0.01,
        el_g_variance = 0.01,
        k = 1,
        max_deletions = 10,
        repetitions = 5,
        force_monoclonal = False,
        start_temp = 10**4,
        cooling_rate =  10**-2,
        cores = 12,
        get_cells = get_cells
    )

    [print(node) for node in out[0].as_digraph().nodes(data = True)]
    [print(edge) for edge in out[0].as_digraph().edges]
    [print(line) for line in out[1].matrix()]
    print(out[2])
    print(out[3])
    print(out[4])


if __name__ == '__main__':
    test_dummy_input()

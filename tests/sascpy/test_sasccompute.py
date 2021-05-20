from phylo.sasc import infer_tree, GenotypeMatrix

def test_dummy_input(get_cells = False):

    parsed = GenotypeMatrix.from_files('matrepr_files/clustered_example')
    with_str_labels = GenotypeMatrix(
        genotype_matrix = parsed.matrix(),
        mutation_labels = [f'm_{lb}' for lb in parsed.mutation_labels]
    )
    
    out = infer_tree(
        labeled_genotype_matrix = with_str_labels,
        alphas = [0.1] * len(with_str_labels.mutations()),
        beta = 0.00001,
        gammas = 0.2,
        el_a_variance = 0.01,
        el_b_variance = 0.01,
        el_g_variance = 0.01,
        k = 1,
        max_deletions = 10,
        repetitions = 5,
        monoclonal = False,
        start_temp = 10**4,
        cooling_rate =  10**-2,
        cores = 12,
        get_cells = get_cells
    )

    [print(node) for node in out['inferred_tree'].as_digraph().nodes(data = True)]
    [print(edge) for edge in out['inferred_tree'].as_digraph().edges]
    [print(line) for line in out['expected_genotype_matrix'].matrix()]
    print(out['inferred_alphas'])
    print(out['inferred_beta'])
    print(out['inferred_gammas'])

if __name__ == '__main__':
    test_dummy_input()

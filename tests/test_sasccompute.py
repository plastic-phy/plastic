from sasc import compute

def test_dummy_input():
    out = compute(
        mutations_matrix = [
            [1, 0, 0],
            [1, 2, 0],
            [1, 0, 1]
        ],
        mutation_labels = ['m1', 'm2', 'm3'],
        cell_labels = ['c1', 'c2', 'c3'],
        alphas = [0.1, 0.1, 0.1],
        beta = 0.00001,
        gammas = [0.2, 0.2, 0.2],
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
        cooling_rate =  10**-2,
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

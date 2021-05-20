from phylo import sasc as sc, celluloid as cd, mp3
from multiprocessing import cpu_count

# Let's parse a matrix and cluster it to different degrees
matrix = cd.GenotypeMatrix.from_files('example_in')
print('starting clustering process')
clustered1 = cd.cluster_mutations(matrix, k=100, number_of_iterations=1, verbose=True)
clustered1.to_files('clust/clustered1', mutations_file='mut_clusters')
print('starting SASC')

out1 = sc.infer_tree(
    clustered1.with_automatic_mutation_labels(),
    alphas=0.1,
    beta=0.001,
    gammas=0.2,
    k=1,
    max_deletions=200000,
    repetitions=3,
    monoclonal=False,
    start_temp=10 ** 4,
    cooling_rate=10 ** -2,
    cores=cpu_count()
)

out2 = sc.infer_tree(
    clustered1.with_automatic_mutation_labels(),
    alphas=0.1,
    beta=0.001,
    gammas=0.2,
    k=1,
    max_deletions=200000,
    repetitions=3,
    monoclonal=False,
    start_temp=10 ** 4,
    cooling_rate=10 ** -2,
    cores=cpu_count(),
    get_cells=True
)

out1['inferred_tree'].draw_to_file('trees/tree1.pdf')
out2['inferred_tree'].to_file('trees/tree2.gv')
out2['inferred_tree'].with_visualization_features(collapse_simple_paths=True) \
    .draw_to_file('trees/tree2.png', show_color=False)

similarity = mp3.tree_similarity(
    out1['inferred_tree'],
    out2['inferred_tree'],
    ignore_unlabeled_nodes=True,
    cores=cpu_count()
)

print(similarity)

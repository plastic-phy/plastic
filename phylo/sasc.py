"""
For more info about this module's functionalities:
- help(sasc.SASCPhylogeny)
- help(sasc.GenotypeMatrix)
- help(sasc.infer_tree)

---------

A module made to use the SASC tree inference algorithm. 
An example workflow with this module could be:

from phylo import sasc as sc
from multiprocessing import cpu_count

# load the matrix and the mutation labels
gmat = sc.GenotypeMatrix.from_files(
    matrix_file = 'genotype.txt', 
    mutations_file = 'mutation_labels.txt', 
    matstring_format = 'SCITE', #SASC, SCITE and SPHYR formats are supported 
    get_cells = True
)

# infer a tree and get the expected genotype matrix for it, then store these to files.
out = sc.infer_tree(gmat, alphas = 0.3, beta = 0.1, k = 1, cores = cpu_count())
out['inferred_tree'].to_file('tree.gv')
out['inferred_tree'].with_visualization_features() \
    .draw_to_file('tree.png', show_support = True, show_color = False)
out['expected_genotype_matrix'].to_files('expected_genotype.txt')
"""

from phylo.sascpy import compute as infer_tree
from .core.genotypematrix import GenotypeMatrix
from .core.phylogenytree import SASCPhylogeny


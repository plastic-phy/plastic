"""

For more info about the components:
- help(mp3.PhylogenyTree)
- help(mp3.tree_similarity)

---------

Module that exposes the tree similarity measure presented at
https://github.com/AlgoLab/mp3treesim

Simple example workflow:

from phylo import mp3
import networkx as nx

# load a tree from a file in DOT format.
tree1 = mp3.PhylogenyTree.from_file('tree1.gv')

# build the other from a networkx graph.
tree2 = nx.DiGraph()
tree2.add_node('0', label = 'some,comma,separated,label,list')
...
tree2 = mp3.PhylogenyTree(tree2)

# compute similarity
similarity = tree_similarity(tree1, tree2, ...)

# store the second tree given that it wasn't stored and then
# visualize the results
tree2.to_file('tree2.gv')
tree1.draw_to_file('tree1.pdf'); tree2.draw_to_file('tree2.pdf')
with open('similarity.txt', w+) as f:
    f.write(similarity)
"""

import mp3treesim as mp3
from ._core.phylogenytree import PhylogenyTree


def tree_similarity(
        tree1,
        tree2,
        excluded1='',
        excluded2='',
        excluded_global='',
        ignore_unlabeled_nodes=False,
        mode='sigmoid',
        cores=1
):
    """
    Computes the similarity score between two trees. The score ranges from 0 to 1, with
    higher scores implying higher similarity.

    Parameters:
        tree1(PhylogenyTree),
        tree2(PhylogenyTree):
            The trees that need to be compared.
        excluded1(str) by default '',
        excluded2(str) by default '',
        excluded_global(str) by default '':
            Comma-separated sets of labels that will be ignored while reading the trees.
            excluded_1 and excluded_2 refer to the respective trees, excluded_global to both of them.
        ignore_unlabeled_nodes(bool), by default False:
            If this is true, then nodes without a 'label' attribute will be ignored while reading the trees.
        mode(str), by default 'sigmoid':
            The mode with which the trees will be compared. Must be one of 'sigmoid', 'union', 'intersection'
            or 'geometric'.
        cores(int), by default 1:
            The number of cores that will be used in the computation. Enter a value less than 1 to use all cores.

    Returns:
        float:
            The similarity score for the trees.
    """

    if int(cores) != cores:
        raise ValueError(f'the number of cores must be a positive integer, but {cores} is not.')

    if not isinstance(tree1, PhylogenyTree):
        raise TypeError(
            "tree1 needs to be a valid PhylogenyTree. Load it from a dot file with mp3.load()"
            "or make one from a networkx graph."
        )
    if not isinstance(tree2, PhylogenyTree):
        raise TypeError(
            "tree2 needs to be a valid PhylogenyTree. Load it from a dot file with mp3.load()"
            "or make one from a networkx graph."
        )

    excluded1 = list(set(excluded1.split(',') + excluded_global.split(',')))
    excluded2 = list(set(excluded2.split(',') + excluded_global.split(',')))

    tree1_as_mp3 = mp3.build_tree(tree1.as_digraph(), ignore_unlabeled_nodes, excluded1)
    tree2_as_mp3 = mp3.build_tree(tree2.as_digraph(), ignore_unlabeled_nodes, excluded2)

    # mp3 will handle bad values for the mode
    return mp3.similarity(tree1_as_mp3, tree2_as_mp3, mode=mode, cores=cores)

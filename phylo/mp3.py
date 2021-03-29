import mp3treesim as mp3
from .core.phylogenytree import PhylogenyTree


def tree_similarity(
    tree1,
    tree2,
    excluded1 = '',
    excluded2 = '',
    excluded_global = '',
    ignore_unlabeled_nodes = False,
    mode = 'sigmoid',
    cores = 1
    ):
    """
    Computes the similarity score between two trees.

    Parameters:
    - tree_1 (PhylogenyTree),
    - tree_2 (PhylogenyTree): The trees that need to be compared. 
    - excluded_1 (str) by default '',
    - excluded_2 (str) by default '',
    - excluded_global (str) by default '': Comma-separated sets of labels that will be
      ignored while reading the trees. excluded_1 and excluded_2 refer to the respective
      trees, excluded_global to both of them
    - ignore_unlabeled_nodes(bool), by default False: If this is true, then nodes without 
      a 'label' attribute will be ignored while reading the trees.
    - mode (str): The mode with which the trees will be compared. Must be one of:
      * 'sigmoid' (the default mode)
      * 'union'
      * 'intersection'
      * 'geometric'
    - cores (int), by default 1: The number of cores that will be used in the computation.
      set to 0 to use all the available cores.

    Returns:
    - float: The similarity score for the tree.
    """

    if not isinstance(tree1, PhylogenyTree):
        raise TypeError("tree1 needs to be a valid PhylogenyTree. Load it from a dot file with mp3.load().")
    if not isinstance(tree2, PhylogenyTree):
        raise TypeError("tree2 needs to be a valid PhylogenyTree. Load it from a dot file with mp3.load().")

    excluded1 = excluded1.split(',') + excluded_global.split(',')
    excluded2 = excluded2.split(',') + excluded_global.split(',')
    tree1_as_mp3 = mp3.build_tree(tree1.as_digraph(), ignore_unlabeled_nodes, excluded1)
    tree2_as_mp3 = mp3.build_tree(tree2.as_digraph(), ignore_unlabeled_nodes, excluded2)

    return mp3.similarity(tree1_as_mp3, tree2_as_mp3, mode = mode, cores = cores)

load_tree = PhylogenyTree.from_file

dump_tree = PhylogenyTree.to_file
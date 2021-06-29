
from ._core.genotypematrix import GenotypeMatrix
from ._core.phylogenytree import PhylogenyTree, SASCPhylogeny
from matplotlib import pyplot as plt

def _build_colormap(unclustered, clustered):
    mapping = dict()
    for ix, l in enumerate(clustered):
        for ll in l.split(','):
            mapping[ll] = ix
            
    colors = list()
    for l in unclustered:
        colors.append(mapping[l])
        
    return colors

def clusters(unclustered, clustered, ax=plt, **kwargs):
    if not isinstance(unclustered, GenotypeMatrix):
        raise TypeError(
            "unclustered needs to be a valid GenotypeMatrix."
        )
    if not isinstance(clustered, GenotypeMatrix):
        raise TypeError(
            "clustered needs to be a valid GenotypeMatrix."
        )

    from sklearn.decomposition import PCA

    pca = PCA(n_components=2)
    red = pca.fit_transform(unclustered.matrix().transpose())

    ax.scatter(x=red[:,0], y=red[:,1], c=_build_colormap(unclustered.mutation_labels, clustered.mutation_labels), **kwargs)


def _get_label_to_id_map(tree):
    return {
        node: '' if 'label' not in tree.nodes[node] else tree.nodes[node]['label']
        for node in tree
    }


def phylogeny(tree, show_labels=False, show_support=False, graphviz_positioning=True, **kwargs):
    import networkx as nx
    from colour import Color

    if not isinstance(tree, PhylogenyTree):
        raise TypeError(
            "tree needs to be a valid PhylogenyTree."
        )

    t = tree.as_digraph()
    
    c_gradient = list(Color("#3270FC").range_to(Color("#397D02"), 101))

    if show_labels:
        kwargs['labels'] = _get_label_to_id_map(t)
    
    if show_support:
        kwargs['node_color'] = [c_gradient[int(v)].hex for _, v in nx.get_node_attributes(t, 'support').items()]

    if graphviz_positioning:
        kwargs['pos'] = nx.nx_agraph.graphviz_layout(t, prog="dot")

    nx.draw(
        t,
        **kwargs
    )

def similarity_matrix(matrix, **kwargs):
    import seaborn as sns

    sns.heatmap(matrix, **kwargs)
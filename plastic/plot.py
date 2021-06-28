
from ._core.genotypematrix import GenotypeMatrix
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
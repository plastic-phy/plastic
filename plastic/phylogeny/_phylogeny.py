from ._sasc import sasc as _sasc
from .._core.genotypematrix import GenotypeMatrix 

def infer_phylogeny(genotype, **kwargs):

    if not isinstance(genotype, GenotypeMatrix):
        raise TypeError('genotype must be a valid GenotypeMatrix.')

    return _sasc(genotype, **kwargs)

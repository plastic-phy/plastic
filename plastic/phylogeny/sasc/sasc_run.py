from .._sasc import sasc as _sasc
from plastic import GenotypeMatrix 

def inference(genotype_matrix, **kwargs):

    if not isinstance(genotype_matrix, GenotypeMatrix):
        raise TypeError('genotype must be a valid GenotypeMatrix.')

    return _sasc(genotype_matrix, **kwargs)

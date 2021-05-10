from sascpy import compute
from .core.genotypematrix import GenotypeMatrix
from .core.phylogenytree import PhylogenyTree

infer_tree = compute

load_matrix = GenotypeMatrix.from_files
dump_matrix = GenotypeMatrix.to_files

load_tree = PhylogenyTree.from_file
dump_tree = PhylogenyTree.to_file

draw_tree = PhylogenyTree.draw_to_file

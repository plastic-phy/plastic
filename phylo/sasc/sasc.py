from .sascpy import compute
from .core.labeledmutationmatrix import LabeledMutationMatrix
from .core.phylogenytree import PhylogenyTree

infer_tree = sascpy.compute

load_matrix = LabeledMutationMatrix.from_files
dump_matrix = LabeledMutationMatrix.to_files

load_tree = PhylogenyTree.from_file
dump_tree = PhylogenyTree.to_file

draw_tree = PhylogenyTree.draw_to_file
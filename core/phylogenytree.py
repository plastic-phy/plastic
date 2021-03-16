import networkx as nx
from copy import deepcopy

class PhylogenyTree():

    class NotATreeError(Exception): pass

    def __init__(self, tree_as_nx_graph):

        if not isinstance(tree_as_nx_graph, nx.DiGraph):
            raise TypeError('the input must be a networkx graph')

        if any([not isinstance(label, str) for label in nx.get_node_attributes(drawtree, 'label').values()]):
            raise TypeError('if a node has a "label" attribute, then it must be a string.')

        if not nx.is_tree(tree_as_nx_graph):
            raise NotATreeError('the graph must be a tree.')

        self._tree = deepcopy(tree_as_nx_graph)

    def as_digraph(self, with_terminal_specimen = True):

        out = deepcopy(self._tree)
        if not with_terminal_specimen:
            leaves = [node for node in out if out.out_degree(node) == 0]
            out.remove_nodes_from(leaves)

        return out

    def draw(with_terminal_specimen = True):

        # Should this depend from mp3 instead? I think it would be a fine idea to let it be 
        # available from SASC as well.
        drawtree = nx.convert_node_labels_to_integers(self.as_digraph(with_terminal_specimen))

        labels = nx.get_node_attributes(drawtree, 'label')

        #Why does the original source do this?
        for node in drawtree.nodes(data = True):
            del node[1]['label']
        
        try:
            pos = nx.nx_pydot.pydot_layout(drawtree, prog = 'dot')
        except:
            pos = None

        nx.draw_networkx(drawtree, pos=pos, labels=labels)

    def from_dotstring(dot_string):
        
        nx_representation = nx.DiGraph(nx.drawing.nx_agraph.from_agraph(AGraph(string=dot_string)))
        return PhilogenyTree(nx_repr)

    def to_dotstring(self):

        pydot_repr = networkx.drawing.nx_pydot.to_pydot(self._tree)
        return pydot_repr.to_string()
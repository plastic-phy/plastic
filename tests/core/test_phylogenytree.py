from phylo.core.phylogenytree import *
import networkx as nx
import pytest as pt

def dummy_graph(with_unlabeled_node = True):
    G = nx.DiGraph()
    G.add_node('0', label = 'sabbia')
    G.add_node('1', label = 'pollo,boh')
    G.add_node('2')
    if not with_unlabeled_node:
        G.nodes['2']['label'] = 'sBin'
    G.add_edge('0', '1')
    G.add_edge('0', '2')
    return G

def dummy_tree(with_unlabeled_node = True, fully_labeled = False):
    
    return PhylogenyTree(dummy_graph(with_unlabeled_node = with_unlabeled_node), fully_labeled)

class TestInit:

    def test_standard_init(self):
        tree = dummy_tree()

    def test_init_with_null_node_but_fully_labeled_tree(self):
        with pt.raises(NotFullyLabeled):
            tree = dummy_tree(with_unlabeled_node = True, fully_labeled = True)

    def test_init_with_bad_label_type(self):
        G = nx.DiGraph()
        G.add_node('0', label = 1337)

        with pt.raises(TypeError):
            tree = PhylogenyTree(G)

    def test_init_with_empty_label(self):
        G = nx.DiGraph()
        G.add_node('0', label = '')
        
        with pt.raises(ValueError):
            tree = PhylogenyTree(G)

    def test_init_with_non_tree(self):
        not_a_tree = dummy_graph()
        not_a_tree.add_edge('2', '0')
        
        with pt.raises(NotATreeError):
            tree = PhylogenyTree(not_a_tree)

    def test_init_with_bad_node(self):
        integer_node = nx.DiGraph()
        integer_node.add_node(0)
        with pt.raises(TypeError):
            tree = PhylogenyTree(integer_node)

    def test_init_with_bad_graph_attribute(self):
        G = dummy_graph()
        G.graph['edge'] = 'this is not allowed!'
        with pt.raises(ValueError):
            tree = PhylogenyTree(G)

class TestAccessors:

    def test_graph_accessor(self):
        G = dummy_graph()
        T = dummy_tree()
        assert list(T.as_digraph().nodes(data = True)) == list(G.nodes(data = True))
        assert list(T.as_digraph().edges(data = True)) == list(G.edges(data = True))

class TestImmutability:

    def test_immutability_after_init(self):
        G = dummy_graph()
        T = PhylogenyTree(G)
        G.add_edge('2', '0')
        assert ('2', '0') not in T.as_digraph().edges

    def test_immutability_after_graph_access(self):
        T = dummy_tree()
        G = T.as_digraph()
        G.add_edge('2', '0')
        assert ('2', '0') not in T.as_digraph().edges

class TestSerialization:

    def test_serialization_roundtrip(self):
        T = dummy_tree()
        T_as_string = T.to_dotstring()
        T_after_roundtrip = PhylogenyTree.from_dotstring(T_as_string)

        assert list(T.as_digraph().graph) == list(T_after_roundtrip.as_digraph().graph)
        assert list(T.as_digraph().nodes(data = True)) == list(T_after_roundtrip.as_digraph().nodes(data = True))
        assert list(T.as_digraph().edges(data = True)) == list(T_after_roundtrip.as_digraph().edges(data = True))

class TestDrawing():

    def test_rendering_to_png_and_pdf(self):
        T = dummy_tree()
        T.draw_to_file('sabbia.png')
        T.draw_to_file('pollo.pdf')

    
        
    
        

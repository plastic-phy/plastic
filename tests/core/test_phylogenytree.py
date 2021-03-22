from phylo.core.phylogenytree import *
import networkx as nx
import pytest as pt

def dummy_graph(with_null_label = True):
    G = nx.DiGraph()
    G.add_node(0, label = 'sabbia', is_terminal = False)
    G.add_node(1, label = 'pollo', is_terminal = False)
    G.add_node(2, label = 'tasso' if not with_null_label else None, is_terminal = True)
    G.add_edge(0, 1)
    G.add_edge(0, 2)
    return G

def dummy_tree(with_null_label = True, fully_labeled = False):
    
    return PhylogenyTree(dummy_graph(with_null_label = with_null_label), fully_labeled)

class TestInit:

    def test_standard_init():
        tree = dummy_tree()

    def test_init_with_null_node_but_fully_labeled_tree():
        with pt.raises(ValueError):
            tree = dummy_tree(with_null_label = True, fully_labeled = True)

    def test_init_with_bad_label():
        G = nx.DiGraph()
        G.add_node(0, label = 1337, is_terminal = True)

        with pt.raises(TypeError):
            tree = PhylogenyTree(G)

    def test_init_with_empty_label():
        G = nx.DiGraph()
        G.add_node(0, label = '', is_terminal = True)
        
        with pt.raises(ValueError):
            tree = PhylogenyTree(G)

    def test_init_with_missing_attribute():
        G = nx.DiGraph()
        G.add_node(0, is_terminal = False)

        with pt.raises(ValueError):
            tree = PhylogenyTree(G)

    def test_init_with_non_tree():
        not_a_tree = dummy_graph()
        not_a_tree.add_edge(2, 0)
        
        with pt.raises(NotATreeError):
            tree = PhylogenyTree(not_a_tree)

class TestAccessors:

    def test_graph_accessor():
        G = dummy_graph()
        T = dummy_tree()
        assert T.as_digraph() == G

    def test_graph_without_leaves_accessor():
        T = dummy_tree()
        assert T.as_digraph(with_terminals = False).nodes == [0, 1]

    def test_terminals_accessor():
        T = dummy_tree()
        assert T.terminal_nodes == [(2, {label : None, is_terminal : True})]

class TestImmutability:

    def test_immutability_after_init():
        G = dummy_graph()
        T = PhylogenyTree(G, fully_labeled = False)
        G.add_edge(2, 0)
        assert (2, 0) not in T.as_digraph().edges

    def test_immutability_after_graph_access():
        T = dummy_tree()
        G = T.as_digraph()
        G.add_edge(2, 0)
        assert (2, 0) not in T.as_digraph().edges

    def test_immutability_after_terminals_access():
        T = dummy_tree()
        terminals = T.terminal_nodes
        terminals[0][1]['label'] = 'pollozorg'
        assert T.terminal_nodes[0][1]['label'] != 'pollozorg'

class TestSerialization:

    def test_serialization_roundtrip():
        T = dummy_tree()
        T_as_string = T.to_dotstring()
        T_after_roundtrip = PhylogenyTree.from_dotstring()

        assert T.as_digraph() == T_after_roundtrip.as_digraph()

class TestDrawing():

    def test_rendering_to_png_and_pdf():
        T = dummy_tree()
        T.draw_to_file('sabbia.png')
        T.draw_to_file('pollo.pdf', with_terminals = False)

    
        
    
        

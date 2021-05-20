import networkx as nx
from phylo._core.phylogenytree import SASCPhylogeny
import pytest as pt

class TestInit:

    def test_init_with_cells(self):
        a = nx.DiGraph()
        a.add_node('0', label = 'pollo,zorg')
        a.add_node('1', shape = 'box')
        a.add_edge('0', '1')
        a = SASCPhylogeny(a)

        assert a.__dict__['_has_cells'] == True

    def test_init_without_cells(self):
        a = nx.DiGraph()
        a.add_node('0', label = 'pollo,zorg')
        a.add_node('1', label = 'sabbia')
        a.add_edge('0', '1')
        a = SASCPhylogeny(a)

        assert a.__dict__['_has_cells'] == False

    def test_init_with_unlabeled_inner_node(self):
        a = nx.DiGraph()
        a.add_node('this_is_an_unlabeled_inner_node')
        a.add_node('1', label = 'sabbia')
        a.add_edge('this_is_an_unlabeled_inner_node', '1')
        
        with pt.raises(ValueError):
            a = SASCPhylogeny(a)

    def test_init_with_unlabeled_noncell(self):
        a = nx.DiGraph()
        a.add_node('0', label = 'pollo,zorg')
        a.add_node('unlabeled_noncell')
        a.add_edge('0', 'unlabeled_noncell')

        with pt.raises(ValueError):
            a = SASCPhylogeny(a)


        

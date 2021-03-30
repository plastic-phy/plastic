import networkx as nx
from copy import deepcopy
import pygraphviz
from io import StringIO


class NotATreeError(Exception): pass
class NotFullyLabeled(Exception): pass

class PhylogenyTree():

    def __init__(self, tree_as_nx_graph, fully_labeled = False):
        """
        Verifies that a networkx graph is a tree and validates its attributes, then
        it creates a PhylogenyTree instance for that graph.

        Parameters:
        - tree_as_nx_graph(networkx.DiGraph): a graph that represents a tree (a completely connected,
          acyclic and directed graph), and where for each node that has a label attribute the label is
          a non-empty string. 
        - fully_labeled(boolean), by default False: if this is true, initialization will fail if 
          the tree has unlabeled nodes.

        Returns:
          PhylogenyTree: the object that has been initialized with tree_as_nx_graph if the tree is valid.
        """
        if not isinstance(tree_as_nx_graph, nx.DiGraph):
            raise TypeError('the input must be a networkx graph.')

        if not nx.is_tree(tree_as_nx_graph):
            raise NotATreeError('the graph must be a tree.')

        # to allow correct serialization in DOT format, all attributes and node IDS must be strings,
        # and some graph attribute names must be reserved.
        # All of coupled with the ugly hack done during deserialization is making me think we should
        # think about using a different format.

        if any([attribute in tree_as_nx_graph.graph for attribute in {'edge', 'node', 'graph'}]):
            bad_keys = set(tree_as_nx_graph.graph.keys()).intersection({'edge', 'node', 'graph'})
            raise ValueError(f'graph attributes with keys "edge", "node" or "graph", are not allowed, but {bad_keys} were present')
        for (key, value) in tree_as_nx_graph.graph.items():
            if not isinstance(value, str):
                raise TypeError(f'graph attribute {key} has non-string value {value}')
            
        for (u, v, attributes) in tree_as_nx_graph.edges(data = True):
            for (key, value) in attributes.items():
                if not isinstance(value, str):
                    raise TypeError(f'{key} attribute of edge {u} -> {v} has non-string value {value}')
                
        for (node, attributes) in tree_as_nx_graph.nodes(data = True):
            if not isinstance(node, str):
                raise TypeError(f'all nodes must be identified by strings, but {node} is not')
            for (key, value) in attributes.items():
                if not isinstance(value, str):
                    raise TypeError(f'{key} attribute of node {node} has non-string value {value}')
            if 'label' in attributes:
                label = attributes['label']
                if len(label) == 0:
                    raise ValueError(f'the node {node} has an empty label list.')
                if any([len(single_label) == 0 for single_label in label.split(',')]):
                    raise ValueError(f'the node {node} with label list {label} has empty labels in its label list.')
            elif fully_labeled:
                raise NotFullyLabeled()
            

        # In order to ensure the object's immutability, everything is copied over.
        self._tree = deepcopy(tree_as_nx_graph)

    # the external representation gives an independent copy. 
    def as_digraph(self):
        """
        Returns a networkx representation of the tree as a networkx digraph.
        The returned tree is independent from the internal representation.
        """
        out = deepcopy(self._tree)
        
        return out

    def draw_to_file(self, file_path):
        """
        Draws the tree to a file using a dot layout. Requires a Graphviz installation.

        Parameters: 
        - file_path(string): The file in which the tree will be drawn. The tested use cases are drawing the tree
          as an image (file with .png extension) or as a PDF (file with .pdf extension).

        Returns: nothing.

        Side effects: 
        - Draws the tree to the file using a dot layout. Reserved dot attributes will work as specified.
          If a node is not labeled, then its ID will be used as a label (alongside with a warning). If the file 
          specified by file_path doesn't exist, it will be created; if it exists, **its content will be overwritten**.
        """

        # Should this depend from mp3 instead? I think it would be a fine idea to let it be 
        # available from SASC as well.
        drawtree = self.as_digraph()

        # The nodes will be labeled with their numerical ID if a label isn't present.
        for [node, attributes] in drawtree.nodes(data = True):
            if 'label' not in attributes:
                attributes['label'] = f'no label for node with ID: {node}'

        drawtree = nx.nx_agraph.to_agraph(drawtree)
        try:
            drawtree.layout(prog = 'dot')
        except:
            raise Exception('This feature requires a GraphViz installation.')
        drawtree.draw(file_path)

    def from_dotstring(dot_string):
        """
        Validates the content of a string representing a graph in dot format, then
        initializes and returns a PhylogenyTree representation if the string is valid.
        """
        
        nx_representation = nx.nx_agraph.from_agraph(pygraphviz.AGraph(string = dot_string))

        # Converting to an AGraph adds its own set of default attributes that need
        # to be removed before conversion to a NetworkX graph.
        # Still, this is quite the weird hack and I'd much rather not do it.
        del nx_representation.graph['graph']
        del nx_representation.graph['node']
        del nx_representation.graph['edge']
        return PhylogenyTree(nx_representation)

    def to_dotstring(self):
        """
        Dumps the tree to a string with its dot representation.
        """
        gv_representation = nx.drawing.nx_agraph.to_agraph(self._tree)
        return gv_representation.to_string()

    def from_file(dotfile_path):
        """
        Validates and loads a PhylogenyTree from the specified file path, if the file exists.
        """
        # Please introduce me to the art of actually checking files. Thank you!
        with open(dotfile_path, 'r') as f:
            return PhylogenyTree.from_dotstring(f.read())

    def to_file(self, dotfile_path):
        """
        Dumps a PhylogenyTree to the specified file. If the file doesn't exist, it will be created; if
        it exists, **it will be overwritten**.
        """
        with open(dotfile_path, 'w+') as f:
            tree_as_dot = self.to_dotstring()
            f.write(tree_as_dot)

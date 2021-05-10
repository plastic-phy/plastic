import networkx as nx
from copy import deepcopy
import pygraphviz


class NotATreeError(Exception): pass
class NotFullyLabeled(Exception): pass


class PhylogenyTree:

    def __init__(self, tree_as_nx_graph, fully_labeled = False):
        """
        Verifies that a networkx graph is a tree and validates its attributes, then
        it creates a PhylogenyTree instance for that graph.

        Parameters:
            tree_as_nx_graph(networkx.DiGraph):
                A networkx graph that represents a directed tree.
                Each node attribute must be a string, and graph attribute keys "node", "edge" and
                "graph" are reserved.
                If a node has the "label" attribute, then it must be a comma-separated list of
                non-empty strings.
            fully_labeled(bool), by default False:
                If this is true, initialization will fail if the tree has unlabeled nodes.

        Returns:
            PhylogenyTree:
                The object that has been initialized with tree_as_nx_graph if the tree is valid.
                Its internal state is fully independent from the tree that was used for the
                initialization, so subsequent alteration to it won't alter this object.
        """
        if not isinstance(tree_as_nx_graph, nx.DiGraph):
            raise TypeError('the input must be a networkx graph.')

        if not nx.is_arborescence(tree_as_nx_graph):
            raise NotATreeError('the graph must be a tree.')

        # To allow correct serialization in DOT format, all attributes and node IDs must be strings.
        # The graph attributes "node", "edge" and "graph" are reserved due to oddities with networkx
        # to AGraph conversion.
        if any([attribute in tree_as_nx_graph.graph for attribute in {'edge', 'node', 'graph'}]):
            bad_keys = set(tree_as_nx_graph.graph.keys()).intersection({'edge', 'node', 'graph'})
            raise ValueError(f'graph attributes with keys "edge", "node" or "graph" are not allowed, but {bad_keys} are present')
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
        The returned tree is independent from the internal representation, so subsequent
        changes to it don't alter the state of this object.
        """
        out = deepcopy(self._tree)
        
        return out

    def draw_to_file(self, file_path):
        """
        Draws the tree to a file using a dot layout. Requires a Graphviz installation.

        Parameters:
            file_path(string):
                The file in which the tree will be drawn. The tested use cases are drawing the tree
                as an image or as a PDF (file with .pdf extension).

        Returns: nothing.

        Side effects:
            Draws the tree to the file using a dot layout. Reserved dot attributes will work as specified.
            If a node is not labeled, then its ID will be used as a label (alongside with a warning).
            If the file specified by file_path doesn't exist, it will be created; if it exists,
            its content will be overwritten.
        """
        drawtree = self.as_digraph()

        # The nodes will be labeled with their numerical ID if a label isn't present.
        for [node, attributes] in drawtree.nodes(data = True):
            if 'label' not in attributes:
                attributes['label'] = f'no label for node with ID: {node}'

        drawtree = nx.nx_agraph.to_agraph(drawtree)
        drawtree.layout(prog = 'dot')
        drawtree.draw(file_path)

    def from_dotstring(dot_string):
        """
        Validates the content of a string representing a graph in dot format, then
        initializes and returns a PhylogenyTree representation if the string is valid.
        """
        
        nx_representation = nx.nx_agraph.from_agraph(pygraphviz.AGraph(string = dot_string))

        # Converting to an AGraph adds its own set of default attributes that need
        # to be removed before conversion to a networkx graph.
        del nx_representation.graph['graph']
        del nx_representation.graph['node']
        del nx_representation.graph['edge']
        return PhylogenyTree(nx_representation)

    def to_dotstring(self):
        """
        Dumps the tree to a DOT representation string.
        """
        gv_representation = nx.drawing.nx_agraph.to_agraph(self._tree)
        return gv_representation.to_string()

    def from_file(file_path):
        """
        Validates and loads a PhylogenyTree from the specified file path, if the file exists.
        """
        # Please introduce me to the art of actually checking files. Thank you!
        with open(file_path, 'r') as f:
            return PhylogenyTree.from_dotstring(f.read())

    def to_file(self, file_path):
        """
        Dumps a PhylogenyTree to the specified file in DOT format.
        If the file doesn't exist, it will be created; if it exists, **it will be overwritten**.
        """
        with open(file_path, 'w+') as f:
            tree_as_dot = self.to_dotstring()
            f.write(tree_as_dot)

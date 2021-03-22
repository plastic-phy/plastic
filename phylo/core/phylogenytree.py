import networkx as nx
from copy import deepcopy


class NotATreeError(Exception): pass
class NotFullyLabeled(Exception): pass

class PhylogenyTree():

    def __init__(self, tree_as_nx_graph, fully_labeled = True):
        """
        Verifies that a networkx graph is a tree and validates its attributes, then
        it creates a PhylogenyTree instance for that graph.

        Parameters:
        - tree_as_nx_graph(networkx.DiGraph): a graph that represents a tree (a completely connected,
          acyclic and directed graph), with each node possessing the following attributes:
            * is_terminal(boolean): true if the node represents an individual in the sample that was
              used to infer the phylogeny tree, false if the node represents a set of mutations or a 
              set of deletions of mutations in the inferred phylogeny.
            * label(string), by default True: the name of the set of mutations or deletions represented 
              by the node if the node is not terminal, the name of the individual of the sample if the node
              is terminal. Multiple mutations must be separated by commas, and each label must be a
              non-empty string.
              If fully_labeled is false, then the value None is allowed in place of a string to 
              signify the lack of labels in the node.
          The internal representation of the tree will be independent from the graph that was used to
          initialize it, so subsequent changes to tree_as_nx_graph will not reflect on the PhylogenyTree.
        - fully_labeled(boolean): if this is true, initialization will fail if the tree has unlabeled nodes.

        Returns:
          PhylogenyTree: the object that has been initialized with tree_as_nx_graph if the tree is valid.
        """
        if not isinstance(tree_as_nx_graph, nx.DiGraph):
            raise TypeError('the input must be a networkx graph.')

        if not nx.is_tree(tree_as_nx_graph):
            raise NotATreeError('the graph must be a tree.')

        # Validating the presence of the attributes sets a convention on the tree's representation.
        for (node, node_attributes) in tree_as_nx_graph.nodes(data = True):
            try:
                label = node_attributes['label']
                is_terminal = node_attributes['is_terminal']
            except KeyError as e:
                raise ValueError(f'node with id {node} does not have the necessary attribute {str(e)}')
            if label is None and fully_labeled:
                raise NotFullyLabeled(f'node with id {node} does not have a label, but the tree was specified as fully labeled.')
            if label is not None:
                if not isinstance(label, str):
                    raise TypeError(f'label {label} of node with id {node} must be a string (or None if the node is not labeled).')
                if len(label) == 0:
                    raise ValueError(f'label of node with id {node} must be non-empty.')
            if is_terminal not in {True, False}:
                raise TypeError(f'is_terminal must be a boolean, but node with id {node} has value {is_terminal} for it.')

        # In order to ensure the object's immutability, everything is copied over.
        self._tree = deepcopy(tree_as_nx_graph)
        self.terminal_nodes = [(node, attributes) for (node, attributes) in self._tree.nodes(data = True) if attributes['is_terminal']]

    # the external representation gives an independent copy. 
    def as_digraph(self, with_terminals = True):
        """
        Returns a networkx representation of the tree as a digraph. Gives back the tree without terminal
        nodes if so specified.

        Parameters:
        - with_terminals(boolean), by default True: if this is true, the terminal nodes will be included
          in the representation.

        Returns: 
        - nx.DiGraph: a mutable copy of the tree that the object has been initialized with, 
          with or without the terminal nodes.
        """
        out = deepcopy(self._tree)

        if not with_terminals:
            out.remove_nodes_from(self._terminal_nodes)

        return out

    @property
    def terminal_nodes(self):
        """
        Returns a list of the terminal nodes of the tree.

        Parameters: none.

        Returns: 
        - list: a list in which each element is a tuple of two components: the id of a terminal node, 
          and its attributes as a dictionary. 
        """
        return deecopy(self._terminal_nodes)

    def draw_to_file(file_path, with_terminals = True):
        """
        Draws the tree to a file using a dot layout. Requires a Graphviz installation.

        Parameters: 
        - file_path(string): The file in which the tree will be drawn. The tested use cases are drawing the tree
          as an image (file with .png extension) or as a PDF (file with .pdf extension).
          with_terminals(boolean), by default False: if true, the terminal nodes will be drawn as well.

        Returns: nothing.

        Side effects: 
        - Draws the tree to the file using a dot layout. Reserved dot attributes will work as specified.
          If a node is not labeled, then its ID will be used as a label (alongside with a warning). If the file 
          specified by file_path doesn't exist, it will be created; if it exists, **its content will be overwritten**.
        """

        # Should this depend from mp3 instead? I think it would be a fine idea to let it be 
        # available from SASC as well.
        drawtree = self.as_digraph(with_terminals = with_terminals)

        # The nodes will be labeled with their numerical ID if a label isn't present.
        for (node, attributes) in drawtree:
            if attributes['label'] is None:
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
        nx_representation = nx.nx_agraph.from_agraph(pygraphviz.AGraph(graph))
        return PhilogenyTree(nx_representation)

    def to_dotstring(self):
        """
        Dumps the tree to a string with its dot representation.
        """
        pydot_representation = networkx.drawing.nx_pydot.to_pydot(self._tree)
        return pydot_representation.to_string()
import networkx as nx
from copy import deepcopy
from collections import defaultdict, deque
import pygraphviz
from colour import Color
import phylo.core.sascviz as sv


class NotATreeError(Exception): pass
class NotFullyLabeled(Exception): pass


class PhylogenyTree:

    def __init__(self, tree_as_nx_graph, fully_labeled = False, separator = ','):
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
            separator(str), by default ',':
                the separator for the lists of labels in the tree.

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
        for key in tree_as_nx_graph.graph:
            tree_as_nx_graph.graph[key] = str(tree_as_nx_graph.graph[key])
        if any([attribute in tree_as_nx_graph.graph for attribute in {'edge', 'node', 'graph'}]):
            bad_keys = set(tree_as_nx_graph.graph.keys()).intersection({'edge', 'node', 'graph'})
            raise ValueError(f'graph attributes with keys "edge", "node" or "graph" are not allowed, but {bad_keys} are present')
            
        for [u, v, attributes] in tree_as_nx_graph.edges(data = True):
            for key in attributes:
                attributes[key] = str(attributes[key])
                
        for (node, attributes) in tree_as_nx_graph.nodes(data = True):
            if not isinstance(node, str):
                raise TypeError(f'all nodes must be identified by strings, but {node} is not')
            for key in attributes:
                attributes[key] = str(attributes[key])
            if 'label' in attributes:
                label = attributes['label']
                if len(label) == 0:
                    raise ValueError(f'the node {node} has an empty label list.')
                if any([len(single_label) == 0 for single_label in label.split(separator)]):
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

    @classmethod
    def from_dotstring(cls, dot_string):
        """
        Validates the content of a string representing a graph in dot format, then
        initializes and returns a tree representation if the string is valid.
        """
        
        nx_representation = nx.nx_agraph.from_agraph(pygraphviz.AGraph(string = dot_string))

        # Converting to an AGraph adds its own set of default attributes that need
        # to be removed before conversion to a networkx graph.
        del nx_representation.graph['graph']
        del nx_representation.graph['node']
        del nx_representation.graph['edge']
        return cls(nx_representation)

    def to_dotstring(self):
        """
        Dumps the tree to a DOT representation string.
        """
        gv_representation = nx.drawing.nx_agraph.to_agraph(self._tree)
        return gv_representation.to_string()

    @classmethod
    def from_file(cls, file_path):
        """
        Validates and loads a tree from the specified file path, if the file exists.
        """
        # Please introduce me to the art of actually checking files. Thank you!
        with open(file_path, 'r') as f:
            return cls.from_dotstring(f.read())

    def to_file(self, file_path):
        """
        Dumps the tree to the specified file using DOT format.
        If the file doesn't exist, it will be created; if it exists, **it will be overwritten**.
        """
        with open(file_path, 'w+') as f:
            tree_as_dot = self.to_dotstring()
            f.write(tree_as_dot)


class MergingRootError(Exception): pass


class SASCPhylogeny(PhylogenyTree):

    def __init__(self, tree_as_nx_graph):
        """
        Verifies that a networkx graph uses SASC conventions to represent deletions, mutations and cells
        and in which each inner node of the tree is labeled, and then returns an extended PhylogenyTree with support
        for operations that can only be done with SASC trees.

        Parameters:
            tree_as_nx_graph(networkx.DiGraph):
                A networkx graph that represents a directed tree.
                Each node attribute must be a string, and graph attribute keys "node", "edge" and
                "graph" are reserved.
                If a node has the "label" attribute, then it must be a comma-separated list of
                non-empty strings.
                A node in the tree can be a cell node (unlabeled and with a box shape), or a mutation
                node (labeled). Mutation nodes can be deletion nodes (with fillcolor = indianred1).
                All the inner nodes in the tree must be mutation nodes.
            with_cells(bool), by default False:
                If this is False, then cells cannot be present in the tree.

        Returns:
            SASCPhylogeny:
                The representation of the tree as a SASCPhylogeny; the returned object is also
                a PhylogenyTree instance.
        """

        super().__init__(tree_as_nx_graph, fully_labeled = False, separator = ',')
        self._has_cells = False

        for node in self._tree:
            if 'label' not in self._tree.nodes[node]:
                if self._tree.out_degree[node] != 0:
                    raise ValueError('All the inner nodes of a SASCPhylogeny must be labeled.')
                elif self._tree.nodes[node].get('shape') != 'box':
                    raise ValueError('If a leaf node is unlabeled, then it must be a cell, marked by a box shape attribute.')
                else:
                    self._has_cells = True

    def with_visualization_features(self, support_threshold = 0, collapse_simple_paths = False):
        """
        Creates a modified tree that can be used to get a clearer visualization for the phylogeny.

        Parameters:
            support_threshold(int) by default 0:
                If this is set to anything but 0, then mutation nodes with a lower support
                than the threshold will be collapsed into their parents until only nodes with
                a large enough support remain.
            collapse_simple_paths(bool) by default False:
                If this is true, then simple chains of mutation nodes will be collapsed into
                a single node. If a support threshold was specified, then this will happen as
                a second step.

        Returns:
            SASCPhylogeny:
                A SASCPhylogeny in which the support has been computed for each mutation node, with
                the desired modifications applied to it. The leaves will be removed.
        """

        svtree = self.to_sv_tree()
        sv.collapse_low_support(svtree, svtree.root, support_threshold)
        sv.collapse_simple_paths(svtree, svtree.root)

        sv.calc_supports(svtree.root, defaultdict(int, {svtree.root.id: 1}))
        return SASCPhylogeny.from_sv_tree(svtree)

    def without_cells(self):
        tree = self.as_digraph()
        for node in list(tree):
            if 'label' not in tree.nodes[node]:
                tree.remove_node(node)
        return SASCPhylogeny(tree)

    def draw_to_file(self, file_path, show_support = True, show_color = True):

        """
        Works the same as PhylogenyTree.draw_to_file, but with the options to show/hide node colors
        that visualize how much support a mutation has, and to show/hide the support for each node
        in its label. Said information will only be shown if present in the tree.
        """
        drawtree = self.as_digraph()
        if show_support:
            for node in drawtree:
                if 'support' in drawtree.nodes[node] and 'label' in drawtree.nodes[node]:
                    drawtree.nodes[node]['label'] += '\n[s = {0}]'.format(drawtree.nodes[node]['support'])

        if not show_color:
            for node in drawtree:
                if 'color' in drawtree.nodes[node]:
                    del drawtree.nodes[node]['color']

        PhylogenyTree(drawtree).draw_to_file(file_path)

    def to_sv_tree(self):
        root = [node for node, degree in self._tree.in_degree() if degree == 0][0]

        curr_node = sv.Node(root)
        out = sv.Tree(curr_node)
        out.confidence_score = (
            self._tree.graph['Confidence score'] if 'Confidence score' in self._tree.graph
            else None
        )

        visit_queue = deque([curr_node])
        while visit_queue:
            curr_node = visit_queue.popleft()

            if self._tree.nodes[curr_node.id].get('fillcolor') == 'indianred1':
                out.set_deletion(curr_node.id)
            out.set_mutations(curr_node.id, self._tree.nodes[curr_node.id]['label'].split(','))

            for child in self._tree.successors(curr_node.id):
                if 'label' in self._tree.nodes[child]:
                    to_add = sv.Node(child)
                    visit_queue.append(to_add)

                    out.add_node(to_add)
                    out.add_edge(curr_node.id, to_add.id)

                else:
                    curr_node.support += 1
        
        sv.calc_supports(out.root, defaultdict(int, {root: 1}))
        return out

    @classmethod
    def from_sv_tree(cls, svtree):

        c_red = Color("#FF1919")
        c_green = Color("#397D02")
        c_blue = Color("#3270FC")
        c_gradient = list(c_blue.range_to(c_green, 100))
        
        out = nx.DiGraph()
        out.add_node(svtree.root.id)

        visit_queue = deque([svtree.root])
        while visit_queue:
            curr_node = visit_queue.popleft()

            data = dict()
            if curr_node.deletion:
                data['fillcolor'] = 'indianred1'
                data['style'] = 'filled'
            data['support'] = curr_node.get_s()
            data['label'] = ','.join(curr_node.mutations)
            data['color'] = (
                c_red if data['support'] == 0
                else c_gradient[data['support'] - 1]
            )
            out.nodes[curr_node.id].update(**data)

            for child in curr_node.children:
                visit_queue.append(child)

                out.add_edge(curr_node.id, child.id)

        if svtree.confidence_score is not None:
            out.graph['Confidence score'] = svtree.confidence_score
        out.graph['labelloc'] = 't'
        out.graph['penwidth'] = 2

        return cls(out)


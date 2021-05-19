#     MIT License
#
#     Copyright (c) 2017-2019 Simone Ciccolella
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy
#     of this software and associated documentation files (the "Software"), to deal
#     in the Software without restriction, including without limitation the rights
#     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#     copies of the Software, and to permit persons to whom the Software is
#     furnished to do so, subject to the following conditions:
#
#     The above copyright notice and this permission notice shall be included in all
#     copies or substantial portions of the Software.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#     SOFTWARE.
#

import sys


class Node:
    def __init__(self, id):
        self.id = id
        self.mutations = None
        self.deletion = False

        self.parent = None
        self.children = []

        self.depth = 0

        self.c_grad = None

        self.support = 0
        self.cumulative_support = 0
        self.downstream_support = 0

        self.tot_cells = 0
        self.show_support = False
        self.show_color = False
    
    def get_name(self, sep=','):
        if not self.deletion:
            return sep.join(self.mutations)
        else:
            return sep.join('%s-' % x for x in self.mutations)

    def get_s(self):
        # return int((self.downstream_support / (self.tot_cells - self.parent.cumulative_support))*100)
        # print(self.downstream_support, self.parent.downstream_support - self.parent.support)
        try:
            return int(
                (self.downstream_support / (self.parent.downstream_support - self.parent.support)) * 100
            )
        except:
            return 0

    def calc_cumulative_sup(self):
        if self.parent:
            self.cumulative_support = self.parent.cumulative_support + self.support
        else:
            self.cumulative_support = self.support


class Tree:
    def __init__(self, root):
        self.root = root
        self.id_to_node = {}
        self.mut_to_node = {}
        self.deletions = []
        self.edges = []
        self.tree_label = None

        self.add_node(root)
    
    def add_node(self, node):
        self.id_to_node[node.id] = node

    def remove_node(self, node):
        self.id_to_node.pop(node.id)
        for m in node.mutations:
            self.mut_to_node.pop(m)
        
        node.parent.children.remove(node)
        for child in node.children:
            child.parent = node.parent
            node.parent.children.append(child)
    
    def pop_node(self, node):
        self.id_to_node.pop(node.id)
        for m in node.mutations:
            self.mut_to_node.pop(m)
            
        node.parent.children.remove(node)
        for child in node.children:
            child.parent = node.parent

    def merge_nodes(self, merged, to_merge):
        if not to_merge.parent == merged:
            sys.exit('Merge two non-parent-child nodes')

        for m in to_merge.mutations:
            merged.mutations.append(m)
        self.remove_node(to_merge)
        for m in merged.mutations:
            self.mut_to_node[m] = merged
        
        merged.support += to_merge.support

    def get_node(self, id):
        try:
            return self.id_to_node[id]
        except:
            return None

    def add_edge(self, start_id, end_id):
        s_node = self.get_node(start_id)
        e_node = self.get_node(end_id)

        e_node.parent = s_node
        s_node.children.append(e_node)
        self.edges.append((start_id, end_id))
        e_node.cumulative_support = s_node.cumulative_support

    def set_mutations(self, id, mutations):
        self.get_node(id).mutations = mutations
        for mut in mutations:
            if not self.get_node(id).deletion:
                self.mut_to_node[mut] = self.get_node(id)

    def set_deletion(self, id):
        self.get_node(id).deletion = True
        self.deletions.append(self.get_node(id))

    def is_ancestor(self, anc_mut, node_mut):
        anc = self.mut_to_node[anc_mut]
        node = self.mut_to_node[node_mut]
        
        p = node.parent
        while p:
            if p == anc:
                return True
            
            p = p.parent
        return False

    def get_deletions_name(self):
        ret = []
        for d in self.deletions:
            ret.append(d.get_name())
        return ret


def collapse_simple_paths(tree, node):
    if len(node.children) == 0:
        return
    else:
        if not node.deletion:
            if len(node.children) == 1 and node.parent and not node.children[0].deletion:
                # collapse parent with its child
                only_child = node.children[0]

                tree.merge_nodes(node, only_child)
                collapse_simple_paths(tree, node)
        for child in node.children:
            collapse_simple_paths(tree, child)


def collapse_low_support(tree, node, support_th):
    if node.parent and not node.deletion and not node.parent.deletion:
        if node.get_s() < support_th and not node.parent == tree.root:
            # collapse node with its parent
            tree.merge_nodes(node.parent, node)
            
            node = node.parent
    for child in node.children:
        collapse_low_support(tree, child, support_th)


def delete_subtree(tree, node):
    if len(node.children) == 0:
        tree.pop_node(node)
    else:
        for child in node.children:
            delete_subtree(tree, child)
        tree.pop_node(node)


def calc_supports(node, level_count):
    if node.parent:
        node.depth = node.parent.depth + 1
    level_count[node.depth] += 1
    if len(node.children) == 0:
        node.calc_cumulative_sup()
        node.downstream_support = node.support
    else:
        node.calc_cumulative_sup()
        for child in node.children:
            calc_supports(child, level_count)
        node.downstream_support = sum([x.downstream_support for x in node.children]) + node.support

#!/usr/bin/env python
###############################################################################
#
# graph.py
#
# Copyright (C) 2011 Christopher Davoren
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################
""" Represents a directed graph """

__author__ = "Christopher Davoren"
__email__ = "davoren@wehi.edu.au"
__credits__ = ["Tony Papenfuss", "Arther Hsu"]

class Graph:
    def __init__(self):
        self.nodes = []
        self.node_children = {}
        self.node_parents = {}
    
    def add_node(self, node):
        if node in self.nodes:
            return

        self.nodes.append(node)
        self.node_children[node] = []
        self.node_parents[node] = []

    def remove_node(self, node):
        if not node in self.nodes:
            return

        self.nodes.remove(node)
        del self.node_children[node]
        del self.node_parents[node]

    def add_arc(self, parent_node, child_node, check_cyclic = False):
        if not parent_node in self.nodes:
            raise Exception("Unable to add arc: parent node " + str(parent_node) + " is not in the node list")

        if not child_node in self.nodes:
            raise Exception("Unable to add arc: child node " + str(child_node) + " is not in the node list") 

        # Check if adding this arc would make the graph cyclic
        if check_cyclic:
            visited_nodes = [parent_node,]
            result = self._check_cyclic(child_node, visited_nodes)
            if result:
                raise Exception("Adding an arc from " + str(parent_node) + " to " + str(child_node) + " would create a cyclic graph")

        if not child_node in self.node_children[parent_node]: 
            self.node_children[parent_node].append(child_node)

        if not parent_node in self.node_parents[child_node]:
            self.node_parents[child_node].append(parent_node)

    def remove_arc(self, parent_node, child_node):
        if not parent_node in nodes:
            return

        if not child_node in nodes:
            return

        self.node_children[parent_node].remove(child_node)
        self.node_parents[child_node].remove(parent_node)

    def get_head_nodes(self):
        head_nodes = []

        for node, parents in self.node_parents.items():
            if not len(parents):
                head_nodes.append(node)

        return head_nodes
        
    
    def get_tail_nodes(self):
        tail_nodes = []

        for node, children in self.node_children.items():
            if len(children) == 0:
                tail_nodes.append(node)

        return tail_nodes

    def __str__(self):
        summary = ""

        for parent, children in self.node_children.iteritems():
            for child in children:
                summary += str(parent) + " -> " + str(child) + "\n"

        return summary

    def print_summary(self):
        print "Node list:"
        for node in self.nodes:
            print "\t", node

        summary = str(self)

        print "Node connections:"
        for line in summary.splitlines():
            print "\t" + line

    def _check_cyclic(self, node, visited_nodes):
        for child in self.node_children[node]:
            if child in visited_nodes:
                return True

            visited_nodes.append(child)
            cyclic = self._check_cyclic(child, visited_nodes)

            if cyclic:
                return True

            visited_nodes.remove(child)

        return False

    def is_cyclic(self):
        if len(self.nodes) == 0:
            return False

        head_nodes = self.get_head_nodes()

        if len(head_nodes) == 0:
            return True

        for node in self.nodes:
            cyclic = self._check_cyclic(node, [])
            if cyclic:
                return True

        return False

    def add_paired_list(self, paired_list, check_cyclic = False):
        for pair in paired_list:
            self.add_node(pair[0])
            self.add_node(pair[1])
            self.add_arc(pair[0], pair[1], check_cyclic)

    def topological_sort(self):
        def visit_node(graph, node, node_list):
            if node in node_list:
                return
            for parent in graph.node_parents[node]:
                visit_node(graph, parent, node_list)
            node_list.append(node)

        node_list = []
        for node in self.get_tail_nodes():
            visit_node(self, node, node_list)

        return node_list

    def has_parents(self, node):
        if not self.node_parents[node]:
            raise NodeNotFoundError("The given node is not present in this graph")

        return len(self.node_parents[node]) > 0

    def has_children(self, node):
        if not self.node_children[node]:
            raise NodeNotFoundError("The given node is not present in this graph")

        return len(self.node_children[node]) > 0

    def get_parents(self, node):
        return self.node_parents[node]

class NodeNotFoundError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

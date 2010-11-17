#!/usr/bin/env python
# encoding: utf-8


"""
randomise_normal.py

Created by Moritz Beber on 2010-10-28.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""


import networkx
import random


def _add_edge(graph, src, tar, bunch):
    graph.add_edge(src, tar)
    bunch.append((src, tar))

def _remove_edge(graph, src, tar, bunch):
    graph.remove_edge(src, tar)
    bunch.remove((src, tar))

def _check_conditions(graph, first, second):
    # curiously the conditions for switching single and double edges are
    # are the same for just slightly different reasons
    if first == second:
        return False
    # prevent creation of self-loops
    if first[0] == second[1]:
        return False
    if second[0] == first[1]:
        return False
    # check if we would flip to existing edges
    if graph.has_edge(first[0], second[1]):
        return False
    if graph.has_edge(second[0], first[1]):
        return False
    # check if we would create a bidirectional edge
    # or cover existing reverse edge in double edge switching
    if graph.has_edge(second[1], first[0]):
        return False
    if graph.has_edge(first[1], second[0]):
        return False
    return True

def _switch_double(graph, first, second, group):
    if _check_conditions(graph, first, second):
        # if all of these conditions are met, switch double edge
        # add the forward direction
        _add_edge(graph, first[0], second[1], group)
        _add_edge(graph, second[0], first[1], group)
        # add the reverse direction
        _add_edge(graph, first[1], second[0], group)
        _add_edge(graph, second[1], first[0], group)
        # remove old edges
        _remove_edge(graph, first[0], first[1], group)
        _remove_edge(graph, first[1], first[0], group)
        _remove_edge(graph, second[0], second[1], group)
        _remove_edge(graph, second[1], second[0], group)
        return True
    else:
        return False

def _switch_single(graph, first, second, group):
    if _check_conditions(graph, first, second):
        # all conditions passed
        _add_edge(graph, first[0], second[1], group)
        _add_edge(graph, second[0], first[1], group)
        _remove_edge(graph, first[0], first[1], group)
        _remove_edge(graph, second[0], second[1], group)
        return True
    else:
        return False

def _make_groups(graph):
    uni = list() # group of unidirectional edges
    bi = list() # group of bidirectional edges
    for (src, tar) in graph.edges_iter():
        if graph.has_edge(tar, src):
            bi.extend([(src, tar), (tar, src)])
        else:
            uni.append((src, tar))
    return (uni, bi)

def _flip_edges(graph, num, bunch, method):
    success = 0
    for i in xrange(num):
        first = random.choice(bunch)
        second = random.choice(bunch)
        if method(graph, first, second, bunch):
            success += 1
    return success

def randomise(graph, flip=100):
    """
'randomise' produces a randomised version of the supplied directed
graph 'graph'. This is achieved by switching edges in the graph a number of
times equal to the number of edges present times 'flip'. This function is
intended as a basis for calculating three-node subgraphs according to the
standard mfinder method. As such only degrees of nodes, bi-directional link properties, etc.
(see:
R Milo, S Shen-Orr, S Itzkovitz, N Kashtan, D Chklovskii & U Alon,
Network Motifs: Simple Building Blocks of Complex Networks
Science, 298:824-827 (2002)
) are conserved. For larger subgraph statistics also the smaller subgraph
statistics would have to be conserved.
    """
    rnd_graph = graph.copy()
    success = 0
    sets = _make_groups(rnd_graph) # treat unidirectional and bidirectional
    # edges separately
    w_uni = len(sets[0])
    w_bi = len(sets[1])
    success = 0
    total = 0
    if w_uni > 1:
        num = flip * 2 * w_uni
        total += num
        success += _flip_edges(rnd_graph, num, sets[0], _switch_single)
    if w_bi > 1:
        num = flip * w_bi # only single multiplication as edges are twice in
        # group
        total += num
        success += _flip_edges(rnd_graph, num, sets[1], _switch_double)
    return (rnd_graph, (float(success) / float(total)))

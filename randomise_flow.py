#!/usr/bin/env python
# encoding: utf-8


"""
randomise_flow_networks.py

Created by Moritz Beber on 2010-02-03.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""


from mfinder_wrapper import OptionsManager
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

def _make_groups(graph, nodes_in = 8, nodes_middle = 20, nodes_out = 8):
    nodes_out_begin = nodes_in + nodes_middle
    total = nodes_in + nodes_middle + nodes_out
    
    # simple definitions of categories by node index
    def _is_input(node):
        return (node <= nodes_in)

    def _is_middle(node):
        return (node > nodes_in and node <= nodes_out_begin)

    def _is_output(node):
        return (node > nodes_out_begin and node <= total)

    im = list()
    mm = list()
    mmd = list()
    mo = list()
    for edge in graph.edges_iter():
        if _is_input(edge[0]):
            if _is_middle(edge[1]):
                im.append(edge)
            else:
                networkx.exception.NetworkXError("Edge violates topology!")
        elif _is_middle(edge[0]):
            if _is_middle(edge[1]):
                if graph.has_edge(edge[1], edge[0]):
                    mmd.append(edge)
                else:
                    mm.append(edge)
            elif _is_output(edge[1]):
                mo.append(edge)
            else:
                networkx.exception.NetworkXError("Edge violates topology!")
        else:
            networkx.exception.NetworkXError("Edge violates topology!")
    return (im, mm, mmd, mo)

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
'randomise_bipartite' produces a randomised version of the supplied bipartite
graph 'graph'. This is achieved by switching edges in the graph a number of
times equal to the number of edges present times 'flip'. This function is
intended as a basis for calculating three-node subgraph statistics in metabolic
networks. As such only degrees of nodes, bi-directional link properties, etc.
(see:
R Milo, S Shen-Orr, S Itzkovitz, N Kashtan, D Chklovskii & U Alon,
Network Motifs: Simple Building Blocks of Complex Networks
Science, 298:824-827 (2002)
) are conserved. For larger subgraph statistics also the smaller subgraph
statistics would have to be conserved. The novelty of this function is to
preserve the bipartite nature of the graph. The randomised version is returned.
    """
    options = OptionsManager()
    rnd_graph = graph.copy()
    success = 0
    sets = _make_groups(rnd_graph)
    w_im = len(sets[0])
    w_mm = len(sets[1])
    w_mmd = len(sets[2])
    w_mo = len(sets[3])
    assert w_im and w_mm and w_mo, "Not all network layers are connected!"
    success = 0
    total = 0
    if w_im > 1:
        num = flip * 2 * w_im
        total += num
        success += _flip_edges(rnd_graph, num, sets[0], _switch_single)
    if w_mm > 1:
        num = flip * 2 * w_mm
        total += num
        success += _flip_edges(rnd_graph, num, sets[1], _switch_single)
    if w_mmd > 2:
        num = flip * w_mmd
        total += num
        success += _flip_edges(rnd_graph, num, sets[2], _switch_double)
    if w_mo > 1:
        num = flip * 2 * w_mo
        total += num
        success += _flip_edges(rnd_graph, num, sets[3], _switch_single)
    options.logger.info("Flip success rate: %f", (float(success) / float(total)))
    return rnd_graph

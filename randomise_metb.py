#!/usr/bin/env python
# encoding: utf-8


"""
metb_zscores.py

Created by Moritz Beber on 2009-09-07.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""


import random


m_id = "M"
r_id = "R"
rev_id = "_Rev"


def _add_edge(graph, src, tar, bunch):
    graph.add_edge(src, tar)
    bunch.append((src, tar))

def _remove_edge(graph, src, tar, bunch):
    graph.remove_edge(src, tar)
    bunch.remove((src, tar))

def _get_partners(edge):
    # substrate -> reaction or product -> reverse reaction
    if edge[0].startswith(m_id):
        met = edge[0]
        if edge[1].endswith(rev_id):
            r_rxn = edge[1]
            f_rxn = edge[1].replace(rev_id, "")
        else:
            f_rxn = edge[1]
            r_rxn = edge[1] + rev_id
    # reaction -> product or reverse reaction -> substrate
    else:
        met = edge[1]
        if edge[0].endswith(rev_id):
            r_rxn = edge[0]
            f_rxn = edge[0].replace(rev_id, "")
        else:
            f_rxn = edge[0]
            r_rxn = edge[0] + rev_id
    return (met, f_rxn, r_rxn)

def _switch_reversible_substrates(graph, first, second, bunch):
    if first == second:
        return False
    (met1, f_rxn1, r_rxn1) = _get_partners(first)
    (met2, f_rxn2, r_rxn2) = _get_partners(second)
    if graph.has_edge(met1, f_rxn2) or graph.has_edge(r_rxn2, met1):
        return False
    if graph.has_edge(met2, f_rxn1) or graph.has_edge(r_rxn1, met2):
        return False
    # add the forward direction
    _add_edge(graph, met2, f_rxn1, bunch)
    _add_edge(graph, met1, f_rxn2, bunch)
    # add the reverse direction
    _add_edge(graph, r_rxn1, met2, bunch)
    _add_edge(graph, r_rxn2, met1, bunch)
    # remove old edges
    _remove_edge(graph, met1, f_rxn1, bunch)
    _remove_edge(graph, r_rxn1, met1, bunch)
    _remove_edge(graph, met2, f_rxn2, bunch)
    _remove_edge(graph, r_rxn2, met2, bunch)
    return True

def _switch_reversible_products(graph, first, second, bunch):
    if first == second:
        return False
    (met1, f_rxn1, r_rxn1) = _get_partners(first)
    (met2, f_rxn2, r_rxn2) = _get_partners(second)
    if graph.has_edge(met1, r_rxn2) or graph.has_edge(f_rxn2, met1):
        # double test, one should suffice since we always switch as a unit
        return False
    if graph.has_edge(met2, r_rxn1) or graph.has_edge(f_rxn1, met2):
        return False
    # add the forward direction
    _add_edge(graph, f_rxn1, met2, bunch)
    _add_edge(graph, f_rxn2, met1, bunch)
    # add the reverse direction
    _add_edge(graph, met2, r_rxn1, bunch)
    _add_edge(graph, met1, r_rxn2, bunch)
    # remove old edges
    _remove_edge(graph, f_rxn1, met1, bunch)
    _remove_edge(graph, met1, r_rxn1, bunch)
    _remove_edge(graph, f_rxn2, met2, bunch)
    _remove_edge(graph, met2, r_rxn2, bunch)
    return True

def _switch_irreversible(graph, first, second, bunch):
    if first == second:
        return False
    # check if we would flip to existing edges
    if graph.has_edge(first[0], second[1]):
        return False
    if graph.has_edge(second[0], first[1]):
        return False
    # introduction of double edge is impossible due to categories
    _add_edge(graph, first[0], second[1], bunch)
    _add_edge(graph, second[0], first[1], bunch)
    _remove_edge(graph, first[0], first[1], bunch)
    _remove_edge(graph, second[0], second[1], bunch)
    return True

def _make_groups(graph):
    si = list()
    pi = list()
    sr = list()
    pr = list()
    for edge in graph.edges_iter():
        if edge[0].startswith(m_id):
            if edge[1].endswith(rev_id):
                pr.append(edge)
            elif graph.has_edge(edge[1] + rev_id, edge[0]):
                sr.append(edge)
            else:
                si.append(edge)
        else:
            if edge[0].endswith(rev_id):
                sr.append(edge)
            elif graph.has_edge(edge[1], edge[0] + rev_id):
                pr.append(edge)
            else:
                pi.append(edge)
    return (si, pi, sr, pr)

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
    rnd_graph = graph.copy()
    success = 0
    total = 0
    sets = _make_groups(rnd_graph)
    w_si = len(sets[0])
    w_pi = len(sets[1])
    w_sr = len(sets[2])
    w_pr = len(sets[3])
    assert w_si and w_pi, "Error: No irreversible reactions recognised!"
    assert w_sr and w_pr, "Error: No reversible reactions recognised!"
    if w_si > 1:
        num = flip * 2 * w_si
        total += num
        success += _flip_edges(rnd_graph, num, sets[0], _switch_irreversible)
    if w_pi > 1:
        num = flip * 2 * w_pi
        total += num
        success += _flip_edges(rnd_graph, num, sets[1], _switch_irreversible)
    if w_sr > 2:
        num = flip * w_sr
        total += num
        success += _flip_edges(rnd_graph, num, sets[2],
            _switch_reversible_substrates)
    if w_pr > 1:
        num = flip * w_pr
        total += num
        success += _flip_edges(rnd_graph, num, sets[3],
            _switch_reversible_products)
    print "Flip success rate: %f" % (float(success) / float(total))
    return rnd_graph

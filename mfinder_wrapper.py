#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===============
Mfinder Wrapper
===============

:Author:
    Moritz Emanuel Beber
:Date:
    2011-03-06
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    mfinder_wrapper.py
"""


import sys
import itertools
import os
import mfinder.mfinder as mfinder
import networkx as nx
import random
import numpy
import logging

from meb.utils.classes import BasicOptionsManager
from collections import defaultdict


class MfinderWrapper(object):
    """
    """

    def __init__(self, graph, subgraph_size=3, random_samples=100):
        object.__init__(self)
        self.logger = logging.getLogger("MfinderWrapper")
        self.current_dir = os.getcwd()
        self.graph = graph
        self.name2num = None
        self.num2name = None
        self.mtf_sz = subgraph_size
        self.incl_sym = False
        self.mtf_uniq = None
        self.mtf_counts = None
        self.members = None
        self.rewire = None
        self.rnd_num = random_samples
        self.rnd_counts = None
        self.rnd_graphs = None
        self.zscores = None

    def _make_numbering(self):
        """
        """
        self.name2num = dict()
        self.num2name = dict()
        for (name, num) in itertools.izip(self.graph.nodes_iter(),
                itertools.count(1)):
            self.name2num[name] = num
            self.num2name[num] = name

    def subgraph_adjacency(self, mtf_id):
        """
        """
        mtf_adj = nx.DiGraph()
        # infinitely increasing index 'i'
        for i in itertools.count():
            # if the bits in 'mtf_id' have been exhausted, exit loop
            if not mtf_id:
                break
            # add a link if the least significant bit is '1'
            elif mtf_id & 1:
                # in adjacency matrix row index is the source node
                # and column is the target node
                mtf_adj.add_edge(i // self.mtf_sz, i % self.mtf_sz)
            # remove the least significant bit by shifting to the right
            mtf_id >>= 1
        return mtf_adj

    def _subgraph_symmetries(self, subgraph):
        """
        This is a simple version of finding all symmetries of a given subgraph.
        This method only work when the size of the subgraph is equal to the
        subgraph_size and if the nodes are labelled 0, ..., N.
        """
        return [perm for perm in itertools.permutations(xrange(self.mtf_sz)) if
                all(subgraph.has_edge(perm[src], perm[tar]) for (src, tar) in
                    subgraph.edges_iter())]

    def _check_tuple(self, k_tuple, subgraph, symmetries):
        """
        """
        return [[k_tuple[i] for i in perm] for perm in symmetries if
                all(self.graph.has_edge(k_tuple[perm[src]], k_tuple[perm[tar]])
                    for (src, tar) in subgraph.edges_iter())]

    # wraps up edge consistency check for each ID found
    def ensure_consistency(self, select=random.choice):
        """
        The function select is expected to return one element of an iterable.
        """
        new_members = defaultdict(list)
        for (mtf_id, tuples) in self.members.iteritems():
            mtf_adj = self.subgraph_adjacency(mtf_id)
            symmetries = self._subgraph_symmetries(mtf_adj)
            for k_tuple in tuples:
                # generate a list of node tuples out of 'k_tuple' that are
                # consistent with the network structure
                tmp_tuples = [[k_tuple[k] for k in perm] for perm in symmetries
                        if all(self.graph.has_edge(k_tuple[perm[src]],
                        k_tuple[perm[tar]]) for (src, tar) in mtf_adj.edges_iter())]
                length = len(tmp_tuples)
                if length == 0:
                    self.logger.error("No valid permutation for %s in motif ID %d!",\
                        str(k_tuple), mtf_id)
                    continue
                elif (not self.incl_sym) and length > 1:
                    # we do not include symmetries, choose only one
                    new_members[mtf_id].append(select(tmp_tuples))
                else:
                    new_members[mtf_id].extend(tmp_tuples)
        self.members = new_members

    # translates every vertex tuple for each motif type back into real names
    def _reverse_translate(self):
        for (mtf_id, subgraphs) in self.members.iteritems():
            for member in subgraphs:
                for (index, vertex) in enumerate(member):
                    member[index] = self.num2name[vertex]

    def _make_mfinder_network(self, graph=None):
        if not graph:
            graph = self.graph
        links = mfinder.IntArray(graph.size() * 2)
        for (i, node) in itertools.izip(itertools.count(),
                itertools.chain.from_iterable(graph.edges_iter())):
            links[i] = self.name2num[node]
        return links

    def _extract_real_motifs(self, results):
        """
        """
        mtf_item = results.real.l
        members_arr = mfinder.IntArray(self.mtf_sz)
        self.mtf_counts = dict()
        self.members = dict()
        self.mtf_uniq = dict()
        while (mtf_item != None):
            mtf = mfinder.get_motif(mtf_item.p)
            mtf_id = int(mtf.id)
            self.mtf_counts[mtf_id] = int(mtf.count)
            self.mtf_uniq[mtf_id] = int(mtf.members.size)
            self.members[mtf_id] = list()
            member_item = mtf.all_members.l
            while (member_item != None):
                mfinder.get_members(member_item.p, members_arr, self.mtf_sz)
                self.members[mtf_id].append([int(members_arr[i]) for i in
                    xrange(self.mtf_sz)])
                member_item = member_item.next
            mtf_item = mtf_item.next

    def _extract_rnd_motifs(self, result):
        """
        """
        mtf_item = result.real.l
        while (mtf_item != None):
            mtf = mfinder.get_motif(mtf_item.p)
            self.rnd_counts[int(mtf.id)].append(int(mtf.count))
            mtf_item = mtf_item.next

    def count_subgraphs(self):
        """
        """
        self._make_numbering()
        links = self._make_mfinder_network()
        results = mfinder.subgraphs_interface(links, self.graph.size(),\
            self.mtf_sz, int(min(sys.maxint, len(self.graph) ** self.mtf_sz)))
        self._extract_real_motifs(results)
        mfinder.res_tbl_mem_free(results)
        self._reverse_translate()
        self.ensure_consistency()

    def compute_zscores(self):
        """
        """
        if not self.mtf_counts:
            self.count_subgraphs()
        self.rnd_counts = defaultdict(list)
        self.rnd_graphs = list()
        for i in xrange(self.rnd_num):
            self.logger.info("Generating %dth  random graph...", i + 1)
            (rnd_graph, success) = self.rewire.randomise(self.graph, flip=100)
            self.rnd_graphs.append(rnd_graph)
            self.logger.info("Flip success rate: %f", success)
            links = self._make_mfinder_network(rnd_graph)
            results = mfinder.subgraphs_interface(links, rnd_graph.size(),\
                self.mtf_sz, int(min(sys.maxint, len(rnd_graph) ** self.mtf_sz)))
            self._extract_rnd_motifs(results)
            mfinder.res_tbl_mem_free(results)
        self.zscores = dict()
        total = float(sum(self.mtf_counts.itervalues()))
        for mtf_id in self.mtf_counts:
            mtf_num = float(self.mtf_counts[mtf_id])
            if self.rnd_counts.has_key(mtf_id):
                rnd_mean = numpy.mean(self.rnd_counts[mtf_id])
                rnd_std = numpy.std(self.rnd_counts[mtf_id])
            else:
                rnd_mean = 0.0
                rnd_std = 0.0
            if rnd_std == 0.0:
                z_score = numpy.nan
            else:
                z_score = (mtf_num - rnd_mean) / rnd_std
            c_real = mtf_num / total * 1000.0
            self.zscores[mtf_id] = [self.mtf_counts[mtf_id],\
                rnd_mean, rnd_std, z_score, 0., self.mtf_uniq[mtf_id], c_real]

    def compute_zscores_from_extern(self, original_graph, rnd_graphs):
        """
        Requires mtf_counts, mtf_uniq, and rnd_counts to be set.
        """
        self.graph = original_graph
        self.count_subgraphs()
        self.rnd_graphs = rnd_graphs
        self.rnd_counts = defaultdict(list)
        for rnd_graph in self.rnd_graphs:
            links = self._make_mfinder_network(rnd_graph)
            results = mfinder.subgraphs_interface(links, rnd_graph.size(),\
                self.mtf_sz, min(sys.maxint, len(rnd_graph) ** self.mtf_sz))
            self._extract_rnd_motifs(results)
            mfinder.res_tbl_mem_free(results)
        self.zscores = dict()
        total = float(sum(self.mtf_counts.itervalues()))
        for mtf_id in self.mtf_counts:
            mtf_num = float(self.mtf_counts[mtf_id])
            if self.rnd_counts.has_key(mtf_id):
                rnd_mean = numpy.mean(self.rnd_counts[mtf_id])
                rnd_std = numpy.std(self.rnd_counts[mtf_id])
            else:
                rnd_mean = 0.0
                rnd_std = 0.0
            if rnd_std == 0.0:
                z_score = numpy.nan
            else:
                z_score = (mtf_num - rnd_mean) / rnd_std
            c_real = mtf_num / total * 1000.0
            self.zscores[mtf_id] = [self.mtf_counts[mtf_id],\
                rnd_mean, rnd_std, z_score, 0., self.mtf_uniq[mtf_id], c_real]


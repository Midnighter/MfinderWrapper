#!/usr/bin/env python
# encoding: utf-8


"""
mfinder_wrapper.py

Created by Moritz Beber on 2008-10-18.
Copyright (c) 2008 Jacobs University Bremen. All rights reserved.
"""


import itertools
import os
import logging
import sys
import getopt
import re
import errno
import mfinder.mfinder
import networkx
import random
import numpy


# version information
__version__ = "2.0"


class ArgumentError(StandardError):
    """
    docstring
    """
    def __init__(self, msg, opt=None, arg=None, *args, **kwargs):
        super(ArgumentError, self).__init__(*args, **kwargs)
        self.args = tuple([msg, opt, arg])
        self.errno = errno.EINVAL
        if opt and arg:
            self.strerror = msg % (opt, arg)
        elif opt:
            self.strerror = msg % opt
        elif arg:
            self.strerror = msg % arg
        else:
            self.strerror = msg

    def __str__(self):
        return self.strerror


class FormatError(EnvironmentError):
    """
    docstring
    """
    def __init__(self, line_no, line, filename, *args, **kwargs):
        super(FormatError, self).__init__(*args, **kwargs)
        self.args = tuple([line_no, line, filename])
        self.errno = errno.ENOEXEC
        self.line_no = int(line_no)
        self.strerror = line
        self.filename = filename

    def __str__(self):
        msg = "Bad format at line %d in file '%s'!\nContent: '%s'."\
            % (self.line_no, self.filename, self.strerror)
        return msg


class OptionsManager(object):
    """
    This class unifies some global options. This class is modelled
    as a singleton, i.e., only one instance of this class may exist globally.
    @todo: Consider thread safety, though I would defer that to using this class.
    """
    _singleton = None

    def __new__(cls, *args, **kwargs):
        if not cls._singleton:
            return super(OptionsManager, cls).__new__(cls, *args, **kwargs)
        return cls._singleton

    def __init__(self, *args, **kwargs):
        if self.__class__._singleton:
            return None
        super(OptionsManager, self).__init__(*args, **kwargs)
        self.log_levels = {"debug": logging.DEBUG, "info": logging.INFO,
            "warning": logging.WARNING, "error": logging.ERROR,
            "critical": logging.CRITICAL}
        formatter = logging.Formatter("%(message)s")
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(formatter)
        self.logger = logging.getLogger("MfinderWrapper")
        self.logger.propagate = False
        self.logger.setLevel(logging.ERROR)
        self.logger.addHandler(handler)
        self.current_dir = os.getcwd()
        self.source_dir = ""
        self.target_dir = ""
        self.mtf_sz = 3
        self.incl_sym = False
        self.numbering = None
        self.rnd_method = None
        self.rnd_num = 100
        self.post = False
        self.pattern = None
        self.only_members = False
        self.members_out = True
        self.__class__._singleton = self


# The first column is taken to be a unique identifier, following columns are
# attributed to the 'key', so any rows that have the same value in the first
# column will be replace previous versions
def parse_numbering(input_filename, sep=None):
    input_file = open(input_filename,"r")
    tmp = input_file.readlines()
    input_file.close()
    content = dict()
    for (index, line) in enumerate(tmp):
        if line.startswith("\n") or line.startswith("#"):
            continue
        ltmp = line.split(sep, 2)
        if len(ltmp) >= 2:
            try:
                content[ltmp[0]] = int(ltmp[1])
            except ValueError:
                raise FormatError(index, line, input_filename)
        else:
            raise FormatError(index, line, input_filename)
    return content

# makes a list of strings of each row and appends that list
def parse2list(input_filename, cols, sep=None):
    input_file = open(input_filename,"r")
    tmp = input_file.readlines()
    input_file.close()
    content = list()
    for (index, line) in enumerate(tmp):
        if line.startswith("\n") or line.startswith("#"):
            continue
        ltmp = line.split(sep, cols)
        if len(ltmp) >= cols:
            content.append(ltmp[:cols])
        else:
            raise FormatError(index, line, input_filename)
    return content

# makes a tuple of strings of each row and appends that list
def parse2tuples(input_filename, cols, sep=None):
    options = OptionsManager()
    input_file = open(input_filename,"r")
    tmp = input_file.readlines()
    input_file.close()
    content = list()
    for (index, line) in enumerate(tmp):
        if line.startswith("\n") or line.startswith("#"):
            continue
        ltmp = line.split(sep, cols)
        if len(ltmp) >= cols:
            content.append(tuple(ltmp[:cols]))
        else:
            error = FormatError(index, line, input_filename)
            options.logger.critical(str(error))
            raise error
    return content

# assigns an unique index to each ID appearing in the numbering file as required
# by mfinder later for the analysis
def assign_numbering(filename):
    name2num = parse_numbering(filename)
    num2name = dict()
    for (name, num) in name2num.iteritems():
        num2name[num] = name
    return (name2num, num2name)

# assigns an unique index to each vertex appearing in the edge file as required
# by mfinder later for the analysis
# already numbered vertices will still be regarded as string names and most
# likely reassigned
def make_numbering(edges):
    name2num = dict()
    num2name = dict()
    for (name, num) in itertools.izip(graph.nodes_iter(), itertools.count(1)):
        name2num[name] = num
        num2name[num] = name
    return (name2num, num2name)

# attempt to write a list of strings to a file, strings should be newline
# terminated
def write_lines(lines, filename):
    output = open(filename, 'w')
    output.writelines(lines)
    output.close()

def write_zscores(mtf_counts, mtf_uniq, rnd_counts, filename):
    out = open(filename, 'w')
#    ids = mtf_counts.keys()
#    ids += filter(lambda x: x not in mtf_counts.keys(), rnd_counts.keys())
    total = 0
    for mtf_id in mtf_counts:
        total += mtf_counts[mtf_id]
    for mtf_id in mtf_counts:
        rnd_mean = numpy.average(rnd_counts[mtf_id])
        rnd_std = numpy.std(rnd_counts[mtf_id])
        if rnd_std == 0.:
            z_score = 888888.
        else:
            z_score = (float(mtf_counts[mtf_id]) - rnd_mean) / rnd_std
        c_real = float(mtf_counts[mtf_id]) / float(total) * 1000.
        out.write("%d %d %f %f %f %f %d %f\n" % (mtf_id, mtf_counts[mtf_id],\
            rnd_mean, rnd_std, z_score, 0., mtf_uniq[mtf_id], c_real))
    out.close()

def write_members(members, filename):
    out = open(filename, 'w')
    for ids in members:
        out.write("# id %s\n" % ids)
        out.write("# number %s\n" % str(len(members[ids])))
        for v_tuple in members[ids]:
            out.write("\t".join(v_tuple) + "\n")
    out.close()

# parses a _MEMBERS.txt file and returns its contents sorted by ID
def parse_members(filename):
    options = OptionsManager()
    # read in file
    file = open(filename, "r")
    content = file.readlines()
    file.close()
    members = dict()
    motif_id = 0
    pattern = re.compile("(\\d+)")
    for (index, line) in enumerate(content):
        if line.startswith("\n") or line.startswith("#") or line.startswith("\r\n"):
            continue
        elif line.startswith("subgraph"):
            mobj = pattern.search(line)
            if mobj:
                motif_id = int(mobj.group(0))
                members[motif_id] = list()
        elif line.startswith("Nreal"):
            continue
        elif line.startswith("="):
            continue
        elif line.startswith("Partial"):
            options.logger.warning("Parsing only partial list for subgraph id"\
                " %d in file '%s'!", motif_id, filename)
        elif line.startswith("Full"):
            continue
        else:
            ltmp = line.split()
            if len(ltmp) >= 2 and len(ltmp) <= 6:
                motif = []
                for node in ltmp:
                    try:
                        motif.append(int(node))
                    except ValueError:
                        raise FormatError(index, line, filename)
                members[motif_id].append(motif)
            else:
                raise FormatError(index, line, filename)
    return members

# translates every vertex tuple for each motif type back into real names
def reverse_translate(members, num2name):
    for (mtf_id, subgraphs) in members.iteritems():
        for member in subgraphs:
            for (index, vertex) in enumerate(member):
                member[index] = num2name[vertex]

def find_symmetries(mtf_adj):
    options = OptionsManager()
    return [triple for triple in itertools.permutations(xrange(options.mtf_sz))\
            if all(adj.has_edge(triple[src], triple[tar])\
                    for (src, tar) in adj.edges_iter())]

def motif_adjacency(mtf_id, mtf_sz):
    """
Very fast - due to python, platform-independent - way of computing the motif
adjacency matrix from the binary representation of its integer ID.
    """
    mtf_adj = networkx.DiGraph()
    # infinitely increasing index 'i'
    for i in itertools.count():
        # if the bits in 'mtf_id' have been exhausted, exit loop
        if not mtf_id:
            break
        # add a link if the least significant bit is '1'
        if mtf_id & 1:
            # in adjacency matrix row index is the source node
            # and column is the target node
            mtf_adj.add_edge(i // mtf_sz, i % mtf_sz)
        # remove the least significant bit by shifting to the right
        mtf_id >>= 1
    return mtf_adj

# wraps up edge consistency check for each ID found
def ensure_consistency(members, network):
    options = OptionsManager()
    new_members = dict()
    setdef = new_members.setdefault
    for (mtf_id, tuples) in members.iteritems():
        mtf_adj = motif_adjacency(mtf_id, options.mtf_sz)
        for k_tuple in tuples:
            # generate a list of node tuples out of 'k_tuple' that are
            # consistent with the network structure
            tmp_tuples = [perm for perm in itertools.permutations(k_tuple) if\
                    all(network.has_edge(perm[src], perm[tar])\
                    for (src, tar) in mtf_adj.edges_iter())]
            length = len(tmp_tuples)
            if length == 0:
                options.logger.error("No valid permutation for %s in motif ID %d!",\
                    str(k_tuple), mtf_id)
                continue
            elif (not options.incl_sym) and length > 1:
                # we do not include symmetries, choose only one
                #setdef(mtf_id, []).append(tmp_tuples[0])
                # this could also be done by a random process
                setdef(mtf_id, []).append(random.choice(tmp_tuples))
            else:
                setdef(mtf_id, []).extend(tmp_tuples)
    return new_members

def make_mfinder_network(network):
    links = mfinder.mfinder.IntArray(network.size()* 2)
    index = 0
    for (src, tar) in network.out_edges_iter():
        links[index] = src
        index += 1
        links[index] = tar
        index += 1
    return links

def extract_rnd_motifs(result, mtf_counts):
    mtf_item = result.real.l
    while (mtf_item != None):
        mtf = mfinder.mfinder.get_motif(mtf_item.p)
        mtf_counts.setdefault(int(mtf.id), []).append(int(mtf.count))
        mtf_item = mtf_item.next

def extract_real_motifs(results):
    options = OptionsManager()
    mtf_item = results.real.l
    members = mfinder.mfinder.IntArray(options.mtf_sz)
    mtf_counts = dict()
    members_list = dict()
    mtf_uniq = dict()
    while (mtf_item != None):
        mtf = mfinder.mfinder.get_motif(mtf_item.p)
        mtf_id = int(mtf.id)
        mtf_counts[mtf_id] = int(mtf.count)
        mtf_uniq[mtf_id] = int(mtf.members.size)
        members_list[mtf_id] = []
        member_item = mtf.all_members.l
        while (member_item != None):
            mfinder.mfinder.get_members(member_item.p, members, options.mtf_sz)
            tmp = list()
            for i in xrange(options.mtf_sz):
                tmp.append(int(members[i]))
            members_list[mtf_id].append(tmp)
            member_item = member_item.next
        mtf_item = mtf_item.next
    return (mtf_counts, mtf_uniq, members_list)

def post_file(filename):
    raise NotImplementedError
#    options = OptionsManager()
#    base_file = filename.split(".", 1)[0] + "_mtf_sz%d" % options.mtf_sz
#    base_file = os.path.join(options.target_dir, base_file)
#    print "Starting reverse translation for graph in file '%s'..." % filename
#    edges = parse2tuples(os.path.join(options.source_dir, filename), 2)
#    network = networkx.DiGraph()
#    network.add_edges_from(edges)
#    # parse members file, must not throw
#    members = post_proc(base_file + "_MEMBERS.txt", num2name,\
#        links, mtf_size, sym)
#    m = list()
#    for ids in members:
#        m.append("# id %s\n" % ids)
#        m.append("# number %s\n" % str(len(members[ids])))
#        for v_tuple in members[ids]:
#            m.append("%s\t%s\t%s\n" % (v_tuple[0], v_tuple[1], v_tuple[2]))
#    outfile = filename.rsplit('.', 1)[0] + "_members.tsv"
#    outfile = os.path.join(target, outfile)
#    write_lines(m, outfile)
#    print "Done."

def proc_members(filename):
    # attempt preparation of mfinder analysis
    options = OptionsManager()
    base_file = filename.split(".", 1)[0] + "_mtf_sz%d" % options.mtf_sz
    base_file = os.path.join(options.target_dir, base_file)
    # parse the edges
    edges = parse2tuples(os.path.join(options.source_dir, filename), 2)
    # make an mfinder conform numbering of vertices, no exceptions expected
    if not options.numbering:
        options.numbering = make_numbering(edges)
    network = networkx.DiGraph()
    name2num = options.numbering[0]
    for edge in edges:
        network.add_edge(name2num[edge[0]], name2num[edge[1]])
    links = make_mfinder_network(network)
    results = mfinder.mfinder.subgraphs_interface(links, network.size(),\
        options.mtf_sz, len(network)**options.mtf_sz)
    (mtf_counts, mtf_uniq, members) = extract_real_motifs(results)
    mfinder.mfinder.res_tbl_mem_free(results)
    if options.members_out:
        members = ensure_consistency(members, network)
        reverse_translate(members, options.numbering[1])
        write_members(members, base_file + "_members.tsv")

def process_network(graph):
    options = OptionsManager()
    if not options.numbering:
        options.numbering = make_numbering(graph.edges())
    name2num = options.numbering[0]
    new_graph = networkx.DiGraph()
    for (src, tar) in graph.edges_iter():
        new_graph.add_edge(name2num[src], name2num[tar])
    links = make_mfinder_network(new_graph)
    results = mfinder.mfinder.subgraphs_interface(links, new_graph.size(),\
            options.mtf_sz, len(new_graph)**options.mtf_sz)
    (mtf_counts, mtf_uniq, members) = extract_real_motifs(results)
    mfinder.mfinder.res_tbl_mem_free(results)
    rnd_counts = dict()
    for i in xrange(options.rnd_num):
        (rnd_graph, success) = options.rnd_method.randomise(new_graph, flip=100)
        links = make_mfinder_network(rnd_graph)
        results = mfinder.mfinder.subgraphs_interface(links, rnd_graph.size(),\
                options.mtf_sz, len(rnd_graph)**options.mtf_sz)
        extract_rnd_motifs(results, rnd_counts)
        mfinder.mfinder.res_tbl_mem_free(results)
    return (mtf_counts, rnd_counts)

def proc_file(filename):
    # attempt preparation of mfinder analysis
    options = OptionsManager()
    base_file = filename.split(".", 1)[0] + "_mtf_sz%d" % options.mtf_sz
    base_file = os.path.join(options.target_dir, base_file)
    # parse the edges
    edges = parse2tuples(os.path.join(options.source_dir, filename), 2)
    # make an mfinder conform numbering of vertices, no exceptions expected
    if not options.numbering:
        options.numbering = make_numbering(edges)
    network = networkx.DiGraph()
    name2num = options.numbering[0]
    for edge in edges:
        network.add_edge(name2num[edge[0]], name2num[edge[1]])
    links = make_mfinder_network(network)
    results = mfinder.mfinder.subgraphs_interface(links, network.size(),\
        options.mtf_sz, len(network)**options.mtf_sz)
    (mtf_counts, mtf_uniq, members) = extract_real_motifs(results)
    mfinder.mfinder.res_tbl_mem_free(results)
    if options.members_out:
        members = ensure_consistency(members, network)
        reverse_translate(members, options.numbering[1])
        write_members(members, base_file + "_members.tsv")
    rnd_counts = dict()
    for i in xrange(options.rnd_num):
        options.logger.info("Generating %dth  random graph...", i + 1)
        (rnd_graph, success) = options.rnd_method.randomise(network, flip=100)
        options.logger.info("Flip success rate: %f", success)
        links = make_mfinder_network(rnd_graph)
        results = mfinder.mfinder.subgraphs_interface(links, rnd_graph.size(),\
            options.mtf_sz, len(rnd_graph)**options.mtf_sz)
        extract_rnd_motifs(results, rnd_counts)
        mfinder.mfinder.res_tbl_mem_free(results)
    write_zscores(mtf_counts, mtf_uniq, rnd_counts, base_file + "_mat.tsv")

def process_file(filename):
    options = OptionsManager()
    if options.only_members:
        proc_members(filename)
    elif options.post:
        post_file(filename)
    else:
        proc_file(filename)

def write_usage(exe):
    usage = """
    mfinder_wrapper Version %s
    Created by Moritz Beber on 2008-10-18.
    Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.

    The mfinder_wrapper's principal task is to use mfinder's C functions to
    analyse the motif content of given networks, order the motifs found
    according to the links present in the network (and possibly its symmetric
    rotations), and compute motif statistics according to custom randomisation
    algorithms.

    mfinder_wrapper can be invoked in three different forms:

    1. Take any number of files which represent edgelists. Only two columns are
       are parsed. Empty lines and lines beginning with '#' are ignored.

       %s [options] --randomise <module> file1 [file2] [file3] [...]

    2. Perform the same tasks but on all files in a certain directory. This can
       be combined with a regular expression to restrict the files parsed in
       that directory.

       %s [options] [--regex <pattern>] --randomise <module> -d <path>

    3. From a given mfinder output, verify the node order within motifs and
       translate the file back into actual names. This requires the original
       vertex numbering used. This is compatible with either of the first two
       ways of finding edgelist files.

       %s [options] -n <file> -p file1 [file2] [file3] [...]
       or:
       %s [options] [--regex <pattern>] -n <file> -p -d <path>

    Options:

    -h, --help              print this information.

    -v, --version           print the version information.

    -d, --directory <path>  supply a directory that contains any number of
                            metabolite-centric network graphs only.

    -p, --post-translation <filename>

    -r, --regex <pattern>       supply a regular expression pattern to recognise
                                edgelist files that should be parsed

    --randomise <string>        python module to import that generates random
                                networks as a null model

    --rnd-num                   number of random networks considered for zscore
                                (default is 100)

    -n, --numbering <string>    provide a file with unique node IDs

    -i, --include-symmetric     choose to include all symmetric rotations of a
                                motif into the member file

    -m, --motif-size <int>      specify a motif size to be investigated (default 3),
                                mfinder can detect sizes between 2 and 6

    -t, --target <directory>    supply a directory where the wrapper should write
                                its output

    --no-members                suppress the output of the individual motif
                                members vertex list

    --only-members              generate a list of motif members only, no
                                randomisation is applied

    --log-level <level>         set the verbosity of information printed (default
                                is 'error' other options: 'debug', 'info',
                                'warning', 'error', 'critical')
    """ % (__version__, exe, exe, exe, exe)
    return usage


def main(argv, exe):
    options = OptionsManager()
    if not argv:
        raise ArgumentError(write_usage(exe))
    try:
        (opts, args) = getopt.getopt(argv, "hvd:pit:m:n:r:", ["help", "version",\
            "directory=", "post-translation", "include-symmetric",\
            "target=", "motif-size=", "numbering=", "regex=", "randomise=",\
            "no-members", "only-members", "rnd-num=", "log-level="])
    except getopt.GetoptError, err:
        raise ArgumentError(err.msg, err.opt)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            options.logger.critical(write_usage(exe))
            return 0
        elif opt in ("-v", "--version"):
            options.logger.critical("Version %s.", __version__)
            return 0
        elif opt in ("-d", "--directory"):
            tar = os.path.normpath(os.path.join(options.current_dir, arg))
            if os.path.isdir(tar):
                options.source_dir = tar
            else:
                error = ArgumentError("Option %s's argument %s is not a valid"\
                    " directory!", opt, arg)
                error.errno = errno.ENOENT
                raise error
        elif opt in ("-p", "--post-translation"):
            options.post = True
        elif opt in ("-i", "--include-symmetric"):
            options.incl_sym = True
        elif opt in ("-t", "--target"):
            tar = os.path.normpath(os.path.join(options.current_dir, arg))
            if os.path.isdir(tar):
                options.target_dir = tar
            else:
                error = ArgumentError("Option \"%s\"'s argument \"%s\" is not"\
                    " a valid directory!", opt, arg)
                error.errno = errno.ENOENT
                raise error
        elif opt in ("-m", "--motif-size"):
            try:
                options.mtf_sz = int(arg)
            except ValueError:
                raise ArgumentError("Option \"%s\"'s argument \"%s\" has"\
                    " the wrong value!", opt, arg)
        elif opt in ("-n", "--numbering"):
            tar = os.path.normpath(os.path.join(options.current_dir, arg))
            if os.path.isfile(tar):
                options.numbering = assign_numbering(tar)
            else:
                error = ArgumentError("Option \"%s\"'s argument \"%s\" is not"\
                    " a valid file!", opt, arg)
                error.errno = errno.ENOENT
                raise error
        elif opt in ("-r", "--regex"):
            options.pattern = re.compile(arg)
        elif opt in ("--randomise"):
            try:
                options.rnd_method = __import__(arg)
            except ImportError:
                raise ArgumentError("Option \"%s\"'s argument \"%s\" is not"\
                    " a python module that can be imported!", opt, arg)
        elif opt in ("--no-members"):
            options.members_out = False
        elif opt in ("--only-members"):
            options.only_members = True
        elif opt in ("--rnd-num"):
            try:
                options.rnd_num = int(arg)
            except ValueError:
                raise ArgumentError("Option \"%s\"'s argument \"%s\" has"\
                    " the wrong value!", opt, arg)
        elif opt in ("--log-level"):
            try:
                options.logger.setLevel(options.log_levels[arg])
            except KeyError:
                raise ArgumentError("Option \"%s\"'s argument \"%s\" has"\
                    " the wrong value!", opt, arg)
        else:
            raise ArgumentError("Invalid option \"%s\"!", opt)
    if options.source_dir and not args:
        if not options.target_dir:
            options.target_dir = options.source_dir
        for filename in os.listdir(options.source_dir):
            if os.path.isfile(os.path.join(options.source_dir, filename)):
                if options.pattern:
                    mobj = options.pattern.search(filename)
                    if not mobj:
                        options.logger.debug("'%s' does not match the regular"\
                            " expression!", filename)
                        continue
                options.logger.info("Working on '%s'...", filename)
                try:
                    process_file(filename)
                except StandardError, err:
                    options.logger.debug(str(err))
                    options.logger.error("Aborting motif analysis for the"\
                        " network in file '%s'!",
                        os.path.join(options.source_dir, filename))
                    continue
            else:
                options.logger.warning("'%s' is not a valid file!",\
                    os.path.join(options.source_dir, filename))
    elif args and not  options.source_dir:
        for filename in args:
            if (os.path.isfile(filename)):
                filename = os.path.normpath(os.path.join(options.current_dir,\
                    filename))
                (options.source_dir, filename) = os.path.split(filename)
                if options.pattern:
                    mobj = options.pattern.search(filename)
                    if not mobj:
                        options.logger.debug("'%s' does not match the regular"\
                            " expression!", filename)
                        continue
                if not options.target_dir:
                    options.target_dir = options.source_dir
                options.logger.info("Working on '%s'...", filename)
                try:
                    process_file(filename)
                except IOError, err:
                    options.logger.debug(str(err))
                    options.logger.error("Aborting motif analysis for the"\
                        " network in file '%s'!",
                        os.path.join(options.source_dir, filename))
                    continue
            else:
                options.logger.warning("'%s' is not a valid file!",\
                    os.path.join(options.source_dir, filename))
    else:
        raise ArgumentError("The combination of supplied options and their"\
            " arguments is not recognised!\n%s", " ".join(argv))

if __name__ == '__main__':
    options = OptionsManager()
    rc = 0
    try:
        rc = main(sys.argv[1:], sys.argv[0])
    except StandardError, err:
        options.logger.critical(str(err))
        options.logger.critical("Try %s --help", sys.argv[0])
        if hasattr(err, "errno"):
            rc = err.errno
        else:
            rc = 1
    finally:
        logging.shutdown()
        sys.exit(rc)


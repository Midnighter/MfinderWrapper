#!/usr/bin/env python
# encoding: utf-8


"""
mfinder_wrapper.py

Created by Moritz Beber on 2008-10-18.
Copyright (c) 2008 Jacobs University Bremen. All rights reserved.
"""


import os
import logging
import sys
import getopt
import re
import errno
import mfinder.mfinder
import networkx
#import pdb
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
        self.logger.setLevel(logging.DEBUG)
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
    index = 1
    for edge in edges:
        if not edge[0] in name2num:
            name2num[edge[0]] = index
            num2name[index] = edge[0]
            index += 1
        if not edge[1] in name2num:
            name2num[edge[1]] = index
            num2name[index] = edge[1]
            index += 1
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
        if line.startswith("\n"):
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
    for iden in members:
        for member in members[iden]:
            for (index, vertex) in enumerate(member):
                member[index] = num2name[vertex]

# ensure that a link required by the adjacency matrix of a motif is present
def verify_tuple(n_tuple, mtf_adj, network):
    rc = True
    for (src, tar) in mtf_adj.out_edges_iter():
        rc &= network.has_edge(n_tuple[src], n_tuple[tar])
    return rc

# gives all permutations of a list
# adopted from: http://snippets.dzone.com/posts/show/753
def all_permutations(a):
    if len(a) <=1:
        yield a
    else:
        for perm in all_permutations(a[1:]):
            for i in range(len(perm) + 1):
                yield perm[:i] + a[0:1] + perm[i:]

def test_permutations(mtf_perms, n_tuple, mtf_adj, network):
    options = OptionsManager()
    # tmp list to store permuted ordering of vertices
    tmp_members = list()
    tmp_tuple = list()
    # check all possible permutations
    for perm in mtf_perms:
        # build permuted motif members
        for index in perm:
            tmp_tuple.append(n_tuple[index])
        if verify_tuple(tmp_tuple, mtf_adj, network):
            tmp_members.append(tmp_tuple)
            tmp_tuple = []
            if not options.incl_sym:
                return tmp_members
        else:
            tmp_tuple = []
    return tmp_members

# verifies all occurences of a specific motif
def verify_links(iden, members, network, mtf_perms):
    options = OptionsManager()
    # motif adjacency
    mtf_adj = motif_adjacency(iden, options.mtf_sz)
    # tmp list to store permuted ordering of vertices
    tmp_tuple = list()
    tmp_members = list()
    # check each vertex tuple
    for n_tuple in members:
        tmp_tuple = test_permutations(mtf_perms, n_tuple, mtf_adj, network)
        if not tmp_tuple:
            options.logger.warning("No valid permutation for %s", str(n_tuple))
            continue
        tmp_members.extend(tmp_tuple)
    return tmp_members

# creates the adjacency list for any motif ID
# could be replaced by bin(id) if the width of binary string could be adjusted
# and most importantly the order of the binary string would begin with least
# significant digit
def motif_adjacency(iden, mtf_size):
    num = iden
    motif_adj = networkx.DiGraph()
    for i in xrange(mtf_size**2):
        if num % 2:
            motif_adj.add_edge(i / mtf_size, i % mtf_size)
        num = num>>1
    return motif_adj

# wraps up edge consistency check for each ID found
def ensure_consistency(members, network):
    options = OptionsManager()
    # column switches that may be tried
    # list of all permutations of indeces
    motif_perms = list()
    rotation = range(options.mtf_sz)
    for p in all_permutations(rotation):
        motif_perms.append(p)
    for iden in members:
        members[iden] = verify_links(iden, members[iden], network, motif_perms)

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
        ensure_consistency(members, network)
        reverse_translate(members, options.numbering[1])
        write_members(members, base_file + "_members.tsv")
    rnd_counts = dict()
    for i in xrange(options.rnd_num):
        options.logger.info("Generating %dth  random graph...", i + 1)
        rnd_graph = options.rnd_method.randomise(network, flip=100)
        links = make_mfinder_network(rnd_graph)
        results = mfinder.mfinder.subgraphs_interface(links, rnd_graph.size(),\
            options.mtf_sz, rnd_graph.size()**options.mtf_sz)
        extract_rnd_motifs(results, rnd_counts)
        mfinder.mfinder.res_tbl_mem_free(results)
    write_zscores(mtf_counts, mtf_uniq, rnd_counts, base_file + "_mat.tsv")

def process_file(filename):
    options = OptionsManager()
    if options.post:
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
            "no-members", "--rnd-num="])
    except getopt.GetoptError, err:
        raise ArgumentError(err.msg, err.opt)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            options.logger.critical(write_usage(exe))
            return None
        elif opt in ("-v", "--version"):
            options.logger.critical("Version %s.", __version__)
            return None
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
        elif opt in ("--rnd-num"):
            try:
                options.rnd_num = int(arg)
            except ValueError:
                raise ArgumentError("Option \"%s\"'s argument \"%s\" has"\
                    " the wrong value!", opt, arg)
        else:
            raise ArgumentError("Invalid option \"%s\"!", opt)
    if options.source_dir and not args:
        if not options.target_dir:
            options.target_dir = options.source_dir
        for filename in os.listdir(options.source_dir):
            if os.path.isfile(os.path.join(options.source_dir, filename)):
                try:
                    process_file(filename)
                except StandardError, err:
                    options.logger.debug(str(err))
                    options.logger.error("Aborting motif analysis for the"\
                        " network in file '%s'!",
                        os.path.join(options.source_dir, filename))
                    continue
            else:
                options.logger.error("'%s' is not a valid file!",\
                    os.path.join(options.source_dir, filename))
    elif args and not  options.source_dir:
        for filename in args:
            if (os.path.isfile(filename)):
                filename = os.path.normpath(os.path.join(options.current_dir,\
                    filename))
                (options.source_dir, filename) = os.path.split(filename)
                if not options.target_dir:
                    options.target_dir = options.source_dir
                try:
                    process_file(filename)
                except IOError, err:
                    options.logger.debug(str(err))
                    options.logger.error("Aborting motif analysis for the"\
                        " network in file '%s'!",
                        os.path.join(options.source_dir, filename))
                    continue
            else:
                options.logger.error("'%s' is not a valid file!",\
                    os.path.join(options.source_dir, filename))
    else:
        raise ArgumentError("The combination of supplied options and their"\
            " arguments is not recognised!\n%s", " ".join(argv))

if __name__ == '__main__':
#    pdb.set_trace()
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


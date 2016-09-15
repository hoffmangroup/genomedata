#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
load_seq: DESCRIPTION
"""

# Copyright 2008-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

from re import compile, VERBOSE
import sys
import warnings

from argparse import ArgumentParser

from numpy import frombuffer, uint32
from path import path
from tabdelim import DictReader

from . import (SEQ_ATOM, SEQ_DTYPE, FILE_MODE_CHROMS, 
              FORMAT_VERSION, Genome, __version__)
from ._util import FILTERS_GZIP, LightIterator, maybe_gzip_open

MIN_GAP_LEN = 100000
assert not MIN_GAP_LEN % 2 # must be even for division

# sre_constants.MAXREPEAT is 65535, so I have to break repeats into
# two segments
REGEX_SEGMENT_LEN = MIN_GAP_LEN // 2 # max == MAXREPEAT

DNA_LETTERS_UNAMBIG = "ACGTacgt"

SUPERCONTIG_NAME_FMT = "supercontig_%s"

AGP_FIELDNAMES = ["object", "object_beg", "object_end", "part_number",
                  "component_type", "col6", "col7", "col8", "col9"]

GAP_COMPONENT_TYPES = frozenset("NU")

GenomicPosition = uint32
GENOMIC_POSITION_0 = GenomicPosition(0)

def ignore_comments(iterable):
    return (item for item in iterable if not item.startswith("#"))

def read_agp(iterable):
    """
    converts coordinates from 1-based to 0-based
    """
    reader = DictReader(ignore_comments(iterable), AGP_FIELDNAMES)
    for row in reader:
        row["object_beg"] = int(row["object_beg"]) - 1
        row["object_end"] = int(row["object_end"])
        row["part_number"] = int(row["part_number"])

        if row["component_type"] in GAP_COMPONENT_TYPES:
            row["gap_length"] = int(row["col6"])
            row["gap_type"] = row["col7"]
            if row["col8"] == "yes":
                row["linkage"] = True
            else:
                assert row["col8"] == "no"
                row["linkage"] = False
            row["linkage_evidence"] = row["col9"]
        else:
            row["component_id"] = row["col6"]
            row["component_beg"] = int(row["col7"]) - 1
            row["component_end"] = int(row["col8"])
            row["orientation"] = row["col8"]

        yield row

def create_supercontig(chromosome, index, seq=None, start=None, end=None):
    name = SUPERCONTIG_NAME_FMT % index
    h5file = chromosome.h5file
    where = chromosome.h5group
    supercontig = h5file.create_group(where, name)

    if seq is not None:
        seq_array = frombuffer(seq, SEQ_DTYPE)
        h5file.create_carray(supercontig, "seq", SEQ_ATOM, seq_array.shape)

        # XXXopt: does this result in compression?
        supercontig.seq[...] = seq_array

    if start is None:
        start = chromosome.start
    if end is None:
        end = chromosome.end

    attrs = supercontig._v_attrs
    attrs.start = GenomicPosition(start)
    attrs.end = GenomicPosition(end)

# XXXopt: the all-regex approach is much slower than the hybrid
# approach (3 min to load chr21, 6 min to load chr1), but it is easier
# to understand the code and get it correct
#
# the previous code (r22) might have worked fine. Consider backing down to
# it at some point.
#
# XXXopt: a numpy implementation might be better
re_gap_segment = compile(r"""
(?:([^%s]{%d}[^%s]{%d,})                                  # group(1): ambig
   |                                                      #  OR
   ((?:(?:[%s]+|^)(?:[^%s]{1,%d}[^%s]{,%d}(?![^%s]))*)+)) # group(2): unambig
""" % (DNA_LETTERS_UNAMBIG, REGEX_SEGMENT_LEN,
       DNA_LETTERS_UNAMBIG, REGEX_SEGMENT_LEN,
       DNA_LETTERS_UNAMBIG, DNA_LETTERS_UNAMBIG, REGEX_SEGMENT_LEN,
       DNA_LETTERS_UNAMBIG, REGEX_SEGMENT_LEN-1,
       DNA_LETTERS_UNAMBIG), VERBOSE)

def init_chromosome_start(chromosome):
    """
    doesn't set end, so you can't use create_supercontig(end=None) if you use this
    """
    chromosome.attrs.start = GENOMIC_POSITION_0
    chromosome._file_attrs.genomedata_format_version = FORMAT_VERSION

def init_chromosome(chromosome, size):
    init_chromosome_start(chromosome)
    chromosome.attrs.end = GenomicPosition(size)

def read_seq(chromosome, seq):
    supercontig_index = 0

    init_chromosome(chromosome, len(seq))

    for m_segment in re_gap_segment.finditer(seq):
        seq_unambig = m_segment.group(2)
        if seq_unambig:
            span = m_segment.span()
            create_supercontig(chromosome, supercontig_index,
                               seq_unambig, *span)
            supercontig_index += 1
        else:
            assert m_segment.group(1)

def read_assembly(chromosome, infile):
    supercontig_index = 0

    init_chromosome_start(chromosome)

    for part in read_agp(infile):
        end = part["object_end"]

        if part["component_type"] in GAP_COMPONENT_TYPES:
            continue

        start = part["object_beg"]

        create_supercontig(chromosome, supercontig_index, None, start, end)
        supercontig_index += 1

    chromosome.attrs.end = GenomicPosition(end)

def size_chromosome(chromosome, size):
    init_chromosome(chromosome, size)
    create_supercontig(chromosome, 0)

def load_sizes(filename):
    res = {}

    with maybe_gzip_open(filename) as infile:
        for line in infile:
            row = line.rstrip().split()
            res[row[0]] = int(row[1])

    return res

def get_num_seq(filenames):
    """
    Run through all files once and count the number of sequences
    """
    res = 0

    for filename in filenames:
        with maybe_gzip_open(filename) as infile:
            for line in infile:
                if line.startswith(">"):
                    res += 1

    return res

def create_chromosome(genome, name, mode):
    name = "_".join(name.split())  # Remove any whitespace
    if mode == "dir":
        res = genome[name]
    else: # mode == "file"
        h5file = genome.h5file
        h5file.create_group("/", name, filters=FILTERS_GZIP)
        res = genome[name]

    res.attrs.dirty = True

    return res

def load_seq(gdfilename, filenames, verbose=False, mode=None, seqfile_type="fasta"):
    gdpath = path(gdfilename)

    ## load sizes if necessary to figure out number of chromosomes
    if seqfile_type == "sizes":
        assert len(filenames) == 1
        sizes = load_sizes(filenames[0])
    else:
        sizes = None

    ## decide whether should use dir or file mode
    if mode is None:
        if seqfile_type == "sizes":
            num_seq = len(sizes)
        elif seqfile_type == "agp":
            num_seq = len(filenames)
        else:
            num_seq = get_num_seq(filenames)

        if num_seq < FILE_MODE_CHROMS:
            mode = "dir"
        else:
            mode = "file"

        if verbose:
            msg = ("Implementation unspecified. Found %d"
                   "chromosomes/scaffolds, so using: %s" % (num_seq, mode))
            print(msg, file=sys.stderr)
    if mode == "dir":
        if gdpath.exists():
            assert gdpath.isdir()
        else:
            gdpath.makedirs()
    elif mode == "file":
        assert not gdpath.exists()
        gdpath.touch()  # Create file (then Genomedata archive will be a file)
    else:
        raise ValueError("Unsupported mode: %s" % mode)

    # XXX: ignoring all warnings is probably bad. Why is this here?
    # Can we be more specific?
    #
    # is this because we are going to close chromosomes when they are
    # set dirty?
    warnings.simplefilter("ignore")
    with Genome(gdpath, mode="w", filters=FILTERS_GZIP) as genome:
        if seqfile_type == "sizes":
            for name, size in sizes.items():
                chromosome = create_chromosome(genome, name, mode)
                size_chromosome(chromosome, size)
        else:
            assert seqfile_type in frozenset(["agp", "fasta"])
            for filename in filenames:
                if verbose:
                    print(filename, file=sys.stderr)

                with maybe_gzip_open(filename) as infile:
                    if seqfile_type == "agp":
                        name = path(filename).name.rpartition(".agp")[0]
                        chromosome = create_chromosome(genome, name, mode)
                        read_assembly(chromosome, infile)
                    else:
                        for defline, seq in LightIterator(infile):
                            chromosome = create_chromosome(genome, defline, mode)
                            read_seq(chromosome, seq)
    # XXX: this should be enforced even when there is an exception
    # is there a context manager available?
    warnings.resetwarnings()

def parse_options(args):

    description = ("Start a Genomedata archive at GENOMEDATAFILE with the"
                   " provided sequences. SEQFILEs should be in fasta format,"
                   " and a separate Chromosome will be created for each"
                   " definition line.")

    parser = ArgumentParser(description=description,
                            prog='genomedata-load-seq',
                            version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')

    parser.add_argument('seqfiles', nargs='+',
                        help='sequences in FASTA format')

    parser.add_argument("-a", "--assembly", action="store_const",
                        const="agp", dest="seqfile_type",
                        help="SEQFILE contains assembly (AGP) files instead of"
                        " sequenc")
    parser.add_argument("-s", "--sizes", action="store_const", const="sizes",
                        dest="seqfile_type", default="fasta",
                        help="SEQFILE contains list of sizes instead of"
                        " sequence")
    parser.add_argument("-f", "--file-mode", dest="mode",
                        default=None, action="store_const", const="file",
                        help="If specified, the Genomedata archive will be"
                        " implemented as a single file, with a separate h5 group"
                        " for each Chromosome. This is recommended if there are"
                        " a large number of Chromosomes. The default behavior is"
                        " to use a single file if there are at least %s"
                        " Chromosomes being added." % FILE_MODE_CHROMS)
    parser.add_argument("-d", "--directory-mode", dest="mode",
                        action="store_const", const="dir",
                        help="If specified, the Genomedata archive will be"
                        " implemented as a directory, with a separate file for"
                        " each Chromosome. This is recommended if there are a"
                        " small number of Chromosomes. The default behavior is"
                        " to use a directory if there are fewer than %s"
                        " Chromosomes being added." % FILE_MODE_CHROMS)
    parser.add_argument("--verbose",default=False, action="store_true",
                      help="Print status updates and diagnostic messages")

    args = parser.parse_args(args)

    return args

def main(argv=sys.argv[1:]):
    args = parse_options(argv)

    return load_seq(args.gdarchive, args.seqfiles, verbose=args.verbose,
                    mode=args.mode, seqfile_type=args.seqfile_type)

if __name__ == "__main__":
    sys.exit(main())

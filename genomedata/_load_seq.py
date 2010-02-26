#!/usr/bin/env python
from __future__ import division, with_statement

"""
load_seq: DESCRIPTION
"""

__version__ = "$Revision$"

# Copyright 2008-2009 Michael M. Hoffman <mmh1@washington.edu>

from re import compile, VERBOSE
import sys

from numpy import frombuffer
from path import path
import warnings

from . import SEQ_ATOM, SEQ_DTYPE, FILE_MODE_CHROMS, Genome
from ._util import FILTERS_GZIP, LightIterator, maybe_gzip_open

MIN_GAP_LEN = 100000
assert not MIN_GAP_LEN % 2 # must be even for division

# sre_constants.MAXREPEAT is 65535, so I have to break repeats into
# two segments
REGEX_SEGMENT_LEN = MIN_GAP_LEN // 2 # max == MAXREPEAT

DNA_LETTERS_UNAMBIG = "ACGTacgt"

SUPERCONTIG_NAME_FMT = "supercontig_%s"

def create_supercontig(chromosome, index, seq, start, end):
    name = SUPERCONTIG_NAME_FMT % index
    h5file = chromosome.h5file
    where = chromosome.h5group
    supercontig = h5file.createGroup(where, name)

    seq_array = frombuffer(seq, SEQ_DTYPE)
    h5file.createCArray(supercontig, "seq", SEQ_ATOM, seq_array.shape)
    supercontig.seq[...] = seq_array

    attrs = supercontig._v_attrs
    attrs.start = start
    attrs.end = end

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

def read_seq(chromosome, seq):
    supercontig_index = 0

    for m_segment in re_gap_segment.finditer(seq):
        seq_unambig = m_segment.group(2)
        if seq_unambig:
            span = m_segment.span()
            create_supercontig(chromosome, supercontig_index,
                               seq_unambig, *span)
            supercontig_index += 1
        else:
            assert m_segment.group(1)

def load_seq(gdfilename, filenames, verbose=False, mode=None):
    gdpath = path(gdfilename)

    if mode is None:
        # Run through all files once and count the number of sequences
        num_seq = 0
        for filename in filenames:
            with maybe_gzip_open(filename) as infile:
                for line in infile:
                    if line.startswith(">"):
                        num_seq += 1

        if num_seq < FILE_MODE_CHROMS:
            mode = "dir"
        else:
            mode = "file"

        if verbose:
            print >>sys.stderr, ("Implementation unspecified. Found %d"
                                 " chromosomes/scaffolds, so using: %s"
                                 % (num_seq, mode))

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

    warnings.simplefilter("ignore")
    with Genome(gdpath, mode="w", filters=FILTERS_GZIP) as genome:
        for filename in filenames:
            if verbose:
                print >>sys.stderr, filename

            with maybe_gzip_open(filename) as infile:
                for defline, seq in LightIterator(infile):
                    name = "_".join(defline.split())  # Remove any whitespace
                    if mode == "dir":
                        chromosome = genome[name]
                    else: # mode == "file"
                        h5file = genome._h5file
                        h5group = h5file.createGroup("/", name,
                                                     filters=FILTERS_GZIP)
                        chromosome = genome[name]
                    chromosome.attrs.dirty = True
                    read_seq(chromosome, seq)

    warnings.resetwarnings()

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... GENOMEDATAFILE SEQFILE..."
    version = "%%prog %s" % __version__
    description = ("Start a Genomedata archive at GENOMEDATAFILE with the"
                   " provided sequences. SEQFILEs should be in fasta format,"
                   " and a separate Chromosome will be created for each"
                   " definition line.")
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-v", "--verbose", dest="verbose",
                      default=False, action="store_true",
                      help="Print status updates and diagnostic messages")
    parser.add_option("-f", "--file-mode", dest="mode",
                      default=None, action="store_const", const="file",
                      help="If specified, the Genomedata archive will be"
                      " implemented as a single file, with a separate h5 group"
                      " for each Chromosome. This is recommended if there are"
                      " a large number of Chromosomes. The default behavior is"
                      " to use a single file if there are at least %s"
                      " Chromosomes being added." % FILE_MODE_CHROMS)
    parser.add_option("-d", "--directory-mode", dest="mode",
                      action="store_const", const="dir",
                      help="If specified, the Genomedata archive will be"
                      " implemented as a directory, with a separate file for"
                      " each Chromosome. This is recommended if there are a"
                      " small number of Chromosomes. The default behavior is"
                      " to use a directory if there are fewer than %s"
                      " Chromosomes being added." % FILE_MODE_CHROMS)

    options, args = parser.parse_args(args)

    if not len(args) >= 2:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    gdfilename = args[0]
    seqfiles = args[1:]
    kwargs = {"verbose": options.verbose,
              "mode": options.mode}
    return load_seq(gdfilename, seqfiles, **kwargs)

if __name__ == "__main__":
    sys.exit(main())

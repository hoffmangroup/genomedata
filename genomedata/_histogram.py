#!/usr/bin/env python
from __future__ import division, with_statement

"""
_histogram: prints histogram
"""

__version__ = "$Revision$"

# Copyright 2008, 2013 Michael M. Hoffman <mmh1@washington.edu>

import sys

from collections import defaultdict
from functools import partial
from numpy import (array, concatenate, histogram, iinfo, isfinite, ndarray,
                   NINF, PINF, zeros)

XXX all of these need to be removed
from segway._util import (DTYPE_IDENTIFY, fill_array, iter_chroms_coords,
                          load_coords, walk_continuous_supercontigs)

FIELDNAMES = ["lower_edge", "count"]

IINFO_IDENTIFY = iinfo(DTYPE_IDENTIFY)
MAX_IDENTIFY = IINFO_IDENTIFY.max # sentinel for padding

def calc_range(genome, track_index):
    # not limited to include_coords, so scale is always the same
    return genome.mins[track_index], genome.maxs[track_index]

def calc_histogram(genome, track_index, data_range, num_bins, include_coords):
    histogram_custom = partial(histogram, bins=num_bins, range=data_range,
                               new=True)

    hist, edges = histogram_custom(array([]))

    for chromosome in genome:
        for supercontig, continuous in chromosome.itercontinuous():
            continuous_track = continuous[:, track_index]
            supercontig_hist, new_edges = histogram_custom(continuous_track)

            assert edges == new_edges

            hist += supercontig_hist

    return hist, edges

def print_histogram(hist, edges):
    for row in zip(edges, hist.tolist() + ["NA"]):
        print "\t".join(map(str, row))

def _histogram(genomedataname, trackname, num_bins, include_coords_filename=None,
               include_identify_filelistname=None, identify_label=1):
    print "\t".join(FIELDNAMES) # lower_edge, count

    with Genome(genomedataname) as genome:
        track_index = genome.index_continuous(trackname)

        # go through data in two passes to avoid running out of memory
        # pass 1: calculate just the range
        data_range = calc_range(genome, track_index)

        # pass 2: calculate the histogram
        hist, edges = calc_histogram(genome, track_index, data_range, num_bins,
                                     include_coords)

    print_histogram(hist, edges)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... ARCHIVE TRACKNAME"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("--include-coords", metavar="FILE",
                      help="limit summary to genomic coordinates in FILE")

    parser.add_option("-b", "--num-bins", metavar="BINS", type=int,
                      default=100, help="use BINS bins")

    options, args = parser.parse_args(args)

    if not len(args) >= 1:
        parser.print_usage()
        sys.exit(1)

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    return _histogram(*args, options.num_bins,
                      options.include_coords, options.include_identify,
                      options.identify_label)

if __name__ == "__main__":
    sys.exit(main())

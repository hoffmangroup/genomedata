#!/usr/bin/env python
from __future__ import division, with_statement

"""
_histogram: prints histogram
"""

__version__ = "$Revision$"

# Copyright 2008-2014 Michael M. Hoffman <mhoffman@uhnresearch.ca>

import sys

from functools import partial
from numpy import array, histogram

from . import Genome

FIELDNAMES = ["lower_edge", "count"]


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


def _histogram(genomedataname, trackname, num_bins, include_coords):
    print "\t".join(FIELDNAMES)  # lower_edge, count

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

    if options.include_coords:
        raise NotImplementedError

    return _histogram(*args, num_bins=options.num_bins,
                      include_coords=options.include_coords)

if __name__ == "__main__":
    sys.exit(main())

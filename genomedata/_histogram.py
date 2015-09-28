#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
_histogram: prints histogram
"""

# Copyright 2008-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

import sys

from argparse import ArgumentParser
from functools import partial

from numpy import array, histogram

from . import Genome, __version__

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
        print(*row, sep="\t")


def _histogram(genomedataname, trackname, num_bins, include_coords):
    print(FIELDNAMES, sep="\t")  # lower_edge, count

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

    description = ("Print a histogram of values from a genomedata"
                   " archive")

    parser = ArgumentParser(description=description,
                            prog='genomedata-histogram',
                            version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')
    parser.add_argument('trackname', help='track name')

    parser.add_argument("--include-coords", metavar="FILE",
                        help="limit summary to genomic coordinates in FILE")

    parser.add_argument("-b", "--num-bins", metavar="BINS", type=int,
                        default=100, help="use BINS bins")

    args = parser.parse_args(args)

    return args


def main(argv=sys.argv[1:]):
    args = parse_options(argv)

    if args.include_coords:
        raise NotImplementedError

    return _histogram(args.gdarchive, args.trackname, args.num_bins,
                      args.include_coords)

if __name__ == "__main__":
    sys.exit(main())

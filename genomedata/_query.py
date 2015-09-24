#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
_query: DESCRIPTION
"""

# Copyright 2012, 2013 Michael M. Hoffman <mmh1@uw.edu>

import sys

from argparse import ArgumentParser

from . import Genome, __version__

def _query(filename, trackname, chromosome_name, begin, end):

    with Genome(filename) as genome:
        chromosome = genome[chromosome_name]
        track_index = genome.index_continuous(trackname)
        data = chromosome[begin:end, track_index]
        for index, point in enumerate(data):
            print(point)


def parse_options(args):

    description = ('print data from genomedata archive in specified '
                   ' trackname and coordinates')

    parser = ArgumentParser(prog='genomedata-query',
                            description=description,
                            version=__version__)
    
    parser.add_argument('gdarchive', help='genomedata archive')
    parser.add_argument('trackname', help='track name')
    parser.add_argument('chrom', help='chromosome name')
    parser.add_argument('begin', help='chromosome start', type=int)
    parser.add_argument('end', help='chromosome end', type=int)

    args = parser.parse_args(args)

    return args


def main(argv=sys.argv[1:]):
    args = parse_options(argv)

    return _query(args.gdarchive, args.trackname, args.chrom,
                  args.begin, args.end)

if __name__ == "__main__":
    sys.exit(main())

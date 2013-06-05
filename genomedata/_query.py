#!/usr/bin/env python
from __future__ import division, with_statement

"""
_query: DESCRIPTION
"""

__version__ = "$Revision$"

# Copyright 2012, 2013 Michael M. Hoffman <mmh1@uw.edu>

import sys

from . import Genome


def die(msg="Unexpected error."):
    print >>sys.stderr, msg
    sys.exit(1)


def _query(filename, trackname, chromosome_name, begin, end):
    begin = int(begin)
    end = int(end)
    with Genome(filename) as genome:
        chromosome = genome[chromosome_name]
        track_index = genome.index_continuous(trackname)
        data = chromosome[begin:end, track_index]
        for index, point in enumerate(data):
            print point


def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... ARCHIVE TRACKNAME CHROM BEGIN END"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    options, args = parser.parse_args(args)

    if not len(args) == 5:
        parser.error("incorrect number of arguments")

    return options, args


def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    return _query(*args)

if __name__ == "__main__":
    sys.exit(main())

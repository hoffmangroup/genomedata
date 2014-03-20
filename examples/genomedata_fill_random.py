#!/bin/env python

"""
Fills all data tracks in the Genomedata archive with random numbers
between [0, 1).
"""

from __future__ import division, with_statement

__version__ = "$Revision: "

import sys

from genomedata import Genome
from numpy.random import rand

def genomedata_fill_random(gdfilename):
    with Genome(gdfilename, mode="r+") as genome:
        print >>sys.stderr, "Opening %s..." % gdfilename
        for chromosome in genome:
            print >>sys.stderr, "Overwriting %s with random data" % chromosome
            for supercontig, continuous in chromosome.itercontinuous():
                continuous[...] = rand(*continuous.shape)

def parse_args(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] GENOMEDATAFILE"
    parser = OptionParser(usage=usage, version=__version__,
                          description=__doc__.strip())

    options, args = parser.parse_args(args)

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    genomedata_fill_random(args[0])

if __name__ == "__main__":
    sys.exit(main())

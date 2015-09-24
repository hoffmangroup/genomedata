#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
report: report some summary statistics on a genomedata archive that
already has save_metadata() run
"""

# Copyright 2009-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

import sys

from argparse import ArgumentParser

from . import Genome, __version__
from tabdelim import ListWriter

def report(gdarchive):
    writer = ListWriter()

    with Genome(gdarchive) as genome:
        writer.writerow(["measurement"] + genome.tracknames_continuous)
        writer.writerow(["mean"] + list(genome.means))
        writer.writerow(["var"] + list(genome.vars))

def parse_options(args):

    parser = ArgumentParser(prog='genomedata-report',
                            version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')

    args = parser.parse_args(args)

    return args

def main(argv=sys.argv[1:]):
    args = parse_options(argv)

    return report(args.gdarchive)

if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python

from __future__ import absolute_import, division, print_function
from future_builtins import ascii, filter, hex, map, oct, zip

"""
report: report some summary statistics on a genomedata archive that
already has save_metadata() run
"""

# Copyright 2009-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

import sys

from genomedata import Genome
from tabdelim import ListWriter

def report(gdarchive):
    writer = ListWriter()

    with Genome(archive) as genome:
        writer.writerow(["measurement"] + genome.tracknames_continuous)
        writer.writerow(["mean"] + list(genome.means))
        writer.writerow(["var"] + list(genome.vars))

def parse_options(args):

    from argparse import ArgumentParser
    from . import __version__

    parser = ArgumentParser(prog='genomedata-report',
                            version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')

    args = parser.parse_args(args)

    return args

def main(args=sys.argv[1:]):
    args = parse_options(args)

    return report(args.gdarchive)

if __name__ == "__main__":
    sys.exit(main())

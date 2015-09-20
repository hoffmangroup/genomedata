#!/usr/bin/env python

from __future__ import absolute_import, division, print_function
from future_builtins import ascii, filter, hex, map, oct, zip

"""
Open one or more tracks in the specified Genomedata archive.
These tracks can then be loaded with data using genomedata-load-data.
"""

# Copyright 2008-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

import sys

import warnings

from . import Genome

def open_data(gdarchive, tracknames, verbose):
    warnings.simplefilter("ignore")
    with Genome(gdarchive, "r+") as genome:
        # XXXopt: it would be more efficient to add them all at once
        for trackname in tracknames:
            genome.add_track_continuous(trackname)

    warnings.resetwarnings()

def parse_options(args):

    from argparse import ArgumentParser
    from . import __version__

    description = ("Open one or more tracks in"
                   " the specified Genomedata archive.")

    parser = ArgumentParser(description=description,
                            prog='genomedata-open-data',
                            version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')

    parser.add_argument("--trackname", required=True, nargs='+', 
                        help="tracknames to open")

    parser.add_argument("--verbose", default=False, action="store_true",
                        help="Print status updates and diagnostic messages")

    args = parser.parse_args(args)

    return args

def main(args=sys.argv[1:]):
    args = parse_options(args)
    kwargs = {"verbose": options.verbose}
    return open_data(args.gdarchive, args.tracknames, **kwargs)

if __name__ == "__main__":
    sys.exit(main())

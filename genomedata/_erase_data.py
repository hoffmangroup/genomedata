#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
_erase_data.py: wipe the data for a specified track in such a way
that the data can be reloaded (essentially replacing the track data)
"""

import sys

from argparse import ArgumentParser

from . import Genome, __version__

LINE_WIDTH = 70

def erase_data(gdfilename, trackname, verbose=False):
    if verbose:
        print("Erasing data for track: %s." % trackname, file=sys.stderr)
        print("Each dot is a chromosome erased:", file=sys.stderr)

    with Genome(gdfilename, mode="r+") as genome:
        dot_count = 0
        for chromosome in genome:
            chromosome._erase_data(trackname)
            if verbose:
                sys.stderr.write(".")
                dot_count += 1

                # Make sure this doesn't get rediculous if there are
                # thousands of chromosomes
                if dot_count % LINE_WIDTH == 0:
                    sys.stderr.write("\n")
        if verbose:
            print("\ndone", file=sys.stderr)

def parse_options(args):

    description = ("Erase the specified tracks from the Genomedata archive"
                   " in such a way that the track data can be replaced"
                   " (via genomedata-load-data).")

    parser = ArgumentParser(description=description,
                            prog='genomedata-erase-data',
                            version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')

    parser.add_argument("--trackname", required=True, nargs='+', 
                        help="tracknames to erase")

    parser.add_argument("--verbose", default=False, action="store_true",
                        help="Print status updates and diagnostic messages")

    args = parser.parse_args(args)

    return args

def main(argv=sys.argv[1:]):
    args = parse_options(argv)
    gdarchive = args.archive

    for trackname in args.tracknames:
        erase_data(gdarchive, trackname, verbose=args.verbose)

if __name__ == "__main__":
    sys.exit(main())

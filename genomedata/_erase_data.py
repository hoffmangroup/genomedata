#!/usr/bin/env python
from __future__ import division, with_statement

"""
_erase_data.py: wipe the data for a specified track in such a way
that the data can be reloaded (essentially replacing the track data)
"""

__version__ = "$Revision$"

import sys

from . import Genome

LINE_WIDTH = 70

def erase_data(gdfilename, trackname, verbose=False):
    if verbose:
        print >>sys.stderr, "Erasing data for track: %s." % trackname
        print >>sys.stderr, "Each dot is a chromosome erased:"

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
            print >>sys.stderr, "\ndone"

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... GENOMEDATAFILE TRACKNAME..."
    version = "%%prog %s" % __version__
    description = ("Erase the specified tracks from the Genomedata archive"
                   " in such a way that the track data can be replaced"
                   " (via genomedata-load-data).")
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-v", "--verbose", dest="verbose",
                      default=False, action="store_true",
                      help="Print status updates and diagnostic messages")

    options, args = parser.parse_args(args)

    if not len(args) >= 2:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    gdfilename = args[0]
    tracknames = args[1:]
    kwargs = {"verbose": options.verbose}

    for trackname in tracknames:
        erase_data(gdfilename, trackname, **kwargs)

if __name__ == "__main__":
    sys.exit(main())

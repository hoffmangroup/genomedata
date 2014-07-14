#!/usr/bin/env python

"""
Open one or more tracks in the specified Genomedata archive.
These tracks can then be loaded with data using genomedata-load-data.
"""

from __future__ import division, with_statement

__version__ = "$Revision$"

# Copyright 2008-2014 Michael M. Hoffman <mhoffman@uhnresearch.ca>

import sys

import warnings

from . import Genome

def open_data(gdfilename, tracknames, verbose=False):
    warnings.simplefilter("ignore")
    with Genome(gdfilename, "r+") as genome:
        # XXXopt: it would be more efficient to add them all at once
        for trackname in tracknames:
            genome.add_track_continuous(trackname)

    warnings.resetwarnings()

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... GENOMEDATAFILE TRACKNAME..."
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=__doc__.strip())

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
    return open_data(gdfilename, tracknames, **kwargs)

if __name__ == "__main__":
    sys.exit(main())

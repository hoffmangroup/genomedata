#!/usr/bin/env python
from __future__ import division, with_statement

"""
_open_data: recursively set up HDF5 files to have named tracks
"""

__version__ = "$Revision$"

# Copyright 2008-2009 Michael M. Hoffman <mmh1@washington.edu>

import sys

from numpy import array
import warnings

from . import Genome

def open_data(gdfilename, tracknames, verbose=False):
    warnings.simplefilter("ignore")
    with Genome(gdfilename, "r+") as genome:
        for chromosome in genome:
            attrs = chromosome.attrs
            if "tracknames" in attrs:
                raise ValueError("%s already has named tracks" % gdfilename)

            attrs.dirty = True
            attrs.tracknames = array(tracknames)

    warnings.resetwarnings()

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... GENOMEDATAFILE TRACKNAME..."
    version = "%%prog %s" % __version__
    description = ("Specify the tracks that will be in this Genomedata"
                   " archive")
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
    return open_data(gdfilename, tracknames, **kwargs)

if __name__ == "__main__":
    sys.exit(main())

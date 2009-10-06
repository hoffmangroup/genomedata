#!/usr/bin/env python
from __future__ import division, with_statement

"""
_open_data: recursively set up HDF5 files to have named tracks
"""

__version__ = "$Revision$"

# Copyright 2008-2009 Michael M. Hoffman <mmh1@washington.edu>

import sys

from numpy import array
from path import path
from tables import openFile

def open_data(dirname, tracknames):
    dirpath = path(dirname)
    for filepath in dirpath.walkfiles():
        with openFile(filepath, "r+") as h5file:
            attrs = h5file.root._v_attrs

            if "tracknames" in attrs:
                raise ValueError("%s already has named tracks" % filepath)

            attrs.dirty = True
            attrs.tracknames = array(tracknames)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... GENOMEDATADIR TRACKNAME..."
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    options, args = parser.parse_args(args)

    if not len(args) >= 2:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    genomedatadir = args[0]
    tracknames = args[1:]
    return open_data(genomedatadir, tracknames)

if __name__ == "__main__":
    sys.exit(main())

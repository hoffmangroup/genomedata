#!/usr/bin/env python
from __future__ import division, with_statement

"""
_erase_track.py: wipe the data for a specified track in such a way
that the data can be reloaded (essentially replacing the track data)
"""

__version__ = "$Revision: $"

import sys

from numpy import NAN
from path import path
from tables import NoSuchNodeError, openFile

from ._util import fill_array, get_tracknames, walk_supercontigs

def clear_track_data(chromosome, trackname, verbose=False):
    # XXX: should replace with Chromosome.index_continuous
    tracknames = get_tracknames(chromosome)

    if trackname in tracknames:
        col_index = tracknames.index(trackname)
    else:
        raise ValueError("Could not find track: %s in chromosome: %s" % \
                             (trackname, chromosome))

    chromosome.root._v_attrs.dirty = True

    # XXX: should replace with iter(Chromosome)
    for supercontig in walk_supercontigs(chromosome):
        try:
            new_data = fill_array(NAN, (supercontig.continuous.shape[0], ),
                                  dtype = supercontig.continuous.atom.dtype)

            supercontig.continuous[:, col_index] = new_data
        except NoSuchNodeError:
            print >>sys.stderr, "  No continuous node found for %s in %s" % \
                (supercontig._v_name, chromosome.title)
            continue

def erase_track(dirname, trackname, verbose=False):
    if verbose:
        print >>sys.stderr, "Erasing data for track: %s" % trackname

    dirpath = path(dirname)
    for filepath in dirpath.walkfiles():
        with openFile(filepath, "r+") as h5file:
            clear_track_data(h5file, trackname, verbose=verbose)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... GENOMEDATADIR TRACKNAME..."
    version = "%%prog %s" % __version__
    description = ("Erase the specified tracks from the Genomedata collection"
                   " in such a way that the track can be replaced.")
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
    genomedatadir = args[0]
    tracknames = args[1:]
    kwargs = {"verbose": options.verbose}

    for trackname in tracknames:
        erase_track(genomedatadir, trackname, **kwargs)

if __name__ == "__main__":
    sys.exit(main())

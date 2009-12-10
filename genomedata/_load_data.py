#!/usr/bin/env python
from __future__ import division, with_statement

"""
_load_data.py: A python interface for genome_load_data.c
"""

__version__ = "$Revision: $"

import sys

from subprocess import PIPE, Popen

from ._util import EXT_GZ

DEFAULT_CHUNK_SIZE = 10000
LOAD_DATA_CMD = "genomedata-load-data"

def die(msg="Unexpected error!"):
    print >>sys.stderr, msg
    sys.exit(1)

def load_data(genomedatadir, trackname, datafile,
              chunk_size=None, verbose=False):
    """Loads data from datafile into specific track of Genomedata collection

    genomedatadir: genomedata collection path
    trackname: name of track (as specified in open_data) to load data for
    datafile: file to read data from
    chunk_size: number of rows in each hdf5 data chunk

    """
    if verbose:
        print ">> Loading data for track: %s" % trackname

    if datafile.endswith(EXT_GZ):
        read_cmd = ["zcat"]
    else:
        read_cmd = ["cat"]
    read_cmd.append(datafile)

    load_cmd = [LOAD_DATA_CMD]
    if verbose:
        load_cmd.append("-v")
    if chunk_size:
        load_cmd.append("--chunk-size=%d" % chunk_size)
    load_cmd.extend([genomedatadir, trackname])

    # Pipe read command into load command
    reader = Popen(read_cmd, stdout=PIPE)
    loader = Popen(load_cmd, stdin=reader.stdout)
    loader.communicate()
    retcode = loader.poll()
    if retcode != 0:
        die("Error loading data from track file: %s" % datafile)

def parse_args(args):
    from optparse import OptionParser
    usage = "%prog [OPTION]... GENOMEDATADIR TRACKNAME DATAFILE"
    version = "%%prog %s" %__version__
    description = ("Load data from DATAFILE into the specified TRACKNAME"
                   " of the Genomedata collection")
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-c", "--chunk-size", dest="chunk_size",
                      metavar="NROWS", type="int",
                      default=DEFAULT_CHUNK_SIZE,
                      help="Chunk hdf5 data into blocks of NROWS."
                      " A higher value increases compression but slows"
                      " random access. Must always be smaller than the"
                      " max size for a dataset. [default: %default]")
    parser.add_option("-v", "--verbose", dest="verbose",
                      default=False, action="store_true",
                      help="Print status and diagnostic messages")

    options, args = parser.parse_args(args)

    if not len(args) == 3:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    kwargs = {"verbose": options.verbose,
              "chunk_size": options.chunk_size}
    load_data(*args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())

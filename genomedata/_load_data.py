#!/usr/bin/env python
from __future__ import division, with_statement

"""
_load_data.py: A python interface for genome_load_data.c
"""

__version__ = "$Revision$"

import sys

from subprocess import PIPE, Popen

from ._util import SUFFIX_GZ

DEFAULT_CHUNK_SIZE = 10000
LOAD_DATA_CMD = "genomedata-load-data"
MSG_LOAD_ERROR = "Error loading data from track file %%s. %s returned %%s." % LOAD_DATA_CMD

def die(msg="Unexpected error."):
    print >>sys.stderr, msg
    sys.exit(1)

def load_data(gdfilename, trackname, datafile, verbose=False):
    """Loads data from datafile into specific track of Genomedata archive

    gdfilename: genomedata archive path
    trackname: name of track (as specified in open_data) to load data for
    datafile: file to read data from
    """
    if verbose:
        print ">> Loading data for track: %s" % trackname

    if datafile.endswith(SUFFIX_GZ):
        read_cmd = ["zcat"]
    else:
        read_cmd = ["cat"] # XXX: useless use of cat
    read_cmd.append(datafile)

    load_cmd = [LOAD_DATA_CMD]
    if verbose:
        load_cmd.append("-v")
    load_cmd.extend([gdfilename, trackname])

    if verbose:
        read_cmdline = " ".join(read_cmd)
        load_cmdline = " ".join(load_cmd)
        print >>sys.stderr, "%s | %s" % (read_cmdline, load_cmdline)

    # Pipe read command into load command
    reader = Popen(read_cmd, stdout=PIPE)
    loader = Popen(load_cmd, stdin=reader.stdout)
    loader.communicate()
    retcode = loader.poll()
    if retcode != 0:
        die(MSG_LOAD_ERROR % (datafile, retcode))

def parse_args(args):
    from optparse import OptionParser
    usage = "%prog [OPTION]... GENOMEDATAFILE TRACKNAME DATAFILE"
    version = "%%prog %s" %__version__
    description = ("Load data from DATAFILE into the specified TRACKNAME"
                   " of the Genomedata archive")
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

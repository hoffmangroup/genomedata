#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
_load_data.py: A python interface for genome_load_data.c
"""

__version__ = "$Revision$"

import sys

from subprocess import PIPE, Popen

from ._util import SUFFIX_GZ, die

DEFAULT_CHUNK_SIZE = 10000
LOAD_DATA_CMD = "genomedata-load-data"
MSG_LOAD_ERROR = "Error loading data from track file %%s. %s returned %%s." % LOAD_DATA_CMD

def load_data(gdfilename, trackname, datafile, verbose=False):
    """Loads data from datafile into specific track of Genomedata archive

    gdfilename: genomedata archive path
    trackname: name of track (as specified in open_data) to load data for
    datafile: file to read data from
    """
    if verbose:
        print(">> Loading data for track: %s" % trackname)

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
        print("%s | %s" % (read_cmdline, load_cmdline), file=sys.stderr)

    # Pipe read command into load command
    reader = Popen(read_cmd, stdout=PIPE)
    loader = Popen(load_cmd, stdin=reader.stdout)
    loader.communicate()
    retcode = loader.poll()
    if retcode != 0:
        die(MSG_LOAD_ERROR % (datafile, retcode))

def parse_args(args):

    from argparse import ArgumentParser
    from . import __version__

    description = ("Load data from DATAFILE into the specified TRACKNAME"
                   " of the Genomedata archive")

    parser = ArgumentParser(description=description,
                            prog='genomedata-load-data',
                            version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')
    parser.add_argument('trackname', help='track name')
    parser.add_argument('datafile', help='data file')

    parser.add_argument("-c", "--chunk-size",
                        metavar="NROWS", type=int,
                        default=DEFAULT_CHUNK_SIZE,
                        help="Chunk hdf5 data into blocks of NROWS."
                        " A higher value increases compression but slows"
                        " random access. Must always be smaller than the"
                        " max size for a dataset. [default: %(default)s]")

    parser.add_argument("--verbose", default=False, action="store_true",
                      help="Print status and diagnostic messages")

    args = parser.parse_args(args)

    return args


def main(args=sys.argv[1:]):
    args = parse_args(args)
    kwargs = {"verbose": args.verbose,
              "chunk_size": args.chunk_size}
    load_data(args.gdarchive, args.trackname, args.datafile, **kwargs)

if __name__ == "__main__":
    sys.exit(main())

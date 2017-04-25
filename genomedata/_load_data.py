#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
_load_data.py: A python interface for genome_load_data.c
"""

import sys

from argparse import ArgumentParser
from pybedtools import cleanup
from subprocess import PIPE, Popen

from . import __version__
from ._util import SUFFIX_GZ, die, get_bed_from_track_file, score_filter

BIG_WIG_READ_CMD = "bigWigToBedGraph"
DEFAULT_CHUNK_SIZE = 10000
LOAD_DATA_CMD = "genomedata-load-data"
MSG_LOAD_ERROR = "Error loading data from track file %%s. %s returned %%s." % LOAD_DATA_CMD


def load_data(gdfilename, trackname, datafile, filterfile=None,
              filter_threshold=None, verbose=False):
    """Loads data from datafile into specific track of Genomedata archive

    gdfilename: genomedata archive path
    trackname: name of track (as specified in open_data) to load data for
    datafile: file to read data from
    """
    if verbose:
        print(">> Loading data for track: %s" % trackname)

    # Check wigFix
    # Get BED from genomic file
    input_bed = get_bed_from_track_file(datafile)
    # If there is a filter being appled
    if filterfile:
        if verbose:
            print(">> Filtering ", trackname, " with ", filterfile)
        # Get BED from filter input file
        filter_bed = get_bed_from_track_file(filterfile)

        # If there is a threshold given
        if filter_threshold:
            # Filter the filter by the threshold given
            filter_bed = filter_bed.filter(score_filter,
                                           threshold=filter_threshold)

        # Interesect the filter with the input data track
        result_bed = input_bed.intersect(filter_bed)
        # Merge results
        result_bed = result_bed.merge()
    # Otherwise just use the input genomic file
    else:
        result_bed = input_bed

    if verbose:
        print(">> Using processed track data from", result_bed.fn)

    # file_is_big_wig = is_big_wig(datafile)

    # if file_is_big_wig:
    #     read_cmd = [BIG_WIG_READ_CMD]
    # elif datafile.endswith(SUFFIX_GZ):
    #     read_cmd = ["zcat"]
    # else:
    #     read_cmd = ["cat"]  # XXX: useless use of cat
    # read_cmd.append(datafile)
    #
    # # If the file is a BigWig
    # if file_is_big_wig:
    #     # set bigWigToBedGraph to output to stdout
    #     read_cmd.append("/dev/stdout")

    load_cmd = [LOAD_DATA_CMD]
    if verbose:
        load_cmd.append("-v")
    load_cmd.extend([gdfilename, trackname])

    if verbose:
        # read_cmdline = " ".join(read_cmd)
        load_cmdline = " ".join(load_cmd)
        print(load_cmdline, result_bed.fn, sep=" < ",  file=sys.stderr)

    with open(result_bed.fn) as input_track_bed_file:
        loader = Popen(load_cmd, stdin=input_track_bed_file)
        retcode = loader.wait()

    if retcode != 0:
        die(MSG_LOAD_ERROR % (datafile, retcode))

    # Force clean up of temporary files only made in this function
    cleanup()


def parse_args(args):

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


def main(argv=sys.argv[1:]):
    args = parse_args(argv)
    load_data(args.gdarchive, args.trackname, args.datafile,
              verbose=args.verbose, chunk_size=args.chunk_size)

if __name__ == "__main__":
    sys.exit(main())

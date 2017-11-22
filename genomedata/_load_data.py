#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
_load_data.py: A python interface for genome_load_data.c
"""

import struct
import sys

from argparse import ArgumentParser
from subprocess import PIPE, Popen

from . import __version__
from ._util import SUFFIX_GZ, die

BIG_WIG_SIGNATURE = 0x888FFC26
BIG_WIG_SIGNATURE_BYTE_SIZE = 4
BIG_WIG_READ_CMD = "bigWigToBedGraph"
BEDTOOLS_CMD = "bedtools"
BEDTOOLS_INTERSECT_CMD = "intersect"
BEDTOOLS_STDIN = "stdin"
DEFAULT_CHUNK_SIZE = 10000
LOAD_DATA_CMD = "genomedata-load-data"
MSG_LOAD_ERROR = "Error loading data from track file %%s. %s returned %%s." % LOAD_DATA_CMD


def load_data(gdfilename, trackname, datafile, maskfile=None, verbose=False):
    """Loads data from datafile into specific track of Genomedata archive

    gdfilename: genomedata archive path
    trackname: name of track (as specified in open_data) to load data for
    datafile: file to read data from
    maskfile: BED file to mask out regions from datafile
    """
    if verbose:
        print(">> Loading data for track: %s" % trackname)

    file_is_big_wig = is_big_wig(datafile)

    if file_is_big_wig:
        read_cmd = [BIG_WIG_READ_CMD]
    elif datafile.endswith(SUFFIX_GZ):
        read_cmd = ["zcat"]
    else:
        read_cmd = ["cat"]  # XXX: useless use of cat
    read_cmd.append(datafile)

    # If the file is a BigWig
    if file_is_big_wig:
        # set bigWigToBedGraph to output to stdout
        read_cmd.append("/dev/stdout")

    # If a mask file was specified
    if maskfile:
        # Construct the mask cmd (using bedtools)
        # e.g. bigWigToBedGraph HISTONE.bigWig /dev/stdout | bedtools
        # intersect -a Umap.bedGraph.gz -b stdin | genomedata-load-data
        # -a and -b specify the file names to bedtools

        # TODO: make masking only applicable on bedGraph files
        # and ensure the mask file is also a bed
        mask_cmd = [BEDTOOLS_CMD, BEDTOOLS_INTERSECT_CMD,
                    "-a", BEDTOOLS_STDIN, "-b", maskfile]
    # Otherwise do not have a mask command
    else:
        mask_cmd = []

    load_cmd = [LOAD_DATA_CMD]
    if verbose:
        load_cmd.append("-v")
    load_cmd.extend([gdfilename, trackname])

    if verbose:
        read_cmdline = " ".join(read_cmd)
        mask_cmdline = " ".join(mask_cmd)
        load_cmdline = " ".join(load_cmd)
        print(" | ".join(filter(None, [read_cmdline,
                                       mask_cmdline,
                                       load_cmdline])),
              file=sys.stderr)

    # Open the read command
    try:
        reader = Popen(read_cmd, stdout=PIPE)
    except OSError as os_exception:
        # If it was a big wig file and the converting program was not found on
        # the path
        if file_is_big_wig:
            # Report that the program could not be found
            raise OSError("A BigWig file was detected but " +
                          BIG_WIG_READ_CMD +
                          " could not be found on your path to convert it for "
                          "use in Genomedata. Please install the utility "
                          "from UCSC and try again.")
        # Otherwise if there was another exception where the file isn't BigWig
        else:
            # Re-raise the exception
            raise os_exception

    # If a maskfile was specified
    if maskfile:
        # Pipe the read command into the mask command
        try:
            masker = Popen(mask_cmd, stdin=reader.stdout, stdout=PIPE)
            loader_input_process = masker
        except OSError:
            # If we could not find the bedtools command
            raise OSError("The bedtools command is necessary for masking out "
                          "the data file and could not be found. Please "
                          "install bedtools and ensure it is on your PATH.")
    # Otherwise have the read command pipe straight to the load command
    else:
        loader_input_process = reader

    loader = Popen(load_cmd, stdin=loader_input_process.stdout)
    loader.communicate()
    retcode = loader.poll()
    if retcode != 0:
        die(MSG_LOAD_ERROR % (datafile, retcode))


def is_big_wig(filename):
    """ Checks that the given filename refers to a valid bigWig file """
    with open(filename, "rb") as big_wig_file:
        signature_string = big_wig_file.read(BIG_WIG_SIGNATURE_BYTE_SIZE)

    # unpack returns a tuple regardless of length
    # the kent reference checks both little endian and big endian packing
    # of the 4 byte signature
    little_endian_signature = struct.unpack("<L", signature_string)[0]
    big_endian_signature = struct.unpack(">L", signature_string)[0]

    if (little_endian_signature == BIG_WIG_SIGNATURE or
       big_endian_signature == BIG_WIG_SIGNATURE):
        return True

    return False


def parse_args(args):

    description = ("Load data from DATAFILE into the specified TRACKNAME"
                   " of the Genomedata archive")

    parser = ArgumentParser(description=description,
                            prog='genomedata-load-data',
                            version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')
    parser.add_argument('trackname', help='track name')
    parser.add_argument('datafile', help='data file')

    parser.add_argument("-m", "--maskfile",
                        help='A BED file containing regions to mask out from'
                        'the data file')

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
              args.maskfile, verbose=args.verbose, chunk_size=args.chunk_size)

if __name__ == "__main__":
    sys.exit(main())

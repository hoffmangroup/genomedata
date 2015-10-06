#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
load_genomedata: DESCRIPTION
"""

# Copyright 2009, 2011 Orion Buske <orion.buske@gmail.com>
# Copyright 2010 Michael Hoffman <mmh1@uw.edu>

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datetime import datetime
from glob import glob
from os import close, extsep
from path import path
from subprocess import call
import sys
from tempfile import mkdtemp, mkstemp

from . import EXT, FILE_MODE_CHROMS, SUFFIX, __version__
from ._load_seq import load_seq
from ._open_data import open_data
from ._load_data import DEFAULT_CHUNK_SIZE, load_data
from ._close_data import close_data
from ._util import die

def print_timestamp(msg=""):
    print(">> %s: %s" % (datetime.now().isoformat(), msg), file=sys.stderr)

def repack(infilename, outfilename, verbose=False):
    if verbose:
        print(">> Repacking: %s -> %s" % (infilename,
                                          outfilename), file=sys.stderr)

    retcode = call(["h5repack", "-f", "GZIP=1", infilename, outfilename])
    if retcode != 0:
        die("HDF5 repacking failed.")

def load_genomedata(gdfilename, tracks=None, seqfilenames=None, mode=None,
                    seqfile_type="fasta", chunk_size=DEFAULT_CHUNK_SIZE,
                    verbose=False):
    """Loads Genomedata collection with given data

    gdfilename: name of Genomedata archive to create
    tracks: a list of tracks to add, with a (track_name, track_filename) tuple
      for each track to add
    seqfilenames: list of filenames containing sequence data to add
    mode: "dir", "file", or None (decide based upon number of sequences)
    """

    gdpath = path(gdfilename).expand()
    try:
        if mode is None:
            if (seqfilenames is not None
                and len(seqfilenames) > FILE_MODE_CHROMS):
                mode = "file"
            else:
                mode = "dir"

        if mode == "dir":
            # Generate hdf5 data in temporary directory and copy out when done
            tempdatadir = mkdtemp(prefix=(EXT + extsep))
            tempdatapath = path(tempdatadir)
            isdir = True
        elif mode == "file":
            tempdatafile, tempdatafilename = mkstemp(suffix=SUFFIX)
            tempdatapath = path(tempdatafilename)
            # Close the file _descriptor_ (not a file object)
            close(tempdatafile)
            # Delete the file to allow load_seq to create it
            tempdatapath.remove()
            isdir = False
        else:
            raise ValueError("Unknown mode: %s" % mode)

        if verbose:
            print(">> Using temporary Genomedata archive: %s" % tempdatapath, file=sys.stderr)

        # Load sequences if any are specified
        if not seqfilenames:
            raise ValueError("No sequence files specified.")

        for seqfilename in seqfilenames:
            if seqfile_type == "fasta":
                seqfile_desc = "sequence"
            else:
                seqfile_desc = "assembly"

            if not path(seqfilename).isfile():
                die("Could not find %s file: %s" % (seqfile_desc, seqfilename))

        if verbose:
            print_timestamp("Loading %s files:" % seqfile_desc)

        load_seq(tempdatapath, seqfilenames, verbose=verbose, mode=mode,
                 seqfile_type=seqfile_type)

        # Load tracks if any are specified
        if tracks is not None and len(tracks) > 0:
            # Open hdf5 with track names
            try:
                track_names = []
                for track_name, track_filename in tracks:
                    if path(track_filename).isfile():
                        if track_name not in track_names:  # No duplicates
                            track_names.append(track_name)
                    else:
                        die("Could not find track file: %s" % track_filename)
            except ValueError:
                die("Error saving data from tracks: %s" % tracks)

            if verbose:
                print_timestamp("Opening Genomedata archive with %d tracks" %
                                len(track_names))

            open_data(tempdatapath, track_names, verbose=verbose)

            # Load track data
            if verbose:
                print_timestamp("Loading data")

            for track_name, track_filename in tracks:
                load_data(tempdatapath, track_name, track_filename,
                          verbose=verbose)

        # Close genomedata
        try:
            close_data(tempdatapath, verbose=verbose)
        except:
            die("Error saving metadata.")

        # Make output directory
        if verbose:
            print_timestamp("Creating Genomedata archive: %s" % gdfilename)

        # Move/repack h5 files to output directory
        if isdir:  # Repack each h5 file separately
            if gdpath.exists():
                assert gdpath.isdir()
            else:
                gdpath.makedirs()

            for tempfilepath in tempdatapath.files():
                tempbasepath = tempfilepath.name
                outfilepath = gdpath / tempbasepath

                repack(tempfilepath, outfilepath, verbose)
        else:
            repack(tempdatapath, gdpath, verbose)
    except:
        print("Error creating genomedata.", file=sys.stderr)
        raise
    finally:
        try:
            # Remove temp directory and all contents
            if verbose:
                print(">> Cleaning up...", end=' ', file=sys.stderr)

            sys.stdout.flush()
            if tempdatapath.isfile():
                tempdatapath.remove()
            else:
                tempdatapath.rmtree()

            if verbose:
                print("done", file=sys.stderr)
        except Exception as e:
            print("\nCleanup failed: %s" % str(e), file=sys.stderr)

    if verbose:
        print("\n===== Genomedata archive successfully created: %s =====\n" % \
              gdfilename, file=sys.stderr)

    return gdfilename

def parse_cmdline(cmdline):

    description = ("Create Genomedata archive named GENOMEDATAFILE by loading\n"
                   " specified track data and sequences. If GENOMEDATAFILE\n"
                   " already exists, it will be overwritten.\n"
                   " --track and --sequence may be repeated to specify\n"
                   " multiple trackname=trackfile pairings and sequence files,\n"
                   " respectively.\n\n"
                   " Example: %(prog)s -t high=signal.high.wig -t"
                   " low=signal.low.bed.gz"
                   " -s chrX.fa -s chrY.fa.gz"
                   " GENOMEDATAFILE")

    citation = ("Citation: Hoffman MM, Buske OJ, Noble WS.\n"
                "2010. The Genomedata format for storing large-scale functional genomics data.\n"
                "Bioinformatics 26 (11):1458-1459.\n"
                "http://dx.doi.org/10.1093/bioinformatics/btq164")

    parser = ArgumentParser(description=description,
                            epilog=citation,
                            formatter_class=RawDescriptionHelpFormatter,
                            prog='genomedata-load',
                            version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive',
                        metavar='GENOMEDATAFILE')

    flags = parser.add_argument_group("Flags")
    flags.add_argument("--verbose",
                      default=False, action="store_true",
                      help="Print status updates and diagnostic messages")

    input_data = parser.add_argument_group("Input data")
    input_data.add_argument("-s", "--sequence", action='append', required=True,
                            default=None,
                            help="Add the sequence data in the specified file or files"
                            " (may use UNIX glob wildcard syntax)")
    input_data.add_argument("-t", "--track", action='append', 
                            default=None,
                            metavar="NAME=FILE", required=True,
                            help="Add data from FILE as the track NAME,"
                            " such as: -t signal=signal.wig")

    input_data_ex = input_data.add_mutually_exclusive_group()
    input_data_ex.add_argument("--assembly", action="store_const",
                            default=None,
                            const="agp", dest="seqfile_type",
                            help="sequence files contain assembly (AGP) files instead of"
                            " sequence")
    input_data_ex.add_argument("--sizes", action="store_const",
                            default=None,
                            const="sizes", dest="seqfile_type", 
                            help="sequence files contain list of sizes instead of"
                            " sequence")
    implementation = parser.add_argument_group("Implementation")
    implementation_ex = implementation.add_mutually_exclusive_group()
    implementation_ex.add_argument("-f", "--file-mode", dest="mode",
                                default='file',
                                action="store_const", const="file",
                                help="If specified, the Genomedata archive will be"
                                " implemented as a single file, with a separate h5 group"
                                " for each Chromosome. This is recommended if there are"
                                " a large number of Chromosomes. The default behavior is"
                                " to use a single file if there are at least %s"
                                " Chromosomes being added." % FILE_MODE_CHROMS)

    implementation_ex.add_argument("-d", "--directory-mode", dest="mode",
                                action="store_const", const="dir",
                                help="If specified, the Genomedata archive will be"
                                " implemented as a directory, with a separate file for"
                                " each Chromosome. This is recommended if there are a"
                                " small number of Chromosomes. The default behavior is"
                                " to use a directory if there are fewer than %s"
                                " Chromosomes being added." % FILE_MODE_CHROMS)

    args = parser.parse_args(cmdline)

    return args

def main(cmdline=sys.argv[1:]):
    args = parse_cmdline(cmdline)

    # default
    if not args.seqfile_type:
        args.seqfile_type = 'fasta'

    # list of lists
    seqfilenames_list = [glob(globname) for globname in args.sequence]
    seqfilenames = sum(seqfilenames_list, [])

    # Parse tracks into list of tuples
    try:
        tracks = []
        for track_expr in args.track:
            track_name, _, track_filename = track_expr.partition("=")
            tracks.append((track_name, track_filename))  # Tuple
    except ValueError:
        die(("Error parsing track expression: %s\Specify tracks"
             "in NAME=FILE form, such as: -t high=signal.high") % track_expr)

    load_genomedata(args.gdarchive, tracks, seqfilenames,
                    seqfile_type=args.seqfile_type, verbose=args.verbose,
                    mode=args.mode)

if __name__ == "__main__":
    sys.exit(main())

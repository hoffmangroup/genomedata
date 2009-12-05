#!/usr/bin/env python
from __future__ import division, with_statement

"""
load_genomedata: DESCRIPTION
"""

__version__ = "$Revision$"

# Copyright 2009 Orion Buske <orion.buske@gmail.com>

import os
import sys

from errno import EEXIST
from shutil import rmtree
from subprocess import call
from tempfile import mkdtemp, NamedTemporaryFile

from ._load_seq import load_seq
from ._open_data import open_data
from ._close_data import close_data
from ._util import EXT_GZ

LOAD_DATA_CMD = "genomedata-load-data"


def die(msg="Unexpected error!"):
    print >>sys.stderr, msg
    sys.exit(1)

def load_genomedata(genomedatadir, tracks=None, seqfiles=None,
                    chunk_size=None):
    """Loads genomedata object with given data

    genomedatadir: name of output directory for genomedata
    tracks: a list of tracks to add, with a (track_name, track_filename) tuple
      for each track to add
    seqfiles: list of filenames containing sequence data to add
    """

    try:
        # Generate hdf5 data in temporary directory and copy out when done
        tempdatadir = mkdtemp(prefix="genomedata.")
        print ">> Using temporary directory: %s" % tempdatadir

        # Load sequences if any are specified
        if seqfiles is not None and len(seqfiles) > 0:
            for seqfile in seqfiles:
                if not os.path.isfile(seqfile):
                    die("Could not find sequence file: %s" % seqfile)

            print ">> Loading sequence files:"
            load_seq(tempdatadir, seqfiles)

        # Load tracks if any are specified
        if tracks is not None and len(tracks) > 0:
            # Open hdf5 with track names
            try:
                track_names = []
                for track_name, track_filename in tracks:
                    if os.path.isfile(track_filename):
                        if track_name not in track_names:  # No duplicates
                            track_names.append(track_name)
                    else:
                        die("Could not find track file: %s" % track_filename)
            except ValueError:
                die("Error saving data from tracks: %s" % tracks)

            print ">> Opening genomedata with %d tracks" % len(track_names)
            open_data(tempdatadir, track_names)

            # Load track data
            for track_name, track_filename in tracks:
                print ">> Loading data for track: %s" % track_name
                cmd_args = [track_filename, "|", LOAD_DATA_CMD]
                if chunk_size:
                    cmd_args.append("--chunk-size=%d" % chunk_size)
                cmd_args.extend([tempdatadir, track_name])
                if track_filename.endswith(EXT_GZ):
                    cmd_args.insert(0, "zcat")
                else:
                    cmd_args.insert(0, "cat")

                # Need call in shell for piping
                retcode = call(" ".join(cmd_args), shell=True)
                if retcode != 0:
                    die("Error loading data from track file: %s" % \
                            track_filename)


        # Close hdf5 files
        try:
            close_data(tempdatadir)  # Close genomedata
        except:
            die("Error saving metadata!")

        # Make output directory
        if not os.path.isdir(genomedatadir):
            print ">> Creating directory: %s" % genomedatadir
            os.makedirs(genomedatadir)

        # Close and repack hdf5 files to output directory
        tempfiles = os.listdir(tempdatadir)
        for tempfilename in tempfiles:
            # Repack file to output dir
            print ">> Repacking: %s into genomedata" % tempfilename
            tempfilepath = os.path.join(tempdatadir, tempfilename)
            outfilepath = os.path.join(genomedatadir, tempfilename)
            cmd_args = ["h5repack", "-f", "GZIP=1", tempfilepath, outfilepath]
            retcode = call(cmd_args)
            if retcode != 0:
                die("HDF5 repacking failed!")
    except:
        print >>sys.stderr, "Error creating genomedata!"
        raise
    finally:
        try:
            # Remove temp directory and all contents
            print ">> Cleaning up...",
            sys.stdout.flush()
            rmtree(tempdatadir)
            print "done"
        except Exception, e:
            print >>sys.stderr, "\nCleanup failed: %s" % str(e)

    print "\n===== Genomedata successfully created in: %s =====\n" % \
        genomedatadir


def parse_options(args):
    from optparse import OptionParser

    usage = ("%prog [OPTIONS] GENOMEDATADIR"
             "\ne.g. %prog -t high=signal.high.wig -t low=signal.low.bed.gz"
             " -s chrX.fa -s chrY.fa.gz mygenomedata")
    version = "%%prog %s" % __version__
    description = ("--track and --sequence may be repeated to specify multiple"
                   " trackname=trackfile pairings and sequence files,"
                   " respectively")
    parser = OptionParser(usage=usage, version=version,
                          description=description)

#     parser.add_option("-c", "--chunk-size", dest="nrows",
#                       type="int", default=10000,
#                       help="Chunk hdf5 data into blocks of NROWS. A higher"
#                       " value increases compression but slows random access."
#                       " Must always be smaller than the max size for a"
#                       " dataset. [default: %default]")
    parser.add_option("-s", "--sequence", action="append",
                      dest="seqfile", default=[],
                      help="Add the sequence data in the specified file")
    parser.add_option("-t", "--track", action="append",
                      dest="track", default=[],
                      help="Add data for the given track. TRACK"
                      " should be specified in the form: NAME=FILE,"
                      " such as: -t signal=signal.wig")

    options, args = parser.parse_args(args)

    if not len(args) == 1:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    genomedatadir = args[0]
    seqfiles = options.seqfile

    # Parse tracks into list of tuples
    try:
        tracks = []
        for track_expr in options.track:
            track_name, track_filename = track_expr.split("=")
            tracks.append((track_name, track_filename))  # Tuple
    except ValueError:
        die(("Error parsing track expression: %s\nMake sure to specify tracks"
             "in NAME=FILE form, such as: -t high=signal.high") % track_expr)

    load_genomedata(genomedatadir, tracks, seqfiles)

if __name__ == "__main__":
    sys.exit(main())

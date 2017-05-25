from __future__ import print_function
import argparse
import sys

from . import __version__
from _util import maybe_gzip_open

from genomedata import Genome
import numpy as np

BED_FILETYPE = "bed"
WIGGLE_FILETYPE = "wig"


def filter_data(gd_filename, track_names, filter_filename,
                filter_threshold=None, is_verbose=False):

    # Get the genomic file type of the filter
    filter_filetype = get_filter_filetype()
    # If no filetype is detected or supported
    if not filter_filetype:
        # Report an error
        raise ValueError("Filter {} file type not "
                         "supported.".format(filter_filename))

    # Open the archive for reading and writing
    with Genome(gd_filename, mode="r+") as genome, \
         maybe_gzip_open(filter_filename) as filter_file:

        # If track names were specified
        if track_names:
            if is_verbose:
                print("Filtering tracks: ", track_names)
            # Get track indexes based on given track names
            track_indexes = [genome.index_continuous(track_name) for track_name
                             in track_names]
        # Otherwise
        else:
            if is_verbose:
                print("Filtering all tracks: ", genome.tracknames_continuous)
            # Get all track indexes
            track_indexes = range(genome.num_tracks_continuous)

        # Create NaN mask based on number of tracks
        nan_mask = np.full(len(track_indexes), np.nan)

        # For every chromosome and filter region
        for chromosome_name, filter_start, filter_end in \
        get_next_genomic_filter_region(filter_file, filter_threshold):

            if is_verbose:
                print("Filtering out region: ", chromosome_name, filter_start,
                      filter_end)

            # XXX: Check if start < end?
            filter_region = slice(filter_start, filter_end)
            chromosome = genome[chromosome_name]

            # Get list of supercontigs from region
            supercontigs = chromosome.supercontigs[filter_region]
            for supercontig in supercontigs:
                # Write in NaNs in region based on track indexes
                supercontig.continuous[filter_region, track_indexes] = nan_mask


def get_filter_filetype(filter_filename):
    if filter_filename.endswith("bed"):
        return BED_FILETYPE
    
    # No filetype detected
    return None


def get_next_genomic_filter_region(filter_file_handle):
    # Assume BED3 for now
    for line in filter_file_handle:
        fields = line.split("\t")
        yield fields[0], int(fields[1]), int(fields[2].rstrip())


def main():
    description = ("Filter TRACKNAME(s) from a genomedata archive with "
                   "FILTERFILE using an optional give threshold value.")
    parser = argparse.ArgumentParser(description=description,
                                     prog="genomedata-filter-data",
                                     version=__version__)

    parser.add_argument("filterfile", help="filter file")
    parser.add_argument('gdarchive', help='genomedata archive')
    parser.add_argument("-t", "--trackname", nargs="+",
                        help="Track(s) to be filtered (default: all)")
    # TODO: Use a threshold string (e.g. ">=0.4") as a parameter and parse
    parser.add_argument("--threshold", default=None, type=float,
                        help="Filter all values less than or equal to a "
                        "specified threshold value (default: no threshold)")
    parser.add_argument("--verbose", default=False, action="store_true",
                        help="Print status and diagnostic messages")

    args = parser.parse_args()

    gd_filename = args.gdarchive
    track_names = args.trackname
    filter_filename = args.filterfile
    filter_threshold = args.threshold
    is_verbose = args.verbose

    filter_data(gd_filename, track_names, filter_filename, filter_threshold,
                is_verbose)


if __name__ == "__main__":
    sys.exit(main())

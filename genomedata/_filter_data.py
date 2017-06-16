from __future__ import print_function
import argparse
from functools import partial
import operator
from re import match
import sys

from . import __version__
from _util import maybe_gzip_open

from genomedata import Genome
import numpy as np

BED_FILETYPE = "bed"
WIGGLE_FILETYPE = "wig"

# NB: operator.op(a,b) is the same as "a op b"
# These operators are partially evaluated with the 'a' value being checked
# against
# E.g. filter out all values <0.5, where the value being considered is "b"
# Return a function where b<0.5 is true or 0.5>b is true
FILTER_OPERATORS = {
    ">=": operator.le,
    ">": operator.lt,
    "<=": operator.ge,
    "<": operator.gt,
    "==": operator.eq,
    "=": operator.eq,
    "!=": operator.ne,
    "~": operator.ne,
}


def filter_data(gd_filename, track_names, filter_filename,
                filter_threshold=None, verbose=False):

    # TODO: Make filter_threshold a function that takes the current value under
    # consideration, where if it evaluates to true, the value given is filtered
    # out

    # Get the genomic file type of the filter
    filter_filetype = get_filter_filetype(filter_filename)
    # If no filetype is detected or supported
    if not filter_filetype:
        # Report an error
        raise ValueError("Filter {} file type not "
                         "supported.".format(filter_filename))

    # Open the archive for reading and writing
    with Genome(gd_filename, mode="r+") as genome, \
            maybe_gzip_open(filter_filename) as filter_file:

        # Create NaN mask based on number of tracks to filter
        num_filter_tracks = len(track_names)
        nan_mask = np.full(num_filter_tracks, np.nan)

        if verbose:
            if track_names:
                print("Filtering tracks: ", track_names)
            else:
                print("Filtering all tracks: ", genome.tracknames_continuous)

        # For every chromosome and filter region
        for chromosome_name, filter_start, filter_end in \
                get_next_genomic_filter_region(filter_file, filter_threshold):

            if verbose:
                print("Filtering out region: ", chromosome_name, filter_start,
                      filter_end)

            # Mask out filter region with NaNs
            chromosome = genome[chromosome_name]
            chromosome[filter_start:filter_end, track_names] = nan_mask


def get_filter_filetype(filter_filename):
    if filter_filename.endswith("bed"):
        return BED_FILETYPE

    # No filetype detected
    return None


def get_next_genomic_filter_region(filter_file_handle, filter_threshold):
    # Assume BED3 for now
    for line in filter_file_handle:
        fields = line.split("\t")
        yield fields[0], int(fields[1]), int(fields[2].rstrip())


def parse_filter_option(filter_option):
    """Gets a operator/value combination as a string and returns a function
    that performs the specified operation (e.g. '>=0.3')
    """
    # Get the value and the operator from the filter string given from options
    re_match = match(r"(?P<operator>\W+)\s*(?P<value>\d+(.\d+)?)",
                     filter_option)
    if not re_match:
        raise ValueError("Could not understand the filter option of "
                         "'{}'".format(filter_option))

    filter_operator = re_match.group("operator")
    filter_value = re_match.group("value")

    try:
        filter_function = partial(FILTER_OPERATORS[filter_operator],
                                  filter_value)
    except KeyError:
        raise ValueError("The operator {} is not understood or "
                         "supported".format(filter_operator))

    return filter_function


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
    # parser.add_argument("--threshold", default=none, type=float,
    #                     help="filter all values less than or equal to a "
    #                     "specified threshold value (default: no threshold)")
    parser.add_argument("--filter",
                        help="Specify a comparison operation on a value to "
                        "filter out (e.g. \"<=0.5\") (default: all values "
                        "filtered)")

    parser.add_argument("--verbose", default=False, action="store_true",
                        help="Print status and diagnostic messages")

    args = parser.parse_args()

    gd_filename = args.gdarchive
    track_names = args.trackname
    filter_filename = args.filterfile
    # filter_threshold = args.threshold
    # If a filter was specified
    if args.filter:
        # Parse the given string to get a filter function
        filter_function = parse_filter_string(args.filter)
    is_verbose = args.verbose

    filter_data(gd_filename, track_names, filter_filename, filter_function,
                is_verbose)


if __name__ == "__main__":
    sys.exit(main())

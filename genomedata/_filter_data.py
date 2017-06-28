from __future__ import print_function
import argparse
from functools import partial
import operator
from re import match
import sys

from . import __version__
from _close_data import write_genome_metadata
from _util import EXT_GZ, maybe_gzip_open
from _filter_data_parsers import get_bed_filter_region, get_wig_filter_region

from genomedata import Genome
import numpy as np

BED_FILETYPE = "bed"
WIGGLE_FILETYPE = "wig"

# Find which filter generator to call based on filetype
FILTER_REGION_GENERATORS = {
    BED_FILETYPE: get_bed_filter_region,
    WIGGLE_FILETYPE: get_wig_filter_region
}

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


def filter_data(gd_filename, filter_filename, track_names=None,
                filter_function=None, verbose=False, keep_archive_open=False,
                dry_run=False):

    # Get the genomic file type of the filter
    filter_filetype = get_filter_filetype(filter_filename)
    get_next_genomic_filter_region = FILTER_REGION_GENERATORS[filter_filetype]

    # Open the archive for reading and writing
    with Genome(gd_filename, mode="r+") as genome, \
            maybe_gzip_open(filter_filename) as filter_file:

        # Create NaN mask based on number of tracks to filter
        if track_names:
            num_filter_tracks = len(track_names)
            if verbose:
                print("Tracks for filtering: ", track_names)
        else:
            num_filter_tracks = genome.num_tracks_continuous
            if verbose:
                print("All tracks selected for filtering: ",
                      genome.tracknames_continuous)

        nan_mask = np.full(num_filter_tracks, np.nan)

        # For every chromosome and filter region
        for chromosome_name, filter_start, filter_end in \
                get_next_genomic_filter_region(filter_file, filter_function):

            if verbose:
                print("Filtering out region: ", chromosome_name, filter_start,
                      filter_end)

            # Mask out filter region with NaNs
            chromosome = genome[chromosome_name]
            if not dry_run:
                chromosome[filter_start:filter_end, track_names] = nan_mask

        # If the archive should be closed
        if not keep_archive_open:
            # Close the archive
            write_genome_metadata(genome, verbose)


def get_filter_filetype(filter_filename):
    filename = filter_filename
    filter_filetype = None

    if filter_filename.endswith("." + EXT_GZ):
        # Find the extension before gzip extension
        filename = filter_filename[:-3]  # Remove ".gz"

    if filename.endswith("bed"):
        filter_filetype = BED_FILETYPE
    elif (filename.endswith("wg") or
          filename.endswith("wig") or
          filename.endswith("wigVar")):
        filter_filetype = WIGGLE_FILETYPE

    # If no known filetype detected
    if not filter_filetype:
        # Report an error
        raise ValueError("Filter {} file type not "
                         "supported.".format(filter_filename))

    # Otherwise return the filetype
    return filter_filetype


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

    # If the given operator is supported
    # And if the the filter can be converted to a number
    try:
        # Attempt to convert the filter value into a number
        filter_value = float(filter_value)
        # Get the function associated with this operator and given value
        filter_function = partial(FILTER_OPERATORS[filter_operator],
                                  filter_value)
    # Otherwise display an error about the operator
    except KeyError:
        raise ValueError("The operator {} is not understood or "
                         "supported".format(filter_operator))
    # Otherwise display an error about the filter value
    except ValueError:
        # XXX: May not be necessary since the regex captures digits only
        raise ValueError("Could not convert filter value {} to a "
                         "number".format(filter_value))

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
    parser.add_argument("--filter",
                        help="Specify a comparison operation on a value to "
                        "filter out (e.g. \"<0.5\" will remove all values less "
                        "than 0.5) (default: all values filtered)")
    parser.add_argument("--no-close", default=False, action="store_true",
                        help="Do not close the genomedata archive after "
                        "filtering.")
    parser.add_argument("--dry-run", default=False, action="store_true",
                        help="Do not perform any filtering. Useful with "
                        "verbosity set to see what regions would be filtered")
    parser.add_argument("--verbose", default=False, action="store_true",
                        help="Print status and diagnostic messages")

    args = parser.parse_args()

    gd_filename = args.gdarchive
    track_names = args.trackname
    filter_filename = args.filterfile
    is_dry_run = args.dry_run
    keep_archive_open = args.no_close
    is_verbose = args.verbose

    # If a filter was specified
    if args.filter:
        # Parse the given string to get a filter function
        filter_function = parse_filter_option(args.filter)
    else:
        # Otherwise use no filter function
        filter_function = None

    filter_data(gd_filename, filter_filename, track_names, filter_function,
                is_verbose, keep_archive_open, is_dry_run)


if __name__ == "__main__":
    sys.exit(main())

from __future__ import print_function
import argparse
from functools import partial
import operator
from os.path import splitext
from re import match
import sys

from . import __version__
from _close_data import write_genome_metadata
from _util import EXT_GZ, maybe_gzip_open
from _filter_data_parsers import (get_bed_filter_region, get_wig_filter_region,
                                  merged_filter_region_generator)
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
HARD_MASK_OPERATORS = {
    "ge": operator.le,
    "gt": operator.lt,
    "le": operator.ge,
    "lt": operator.gt,
    "eq": operator.eq,
    "ne": operator.ne,
}


def hard_mask_data(gd_filename, hard_mask_filename, track_names=None,
                   hard_mask_function=None, verbose=False,
                   keep_archive_open=False, dry_run=False):

    # Get the genomic file type of the mask
    hard_mask_filetype = get_hard_mask_filetype(hard_mask_filename)
    get_next_genomic_mask_region = \
        FILTER_REGION_GENERATORS[hard_mask_filetype]

    # Open the archive for reading and writing
    with Genome(gd_filename, mode="r+") as genome, \
            maybe_gzip_open(hard_mask_filename) as hard_mask_file:

        # Create NaN mask based on number of tracks to mask
        if track_names:
            num_mask_tracks = len(track_names)
            if verbose:
                print("Tracks for masking: ", track_names)
        else:
            num_mask_tracks = genome.num_tracks_continuous
            if verbose:
                print("All tracks selected for masking: ",
                      genome.tracknames_continuous)

        nan_mask = np.full(num_mask_tracks, np.nan)

        merged_filter_regions = merged_filter_region_generator(
                                    get_next_genomic_mask_region,
                                    hard_mask_file,
                                    hard_mask_function)

        # For every chromosome and filter region
        for chromosome_name, mask_start, mask_end in merged_filter_regions:
            if verbose:
                print("Masking out region: ", chromosome_name, mask_start,
                      mask_end)

            # Mask out mask region with NaNs
            chromosome = genome[chromosome_name]
            if not dry_run:
                # If specified track names
                if track_names:
                    # Mask based on given tracks
                    chromosome[mask_start:mask_end, track_names] = nan_mask
                # Otherwise
                else:
                    # Mask out all tracks
                    chromosome[mask_start:mask_end] = nan_mask

        # If the archive should be closed
        if not keep_archive_open:
            # Close the archive
            write_genome_metadata(genome, verbose)


def get_hard_mask_filetype(hard_mask_filename):
    filename = hard_mask_filename.lower()
    hard_mask_filetype = None

    root_filename, file_extension = splitext(filename)
    if file_extension == "." + EXT_GZ:
        # Find the extension before gzip extension
        root_filename, file_extension = splitext(root_filename)

    if file_extension == ".bed":
        hard_mask_filetype = BED_FILETYPE
    elif (file_extension == ".wg" or
          file_extension == ".wig" or
          file_extension == ".wigVar"):
        hard_mask_filetype = WIGGLE_FILETYPE

    # If no known filetype detected
    if not hard_mask_filetype:
        # Report an error
        raise ValueError("Mask {} file type not "
                         "supported.".format(hard_mask_filename))

    # Otherwise return the filetype
    return hard_mask_filetype


def parse_hard_mask_option(mask_option):
    """Gets a operator/value combination as a string and returns a function
    that performs the specified operation (e.g. '>=0.3')
    """
    # Get the value and the operator from the filter string given from options
    # "operator" group matches on exactly two lowercase letters
    # "value" group matches a digital that is an optional float
    # There is an optional space between groups
    re_match = match(r"^(?P<operator>[a-z]{2})\s*(?P<value>\d+(.\d+)?)$",
                     mask_option)
    if not re_match:
        raise ValueError("Could not understand the hard mask option of "
                         "'{}'".format(mask_option))

    mask_operator = re_match.group("operator")
    mask_value = re_match.group("value")

    # If the given operator is supported
    # And if the the filter can be converted to a number
    try:
        # Attempt to convert the filter value into a number
        mask_value = float(mask_value)
        # Get the function associated with this operator and given value
        hard_mask_function = partial(HARD_MASK_OPERATORS[mask_operator],
                                     mask_value)
    # Otherwise display an error about the operator
    except KeyError:
        raise ValueError("The operator {} is not understood or "
                         "supported".format(mask_operator))
    # Otherwise display an error about the filter value
    except ValueError:
        # XXX: May not be necessary since the regex captures digits only
        raise ValueError("Could not convert filter value {} to a "
                         "number".format(mask_value))

    return hard_mask_function


def main():
    description = ("Permanently mask TRACKNAME(s) from a genomedata archive with "
                   "MASKFILE using an optional filter operator.")
    parser = argparse.ArgumentParser(description=description,
                                     prog="genomedata-hard-mask",
                                     version=__version__)

    parser.add_argument("maskfile", help="input mask file")
    parser.add_argument('gdarchive', help='genomedata archive')
    parser.add_argument("-t", "--trackname", nargs="+",
                        help="Track(s) to be filtered (default: all)")
    parser.add_argument("--hardmask",
                        metavar="OPERATOR",
                        help="Specify a comparison operation on a value to "
                        "mask out (e.g. \"lt0.5\" will mask all values less "
                        "than 0.5). See the bash comparison operators for the "
                        "two letter operations (default: all values masked)")
    parser.add_argument("--no-close", default=False, action="store_true",
                        help="Do not close the genomedata archive after "
                        "masking")
    parser.add_argument("--dry-run", default=False, action="store_true",
                        help="Do not perform any masking. Useful with "
                        "verbosity set to see what regions would be filtered")
    parser.add_argument("--verbose", default=False, action="store_true",
                        help="Print status and diagnostic messages")

    args = parser.parse_args()

    gd_filename = args.gdarchive
    track_names = args.trackname
    mask_filename = args.maskfile
    is_dry_run = args.dry_run
    keep_archive_open = args.no_close
    is_verbose = args.verbose

    # If a filter was specified
    if args.hardmask:
        # Parse the given string to get a hard mask function
        hard_mask_function = parse_hard_mask_option(args.hardmask)
    else:
        # Otherwise use no hard mask function
        hard_mask_function = None

    hard_mask_data(gd_filename, mask_filename, track_names,
                   hard_mask_function, is_verbose, keep_archive_open,
                   is_dry_run)


if __name__ == "__main__":
    sys.exit(main())

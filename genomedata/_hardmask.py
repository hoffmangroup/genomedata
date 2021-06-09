#!/usr/bin/env python
from __future__ import absolute_import, print_function

"""
_hardmask.py: A python interface to mask out track data from a genomedata
archive
"""

import argparse
from functools import partial
import operator
from os.path import splitext
from re import match
import sys

from numpy import full, nan

from ._close_data import write_genome_metadata
from ._util import EXT_GZ, maybe_gzip_open
from ._hardmask_parsers import (get_bed_filter_region, get_wig_filter_region,
                                merged_filter_region_generator)
from . import Genome, __version__

BED_FILETYPE = "bed"
BED_SUFFIXES = frozenset({".bed"})

WIGGLE_FILETYPE = "wig"
WIGGLE_SUFFIXES = frozenset({".wg", ".wig", ".wigVar", ".wigFix"})

SUFFIX_GZ = "." + EXT_GZ

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
HARDMASK_OPERATORS = {
    "ge": operator.le,
    "gt": operator.lt,
    "le": operator.ge,
    "lt": operator.gt,
    "eq": operator.eq,
    "ne": operator.ne,
}


def hardmask_data(gd_filename, hardmask_filename, track_names=None,
                  hardmask_function=None, verbose=False,
                  keep_archive_open=False, dry_run=False):

    # Get the genomic file type of the mask
    hardmask_filetype = get_hardmask_filetype(hardmask_filename)
    get_next_genomic_mask_region = \
        FILTER_REGION_GENERATORS[hardmask_filetype]

    # Open the archive for reading and writing
    with Genome(gd_filename, mode="r+") as genome, \
            maybe_gzip_open(hardmask_filename) as hardmask_file:

        # Create NaN mask based on number of tracks to mask
        if track_names:
            num_mask_tracks = len(track_names)
            if verbose:
                print("Tracks for masking: ", track_names, file=sys.stderr)
        else:
            num_mask_tracks = genome.num_tracks_continuous
            if verbose:
                print("All tracks selected for masking: ",
                      genome.tracknames_continuous, file=sys.stderr)

        nan_mask = full(num_mask_tracks, nan)

        merged_filter_regions = merged_filter_region_generator(
                                    get_next_genomic_mask_region,
                                    hardmask_file,
                                    hardmask_function)

        # For every chromosome and filter region
        for chromosome_name, mask_start, mask_end in merged_filter_regions:
            if verbose:
                print("Masking out region: ", chromosome_name, mask_start,
                      mask_end, file=sys.stderr)

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


def get_hardmask_filetype(hardmask_filename):
    filename = hardmask_filename.lower()
    hardmask_filetype = None

    root_filename, file_extension = splitext(filename)
    if file_extension == SUFFIX_GZ:
        # Find the extension before gzip extension
        root_filename, file_extension = splitext(root_filename)

    if file_extension in BED_SUFFIXES:
        hardmask_filetype = BED_FILETYPE
    elif file_extension in WIGGLE_SUFFIXES:
        hardmask_filetype = WIGGLE_FILETYPE
    # If no known filetype detected
    else:
        # Report an error
        raise ValueError("Mask {} file type not "
                         "supported.".format(hardmask_filename))

    # Otherwise return the filetype
    return hardmask_filetype


def parse_hardmask_option(mask_option):
    """Gets a operator/value combination as a string and returns a function
    that performs the specified operation (e.g. 'ge0.3')
    """
    # Get the value and the operator from the filter string given from options
    # "operator" group matches on exactly two lowercase letters
    # "value" group matches a digital that is an optional float
    # There is an optional space between groups
    re_match = match(r"^(?P<operator>[a-z]{2})\s*(?P<value>\d+(.\d+)?)$",
                     mask_option)
    if not re_match:
        raise ValueError("Could not understand the hardmask option of "
                         "'{}'".format(mask_option))

    mask_operator = re_match.group("operator")
    mask_value = re_match.group("value")

    # If the given operator is supported
    # And if the the filter can be converted to a number
    try:
        # Attempt to convert the filter value into a number
        mask_value = float(mask_value)
        # Get the function associated with this operator and given value
        hardmask_function = partial(HARDMASK_OPERATORS[mask_operator],
                                    mask_value)
    # Otherwise display an error about the operator
    except KeyError:
        raise ValueError("The operator {} is not understood or "
                         "supported".format(mask_operator))

    return hardmask_function


def main():
    description = ("Permanently mask TRACKNAME(s) from a genomedata archive "
                   "with MASKFILE using an optional filter operator.")
    parser = argparse.ArgumentParser(description=description,
                                     prog="genomedata-hard-mask")

    parser.add_argument('--version', action='version', version=__version__)

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
        # Parse the given string to get a hardmask function
        hardmask_function = parse_hardmask_option(args.hardmask)
    else:
        # Otherwise use no hardmask function
        hardmask_function = None

    hardmask_data(gd_filename, mask_filename, track_names,
                  hardmask_function, is_verbose, keep_archive_open,
                  is_dry_run)


if __name__ == "__main__":
    sys.exit(main())

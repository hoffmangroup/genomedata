#!/usr/bin/env python
from __future__ import division, with_statement

"""
_close_data: DESCRIPTION
"""

__version__ = "$Revision$"

# Copyright 2008-2014 Michael M. Hoffman <mhoffman@uhnresearch.ca>

import sys

from numpy import (amin, amax, array, diff, hstack, isfinite,
                   NINF, PINF, square)
from tables import NoSuchNodeError

from . import Genome
from ._load_seq import MIN_GAP_LEN
from ._util import fill_array, init_num_obs, new_extrema

def update_extrema(func, extrema, data, col_index):
    extrema[col_index] = new_extrema(func, data, extrema[col_index])

def find_chunks(mask_rows_any_present):
    """
    find chunks that have less than MIN_GAP_LEN missing data
    gaps in a row
    """
    # defaults
    starts = array([0])
    ends = array([mask_rows_any_present.shape[0]])

    # get all of the indices where there is any data
    indices_present = mask_rows_any_present.nonzero()[0]

    if len(indices_present):
        # make a mask of whether the difference from one index to the
        # next is >= MIN_GAP_LEN
        diffs_signif = diff(indices_present) >= MIN_GAP_LEN

        # convert the mask back to indices of the original indices
        # these are the indices immediately before a big gap
        indices_signif = diffs_signif.nonzero()[0]

        if len(indices_signif):
            # start with the indices immediately after a big gap
            starts = indices_present[hstack([0, indices_signif+1])]

            # end with indices immediately before a big gap
            ends_inclusive = indices_present[hstack([indices_signif, -1])]

            # add 1 to ends because we want slice(start, end) to
            # include the last_index; convert from inclusive (closed)
            # to exclusive (half-open) coordinates, as Python needs
            ends = ends_inclusive + 1

    return starts, ends

def write_metadata(chromosome, verbose=False):
    if verbose:
        print >>sys.stderr, "writing metadata for %s" % chromosome

    if chromosome.num_tracks_continuous == 0:
        chromosome.attrs.dirty = False
        return

    h5file = chromosome.h5file
    tracknames = chromosome.tracknames_continuous

    num_obs = len(tracknames)
    row_shape = (num_obs,)
    mins = fill_array(PINF, row_shape)
    maxs = fill_array(NINF, row_shape)
    sums = fill_array(0.0, row_shape)
    sums_squares = fill_array(0.0, row_shape)
    num_datapoints = fill_array(0, row_shape)

    for supercontig in chromosome:
        if verbose:
            print >>sys.stderr, " scanning %s" % supercontig

        try:
            continuous = supercontig.continuous
        except NoSuchNodeError:
            raise NoSuchNodeError("Supercontig found missing continuous")

        # only runs when assertions checked
        if __debug__:
            init_num_obs(num_obs, continuous) # for the assertion

        mask_rows_any_present = fill_array(False, continuous.shape[0])

        # doing this column by column greatly reduces the memory
        # footprint when you have large numbers of tracks. It also
        # simplifies the logic for the summary stats, since you don't
        # have to change the mask value for every operation, like in
        # revisions <= r243
        for col_index, trackname in enumerate(tracknames):
            if verbose:
                print >>sys.stderr, "  %s" % trackname

            ## read data
            col = continuous[:, col_index]

            mask_present = isfinite(col)
            mask_rows_any_present[mask_present] = True
            col_finite = col[mask_present]
            del col  # col not needed anymore (optimization)

            num_datapoints_col = len(col_finite)
            if num_datapoints_col:
                update_extrema(amin, mins, col_finite, col_index)
                update_extrema(amax, maxs, col_finite, col_index)

                sums[col_index] += col_finite.sum(0)
                sums_squares[col_index] += square(col_finite).sum(0)
                num_datapoints[col_index] += num_datapoints_col

        supercontig_attrs = supercontig.attrs

        starts, ends = find_chunks(mask_rows_any_present)
        supercontig_attrs.chunk_starts = starts
        supercontig_attrs.chunk_ends = ends

    chromosome_attrs = chromosome.attrs
    chromosome_attrs.mins = mins
    chromosome_attrs.maxs = maxs
    chromosome_attrs.sums = sums
    chromosome_attrs.sums_squares = sums_squares
    chromosome_attrs.num_datapoints = num_datapoints
    chromosome_attrs.dirty = False

def close_data(gdfilename, verbose=False):
    with Genome(gdfilename, mode="r+") as genome:
        for chromosome in genome:
            if chromosome.attrs.dirty:
                write_metadata(chromosome, verbose=verbose)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... GENOMEDATAFILE"
    version = "%%prog %s" % __version__
    description = ("Compute summary statistics for data in Genomedata archive"
                   " and ready for accessing.")
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-v", "--verbose", dest="verbose",
                      default=False, action="store_true",
                      help="Print status updates and diagnostic messages")

    options, args = parser.parse_args(args)

    if not len(args) == 1:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    gdfilename = args[0]
    kwargs = {"verbose": options.verbose}
    return close_data(gdfilename, **kwargs)

if __name__ == "__main__":
    sys.exit(main())

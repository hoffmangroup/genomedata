#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
_close_data: DESCRIPTION
"""

# Copyright 2008-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

from argparse import ArgumentParser
import sys

from numpy import (amin, amax, argmax, array, diff, hstack, isfinite, inf,
                   square)
from six.moves import zip
from tables import NoSuchNodeError

from . import Genome, __version__
from ._load_seq import MIN_GAP_LEN
from ._util import fill_array, init_num_obs, new_extrema


def update_extrema(func, extrema, data, col_index):
    extrema[col_index] = new_extrema(func, data, extrema[col_index])


def find_chunk_gaps_in_supercontig(mask_rows_any_present):
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
            starts = indices_present[hstack([0, indices_signif + 1])]

            # end with indices immediately before a big gap
            ends_inclusive = indices_present[hstack([indices_signif, -1])]

            # add 1 to ends because we want slice(start, end) to
            # include the last_index; convert from inclusive (closed)
            # to exclusive (half-open) coordinates, as Python needs
            ends = ends_inclusive + 1

    return starts, ends


def find_next_present_coord(chromosome, start_supercontig):
    """
    get next chromosomal coordinate where data is present from starting
    supercontig onwards. Returns None if no data is found or the starting
    supercontig does not exist
    """
    res = None

    # If our starting supercontig exists
    if start_supercontig:
        start_coordinate = chromosome.start + start_supercontig.start

        # For each supercontig between our starting supercontig to end of
        # chromosome
        supercontigs = chromosome.supercontigs[start_coordinate:chromosome.end]
        for supercontig in supercontigs:
            continuous = supercontig.continuous
            tracknames = chromosome.tracknames_continuous

            # For each track in supercontig
            for col_index, trackname in enumerate(tracknames):
                col = continuous[:, col_index]
                mask_present = isfinite(col)
                # If there is any value present
                if any(mask_present):  # mask_present is a bool "present" array
                    # Get the earliest coordinate
                    # NB: argmax gives the index of the maximum value. In the
                    # case of ties (which there will be) argmax returns the
                    # first index with the maximum value. I believe this is the
                    # recommended (and fastest?) way of finding the first
                    # occurrence of a value
                    first_present_coord = (supercontig.start +
                                           argmax(mask_present))
                    if not res:
                        res = first_present_coord
                    else:
                        res = min(res, first_present_coord)

            # If a present value was found
            if res:
                # Stop searching
                break

    return res


def trim_chunks_start(prev_supercontig, prev_end_coord,
                      supercontig, presence_mask,
                      chunk_start_index):
    """Determine the first start index for chunks in this supercontig"""

    # NB: this is unchanged if the distance between the previous
    # end coordinates and the new start coordinate is less than
    # MIN_GAP_LEN
    res = chunk_start_index

    # Find the coordinate where data is first present in this
    # supercontig
    # NB: argmax will give the index of the first element on matching
    # winners of the boolean array. The winners are "present" values
    first_present_value_index = argmax(presence_mask)
    # If a previous supercontig exists
    if prev_supercontig:
        # Get previous end chunk chromosomal coordinate
        first_present_value_coord = supercontig.start + \
            first_present_value_index
        # If the distance between previous and current present data
        # chromsomal positions is greater than MIN_GAP_LEN
        if (first_present_value_coord - prev_end_coord >
                MIN_GAP_LEN):
            # Truncate the start chunk of the supercontig to start
            # where data is present
            res = first_present_value_index
    # Else this is the first supercontig
    else:
        # Truncate the start chunk of the supercontig to start where
        # data is present
        res = first_present_value_index

    return res


def trim_chunks_end(next_supercontig, next_start_coord,
                    supercontig, presence_mask,
                    chunk_end_index):
    """Determine the first start index for chunks in this supercontig"""

    # NB: this is unchanged if the distance between the last
    # present coordinate and the next present coordinate is less
    # than MIN_GAP_LEN
    res = chunk_end_index

    # NB: argmax will give the index of the first element on matching
    # winners of the boolean array. The winners are "present" values
    last_present_index = (len(presence_mask) -
                          argmax(presence_mask[::-1]))

    # If this is not the last supercontig
    if next_supercontig:
        last_present_chrom_coord = (supercontig.start +
                                    last_present_index)
        # If the next present value is greater than MIN_GAP_LEN away
        # or there is no next possible present value
        if (not next_start_coord or
            next_start_coord - last_present_chrom_coord >
                MIN_GAP_LEN):
            # Truncate the chunk end to where data is last present in
            # this supercontig
            res = last_present_index

    # Else this is the last supercontig
    else:
        # Truncate the last end chunk of the supercontig to end where
        # data is last present
        res = last_present_index

    return res


def write_metadata(chromosome, verbose=False):
    if verbose:
        print("writing metadata for %s" % chromosome, file=sys.stderr)

    if chromosome.num_tracks_continuous == 0:
        chromosome.attrs.dirty = False
        return

    tracknames = chromosome.tracknames_continuous

    num_obs = len(tracknames)
    row_shape = (num_obs,)
    mins = fill_array(inf, row_shape)
    maxs = fill_array(-inf, row_shape)
    sums = fill_array(0.0, row_shape)
    sums_squares = fill_array(0.0, row_shape)
    num_datapoints = fill_array(0, row_shape)

    supercontigs = chromosome.supercontigs[chromosome.start:chromosome.end]
    prev_chrom_end = 0  # keeps track of where data was last present

    prev_supercontigs = [None] + supercontigs[:-1]
    next_supercontigs = supercontigs[1:] + [None]

    zipper = zip(prev_supercontigs, supercontigs, next_supercontigs)

    for prev_supercontig, supercontig, next_supercontig in zipper:

        if verbose:
            print(" scanning %s" % supercontig, file=sys.stderr)

        try:
            continuous = supercontig.continuous
        except NoSuchNodeError:
            raise NoSuchNodeError("Supercontig found missing continuous")

        # only runs when assertions checked
        if __debug__:
            init_num_obs(num_obs, continuous)  # for the assertion

        # An array containing True when any data is present across all tracks
        mask_rows_any_present = fill_array(False, continuous.shape[0])

        # doing this column by column greatly reduces the memory
        # footprint when you have large numbers of tracks. It also
        # simplifies the logic for the summary stats, since you don't
        # have to change the mask value for every operation, like in
        # revisions <= r243
        for col_index, trackname in enumerate(tracknames):
            if verbose:
                print("  %s" % trackname, file=sys.stderr)

            # read data
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

        next_chrom_start = find_next_present_coord(chromosome,
                                                   next_supercontig)

        # If any data is present in this supercontig
        if any(mask_rows_any_present):
            # Find all gaps between present data
            starts, ends = find_chunk_gaps_in_supercontig(
                               mask_rows_any_present
                           )

            # Trim the starting chunk index for this supercontig
            starts[0] = trim_chunks_start(prev_supercontig, prev_chrom_end,
                                          supercontig, mask_rows_any_present,
                                          starts[0])

            # Trim the end chunk index for this supercontig
            ends[-1] = trim_chunks_end(next_supercontig, next_chrom_start,
                                       supercontig, mask_rows_any_present,
                                       ends[-1])

            # Update our new previously present coordinate for next iteration
            prev_chrom_end = supercontig.start + ends[-1]
        # Otherwise there isn't any data present for this contig
        else:
            # If this is the last supercontig
            if not next_supercontig:
                # Exclude the entire region
                starts = array([])
                ends = array([])
            # Otherwise there are supercontigs to lookahead at
            else:
                # If the distance between the previously last present value
                # And next possible present value is less than MIN_GAP_LEN
                if (next_chrom_start and
                   next_chrom_start - prev_chrom_end <
                   MIN_GAP_LEN):
                    # Include the entire region
                    starts = array([0])
                    ends = array([mask_rows_any_present.shape[0]])
                # Otherwise
                else:
                    # Exclude the entire region
                    starts = array([])
                    ends = array([])

        supercontig_attrs.chunk_starts = starts
        supercontig_attrs.chunk_ends = ends

    chromosome_attrs = chromosome.attrs
    chromosome_attrs.mins = mins
    chromosome_attrs.maxs = maxs
    chromosome_attrs.sums = sums
    chromosome_attrs.sums_squares = sums_squares
    chromosome_attrs.num_datapoints = num_datapoints
    chromosome_attrs.dirty = False


def write_genome_metadata(genome, verbose):
    for chromosome in genome:
        if chromosome.attrs.dirty:
            write_metadata(chromosome, verbose=verbose)


def close_data(gdfilename, verbose=False):
    with Genome(gdfilename, mode="r+") as genome:
        write_genome_metadata(genome, verbose)


def parse_options(args):
    description = ("Compute summary statistics for data in Genomedata archive"
                   " and ready for accessing.")

    parser = ArgumentParser(
        description=description,
        prog='genomedata-close-data')

    parser.add_argument('--version', action='version', version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')

    parser.add_argument(
        "--verbose",
        default=False,
        action="store_true",
        help="Print status updates and diagnostic messages")

    args = parser.parse_args(args)

    return args


def main(argv=sys.argv[1:]):
    args = parse_options(argv)
    gdarchive = args.gdarchive
    return close_data(gdarchive, verbose=args.verbose)


if __name__ == "__main__":
    sys.exit(main())

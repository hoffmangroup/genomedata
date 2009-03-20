#!/usr/bin/env python
from __future__ import division, with_statement

"""
save_metadata: DESCRIPTION
"""

__version__ = "$Revision$"

# Copyright 2008-2009 Michael M. Hoffman <mmh1@washington.edu>

import sys

from numpy import amin, amax, append, array, diff, isfinite, NINF, PINF, square
from tables import openFile

from .load_seq import MIN_GAP_LEN
from ._util import (fill_array, get_tracknames, init_num_obs, new_extrema,
                    walk_continuous_supercontigs)


def get_nucleotide_index (nuc, nucs):
# Mirela: get the index of the specified nucleotide nuc in the array nucs
    assert len(nucs) == 3    
    if (nuc.tostring().upper() in ('A','T')):
        return 0
    if (nuc.tostring().upper() in ('C','G')):
        return 1
    return 2    
                
def get_dinucleotide_index (dinuc1, dinuc2, dinucs):
# Mirela: get the index of the specified nucleotides dinuc1 and dinuc2 in the array dinucs
    assert len (dinucs) == 11
    d1 = dinuc1.tostring().upper()
    d2 = dinuc2.tostring().upper()
    if ((d1 == 'A' and d2 == 'A') or (d1 == 'T' and d2 == 'T')): return 0
    if ((d1 == 'A' and d2 == 'C') or (d1 == 'G' and d2 == 'T')): return 1
    if ((d1 == 'A' and d2 == 'G') or (d1 == 'C' and d2 == 'T')): return 2
    if  (d1 == 'A' and d2 == 'T'): return 3
    if ((d1 == 'C' and d2 == 'A') or (d1 == 'T' and d2 == 'G')): return 4
    if ((d1 == 'C' and d2 == 'C') or (d1 == 'G' and d2 == 'G')): return 5
    if  (d1 == 'C' and d2 == 'G'): return 6
    if ((d1 == 'G' and d2 == 'A') or (d1 == 'T' and d2 == 'C')): return 7
    if  (d1 == 'G' and d2 == 'C'): return 8
    if  (d1 == 'T' and d2 == 'A'): return 9
    return 10


def update_extrema(func, extrema, data, col_index):
    extrema[col_index] = new_extrema(func, data, extrema[col_index])


def write_metadata(chromosome):
    print >>sys.stderr, "writing metadata for %s" % chromosome.title

    tracknames = get_tracknames(chromosome)
    num_obs = len(tracknames)
    row_shape = (num_obs,)
    mins = fill_array(PINF, row_shape)
    maxs = fill_array(NINF, row_shape)
    sums = fill_array(0.0, row_shape)
    sums_squares = fill_array(0.0, row_shape)
    num_datapoints = fill_array(0, row_shape)
    
    # Mirela: add the counts for nucleotides and dinucleotides
    seq_counts_names = array(['A|T', 'C|G', 'ambig_nucleotide', 'AA|TT', 'AC|GT', 'AG|CT', 'AT', 'CA|TG', 'CC|GG', 'CG', 'GA|TC', 'GC', 'TA', 'ambig_dinucleotide'])
    seq_counts = fill_array(0, len(seq_counts_names))
    seq_total_counts = 0        # the total number of nucleotides in this chromosome (only what is in supercontigs)

    for supercontig, continuous in walk_continuous_supercontigs(chromosome):
        print >>sys.stderr, " scanning %s" % supercontig._v_name

        # Mirela: compute the nucleotide counts
        nuc_categories = seq_counts_names[0:3]        
        for nucleotide_int in supercontig.seq.read():            
            nuc_index = get_nucleotide_index(nucleotide_int, nuc_categories)
            seq_counts[nuc_index] += 1      
        
        # Mirela: compute the dinucleotide counts
        dinuc_categories = seq_counts_names[3:]
        for i in range(len(supercontig.seq.read())-1):
            dinuc_index = get_dinucleotide_index(supercontig.seq.read(i), supercontig.seq.read(i+1), dinuc_categories)
            seq_counts[len(nuc_categories)+dinuc_index] += 1
                
        # Mirela: also store the total number of nucleotides in all supercontigs 
        seq_total_counts += len(supercontig.seq.read())

        # only runs when assertions checked
        if __debug__:
            init_num_obs(num_obs, continuous) # for the assertion

        num_rows = continuous.shape[0]
        mask_rows_any_present = fill_array(False, num_rows)

        # doing this column by column greatly reduces the memory
        # footprint when you have large numbers of tracks. It also
        # simplifies the logic for the summary stats, since you don't
        # have to change the mask value for every operation, like in
        # revisions <= r243
        for col_index, trackname in enumerate(tracknames):
            print >>sys.stderr, "  %s" % trackname

            ## read data
            col = continuous[:, col_index]

            mask_present = isfinite(col)
            mask_rows_any_present[mask_present] = True
            col_finite = col[mask_present]
            # XXXopt: should be able to overwrite col, not needed anymore

            num_datapoints_col = len(col_finite)
            if num_datapoints_col:
                update_extrema(amin, mins, col_finite, col_index)
                update_extrema(amax, maxs, col_finite, col_index)

                sums[col_index] += col_finite.sum(0)
                sums_squares[col_index] += square(col_finite).sum(0)
                num_datapoints[col_index] += num_datapoints_col

        ## find chunks that have less than MIN_GAP_LEN missing data
        ## gaps in a row

        # get all of the indices where there is any data
        indices_present = mask_rows_any_present.nonzero()[0]

        if not len(indices_present):
            # remove continuous of empty supercontigs
            continuous._f_remove()
            continue

        # make a mask of whether the difference from one index to the
        # next is >= MIN_GAP_LEN
        diffs_signif = diff(indices_present) >= MIN_GAP_LEN

        # convert the mask back to indices of the original indices
        indices_signif = diffs_signif.nonzero()[0]

        if len(indices_signif):
            starts = indices_present[indices_signif]

            # finish with the index immediately before each start, and the
            # last index
            ends = indices_present[append(indices_signif[1:]-1, -1)]

            # add 1 because we want slice(start, end) to include the
            # last_index
            ends += 1
        else:
            starts = array([0])
            ends = array([num_rows])

        supercontig_attrs = supercontig._v_attrs
        supercontig_attrs.chunk_starts = starts
        supercontig_attrs.chunk_ends = ends

    chromosome_attrs = chromosome.root._v_attrs
    chromosome_attrs.mins = mins
    chromosome_attrs.maxs = maxs
    chromosome_attrs.sums = sums
    chromosome_attrs.sums_squares = sums_squares
    chromosome_attrs.num_datapoints = num_datapoints
    chromosome_attrs.dirty = False
    chromosome_attrs.seq_counts = seq_counts
    chromosome_attrs.seq_counts_names = seq_counts_names
    chromosome_attrs.seq_total_counts = seq_total_counts
    
    print seq_counts

def save_metadata(*filenames):
    for filename in filenames:
        with openFile(filename, "r+") as chromosome:
            write_metadata(chromosome)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... FILE..."
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    options, args = parser.parse_args(args)

    if not len(args) >= 1:
        parser.print_usage()
        sys.exit(1)

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    return save_metadata(*args)

if __name__ == "__main__":
    sys.exit(main())

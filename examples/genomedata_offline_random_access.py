#!/usr/bin/env python
"""
Looks up the data value for the given tracks at each position on stdin.

Usage: ./genomedata_random_access.py GENOMEDATADIR TRACKNAME...

Expects lines on stdin of the form: chrom<whitespace>index
"""

from __future__ import with_statement, division

import sys
import warnings
from collections import defaultdict
from genomedata import Genome

indices = defaultdict(list)
for line in sys.stdin:
    tokens = line.strip().split()
    if tokens:
        chrom, index = line.strip().split()
        indices[chrom].append(int(index))
indices = dict(indices)  # remove defaultdict behavior

with Genome(sys.argv[1]) as genome:
    tracknames = sys.argv[2:]
    for trackname in tracknames:
        assert trackname in genome.tracknames_continuous

    warnings.simplefilter("ignore")
    count = 0
    for chrom in indices:
        indices[chrom].sort()  # Sort by index ascending
        index_iter = iter(indices[chrom])
        index = index_iter.next()
        chromosome = genome[chrom]
        chrom_done = False
        for supercontig, continuous in chromosome.itercontinuous():
            start = supercontig.start
            end = supercontig.end
            supercontig_indices = []
            if not chrom_done:
                try:
                    while index < start:
                        for trackname in tracknames:
                            print "nan"
                            index = index_iter.next()
                    while index < end:
                        supercontig_indices.append(index - start)
                        index = index_iter.next()
                except StopIteration:
                    chrom_done = True

            for trackname in tracknames:
                col_index = chromosome.index_continuous(trackname)
                data_col = continuous[:, col_index]
                for row_index in supercontig_indices:
                    print "%s" % data_col[row_index]

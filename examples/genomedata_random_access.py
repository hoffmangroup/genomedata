#!/usr/bin/env python
from __future__ import with_statement, division

"""
Usage: ./genomedata_random_access.py GENOMEDATADIR TRACKNAME...

Expects lines on stdin of the form: chrom<whitespace>index

Prints the data value for the each trackname at each position in the format:
<value>TAB<value>TAB<value>
with a value for each trackname provided, in the order provided.
"""

import sys
import warnings
from genomedata import Genome

with Genome(sys.argv[1]) as genome:
    tracknames = sys.argv[2:]
    if tracknames:
        # Make sure all tracks can be found
        for trackname in tracknames:
            assert trackname in genome.tracknames_continuous
    else:
        print >>sys.stderr, "Using all tracks..."
        tracknames = genome.tracknames_continuous

    warnings.simplefilter("ignore")  # Ignore supercontig warnings
    for line in sys.stdin:
        chrom, index = line.strip().split()
        chromosome = genome[chrom]
        values = []
        for trackname in tracknames:
            col_index = chromosome.index_continuous(trackname)
            values.append(str(chromosome[int(index), col_index]))

        print " ".join(values)

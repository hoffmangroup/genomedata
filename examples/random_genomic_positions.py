#!/usr/bin/env python
"""
Generates random genomic positions subject to specified constraints.
Outputs N lines on stdout of the form: chrom<whitespace>index
"""

from __future__ import with_statement, division

import sys
from genomedata import Genome
from random import randint

DEFAULT_N = 1000

def rand_chrom(chroms):
    return chroms[randint(0, len(chroms) - 1)]

def rand_chrom_weighted(chrom_weights, total_weight):
    x = randint(1, total_weight)
    for chrom, weight in chrom_weights.iteritems():
        x -= weight
        if x <= 0:
            return chrom

    print >>sys.stderr, "ERROR!"
    sys.exit(1)

def rand_chromosome_position(chromosome):
    return randint(0, chromosome.end)

def rand_supercontig_position(chromosome, weight):
    x = randint(1, weight)
    for supercontig in chromosome:
        supercontig_length = supercontig.end - supercontig.start
        if x <= supercontig_length:
            return x + supercontig.start
        else:
            x -= supercontig_length

    print >>sys.stderr, "ERROR!"
    sys.exit(1)

def print_random_coordinates(genomedatadir, n=DEFAULT_N, chrom=None,
                             one_chrom=False, in_supercontigs=False):
    with Genome(genomedatadir) as genome:
        # Collect names and total lengths of chromosomes
        # If in_supercontigs, this is the total length of all supercontigs
        chrom_weights = {}
        for chromosome in genome:
            chrom_weight = 0
            if in_supercontigs:
                for supercontig in chromosome:
                    chrom_weight += supercontig.end - supercontig.start
            else:
                chrom_weight = chromosome.end - chromosome.start

            chrom_weights[chromosome.name] = chrom_weight

        total_weight = sum(chrom_weights.values())

        if chrom:
            assert chrom in chrom_weights
            one_chrom = True
            # Use specified chrom
            print >>sys.stderr, "set one_chrom"
        elif one_chrom:
            chrom = rand_chrom(chrom_weights.keys())
            print >>sys.stderr, "using only %s" % chrom

        for i in range(0, n):
            if not one_chrom:
                chrom = rand_chrom_weighted(chrom_weights, total_weight)

            if in_supercontigs:
                index = rand_supercontig_position(genome[chrom],
                                                  chrom_weights[chrom])
            else:
                index = rand_chromosome_position(genome[chrom])

            print "%s\t%d" % (chrom, index)

def parse_args(args):
    from optparse import OptionParser
    usage = "%prog [OPTIONS] GENOMEDATADIR"
    parser = OptionParser(usage=usage, description=__doc__.strip())

    parser.add_option("-n", dest="count", default=DEFAULT_N, type="int",
                      help="Number of random indices to generate"
                      " [default: %default]")
    parser.add_option("-o", "--one-chrom", action="store_true",
                      dest="one_chrom", default=False,
                      help="All random indices will be on the same chromosome")
    parser.add_option("-c", "--chrom", dest="chrom", default=None,
                      help="All random indices will be on this chromosome"
                      " (implies --one-chrom)")
    parser.add_option("-s", "--supercontig", dest="supercontig",
                      default=False, action="store_true",
                      help="All random indices will be within a supercontig")
    options, args = parser.parse_args(args)

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)

    print_random_coordinates(args[0], n=options.count, chrom=options.chrom,
                             one_chrom=options.one_chrom,
                             in_supercontigs=options.supercontig)

if __name__ == "__main__":
    sys.exit(main())

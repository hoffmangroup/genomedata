#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
_info: report specific information about a genomedata archive.
"""

# Copyright 2010, 2012, 2013 Michael M. Hoffman <mmh1@uw.edu>

import sys

from argparse import ArgumentParser

from . import Genome, __version__


def print_sizes(genome):
    # sorted by chrom size, big to small
    chrom_sizes = [(chrom.end, chrom) for chrom in genome]
    for size, chrom in reversed(sorted(chrom_sizes)):
        print(chrom, size, sep="\t")


def print_tracknames_continuous(genome):
    print(genome.tracknames_continuous, sep="\n")


def print_contigs(genome):
    for chrom in genome:
        for contig in chrom:
            print(chrom.name, contig.start, contig.end, sep="\t")


def _info(cmd, filename):

    with Genome(filename) as genome:
        if cmd in ["tracknames_continuous", "tracknames"]:
            print_tracknames_continuous(genome)
        if cmd == "contigs":
            print_contigs(genome)
        if cmd == "sizes":
            print_sizes(genome)


def parse_options(args):

    description = ("Print information about a genomedata archive.")

    parser = ArgumentParser(description=description,
                            prog='genomedata-info',
                            version=__version__)

    choices = ["tracknames", "tracknames_continuous", "contigs", "sizes"]

    parser.add_argument("command", choices=choices,
                        help='available commands')

    parser.add_argument('gdarchive', help='genomedata archive')

    args = parser.parse_args(args)

    return args


def main(argv=sys.argv[1:]):
    args = parse_options(argv)

    return _info(args.command, args.gdarchive)

if __name__ == "__main__":
    sys.exit(main())

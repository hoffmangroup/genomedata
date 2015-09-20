#!/usr/bin/env python

from __future__ import absolute_import, division, print_function
from future_builtins import ascii, filter, hex, map, oct, zip

"""
_info: DESCRIPTION
"""

__version__ = "$Revision$"

# Copyright 2010, 2012, 2013 Michael M. Hoffman <mmh1@uw.edu>

import sys

from . import Genome



def print_tracknames_continuous(genome):
    print("\n".join(genome.tracknames_continuous))


def print_contigs(genome):
    for chrom in genome:
        for contig in chrom:
            print("%s\t%s\t%s" % (chrom.name, contig.start, contig.end))


def print_sizes(genome):
    for chrom in genome:
        # add 1 because chrom sizes files are 1-based, genomedata are
        # 0-based
        print('%s\t%d' % (chrom, chrom.end + 1))


def _info(cmd, gdarchive):

    with Genome(gdarchive) as genome:
        if cmd in ["tracknames_continuous", "tracknames"]:
            print_tracknames_continuous(genome)
        if cmd == "contigs":
            print_contigs(genome)
        if cmd == 'sizes':
            print_sizes(genome)


def parse_options(args):

    from argparse import ArgumentParser
    from . import __version__

    # usage = "%(prog)s [OPTION]... CMD GENOMEDATAFILE"
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


def main(args=sys.argv[1:]):
    args = parse_options(args)

    return _info(args.command, args.gdarchive)


if __name__ == "__main__":
    sys.exit(main())

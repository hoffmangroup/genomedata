#!/usr/bin/env python


"""
_info: report specific information about a genomedata archive.
"""

__version__ = "$Revision$"

# Copyright 2010, 2012, 2013 Michael M. Hoffman <mmh1@uw.edu>

import sys

from . import Genome

from ._util import die

# XXX: add sizes command that produces tab-delimited file of sizes,
# compatible with UCSC bigWig tab-delimited specification file, for
# checking


def print_tracknames_continuous(genome):
    print("\n".join(genome.tracknames_continuous))


def print_contigs(genome):
    for chrom in genome:
        for contig in chrom:
            print("%s\t%s\t%s" % (chrom.name, contig.start, contig.end))


def _info(cmd, filename):

    with Genome(filename) as genome:
        if cmd in ["tracknames_continuous", "tracknames"]:
            print_tracknames_continuous(genome)
        if cmd == "contigs":
            print_contigs(genome)


def parse_options(args):

    from argparse import ArgumentParser
    from . import __version__

    # usage = "%(prog)s [OPTION]... CMD GENOMEDATAFILE"
    description = ("Print information about a genomedata archive.")

    parser = ArgumentParser(description=description,
                            prog='genomedata-info',
                            version=__version__)
                            
    choices = ["tracknames", "tracknames_continuous", "contigs"]

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

#!/usr/bin/env python
from __future__ import division, with_statement

"""
_info: DESCRIPTION
"""

__version__ = "$Revision$"

# Copyright 2010, 2012, 2013 Michael M. Hoffman <mmh1@uw.edu>

import sys

from . import Genome

# XXX: add sizes command that produces tab-delimited file of sizes,
# compatible with UCSC bigWig tab-delimited specification file, for
# checking


def die(msg="Unexpected error."):
    print >>sys.stderr, msg
    sys.exit(1)


def print_tracknames_continuous(genome):
    print "\n".join(genome.tracknames_continuous)


def print_contigs(genome):
    for chrom in genome:
        for contig in chrom:
            print "%s\t%s\t%s" % (chrom.name, contig.start, contig.end)


def _info(cmd, filename):
    choices = ["tracknames", "tracknames_continuous", "contigs"]
    if cmd not in frozenset(choices):
        die("CMD must be one of: %s" % ",".join(choices))

    with Genome(filename) as genome:
        if cmd in ["tracknames_continuous", "tracknames"]:
            print_tracknames_continuous(genome)
        if cmd == "contigs":
            print_contigs(genome)


def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... CMD ARCHIVE"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    options, args = parser.parse_args(args)

    if not len(args) == 2:
        parser.error("incorrect number of arguments")

    return options, args


def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    return _info(*args)

if __name__ == "__main__":
    sys.exit(main())

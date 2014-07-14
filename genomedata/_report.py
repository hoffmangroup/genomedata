#!/usr/bin/env python
from __future__ import division, with_statement

"""
report: report some summary statistics on a genomedata archive that
already has save_metadata() run
"""

__version__ = "$Revision$"

# Copyright 2009-2014 Michael M. Hoffman <mhoffman@uhnresearch.ca>

import sys

from genomedata import Genome
from tabdelim import ListWriter

def report(genomedata):
    writer = ListWriter()

    with Genome(genomedata) as genome:
        writer.writerow(["measurement"] + genome.tracknames_continuous)
        writer.writerow(["mean"] + list(genome.means))
        writer.writerow(["var"] + list(genome.vars))

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]..."
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    options, args = parser.parse_args(args)

    if not len(args) == 1:
        parser.print_usage()
        sys.exit(1)

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    return report(*args)

if __name__ == "__main__":
    sys.exit(main())

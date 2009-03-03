#!/usr/bin/env python
from __future__ import division, with_statement

__version__ = "$Revision$"

# Copyright 2008-2009 Michael M. Hoffman <mmh1@washington.edu>

from contextlib import closing
from gzip import open as _gzip_open
import sys

from numpy import array, empty
from tables import Filters, NoSuchNodeError

FILTERS_GZIP = Filters(complevel=1)

EXT_GZ = "gz"

class LightIterator(object):
    def __init__(self, handle):
        self._handle = handle
        self._defline = None

    def __iter__(self):
        return self

    def next(self):
        lines = []
        defline_old = self._defline

        for line in self._handle:
            if not line:
                if not defline_old and not lines:
                    raise StopIteration
                if defline_old:
                    self._defline = None
                    break
            elif line.startswith(">"):
                self._defline = line[1:].rstrip()
                if defline_old or lines:
                    break
                else:
                    defline_old = self._defline
            else:
                lines.append(line.rstrip())

        if not lines:
            raise StopIteration

        return defline_old, ''.join(lines)

# XXX: suggest as default
def fill_array(scalar, shape, dtype=None, *args, **kwargs):
    if dtype is None:
        dtype = array(scalar).dtype

    res = empty(shape, dtype, *args, **kwargs)
    res.fill(scalar)

    return res

# XXX: suggest as default
def gzip_open(*args, **kwargs):
    return closing(_gzip_open(*args, **kwargs))

# XXX: replace with genomedata.Chromosome.tracknames_continuous
def get_tracknames(chromosome):
    return chromosome.root._v_attrs.tracknames.tolist()

def init_num_obs(num_obs, continuous):
    curr_num_obs = continuous.shape[1]
    assert num_obs is None or num_obs == curr_num_obs

    return curr_num_obs

def new_extrema(func, data, extrema):
    curr_extrema = func(data, 0)

    return func([extrema, curr_extrema], 0)

# XXX: replace with iter(genomedata.Chromosome)
def walk_supercontigs(h5file):
    root = h5file.root

    for supercontig in h5file.walkGroups():
        if supercontig == root:
            continue

        yield supercontig

# XXX: replace with genomedata.Chromosome.itercontinuous
def walk_continuous_supercontigs(h5file):
    for supercontig in walk_supercontigs(h5file):
        try:
            yield supercontig, supercontig.continuous
        except NoSuchNodeError:
            continue

def main(args=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())

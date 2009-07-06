#!/usr/bin/env python
from __future__ import division

__version__ = "$Revision$"

# Copyright 2009 Michael M. Hoffman <mmh1@washington.edu>

from functools import partial
from os import extsep
import sys

from numpy import add, amin, amax
from path import path
from tables import openFile, NoSuchNodeError

FORMAT_VERSION = 0

EXT = "h5" # XXX: change to .genomedata
SUFFIX = extsep + EXT

class InactiveDict(set):
    """
    fake dict that can't be added to
    """
    def __setitem__(self, key, value):
        return

class Genome(object):
    """
    implemented via a file system directory

    if you use this as a context manager, it will keep track of open
    Chromosomes and close them for you later when the context is left
    """
    def __init__(self, dirname):
        self.dirpath = path(dirname)

        # used when the Genome instance is not used as a context
        # manager. replaced by __enter__()
        self.open_chromosomes = InactiveDict()

    def __iter__(self):
        # sorted so that the order is always the same
        for filepath in sorted(self.dirpath.files("*" + SUFFIX)):

            # pass through __getitem__() to allow memoization
            yield self[filepath.namebase]

    def __getitem__(self, name):
        try:
            # memoization
            return self.open_chromosomes[name]
        except KeyError:
            pass

        res = Chromosome(self.dirpath / (name + SUFFIX))

        self.open_chromosomes[name] = res
        return res

    def __enter__(self):
        self.open_chromosomes = {}

        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        for chromosome in self.open_chromosomes:
            chromosome.close()

    def _accum_extrema(self, name, accumulator):
        self.tracknames_continuous # for assertion check

        res = None

        for chromosome in self:
            new_extrema = getattr(chromosome, name)

            if res is None:
                res = new_extrema
            else:
                res = accumulator([res, new_extrema])

        return res

    @property
    def tracknames_continuous(self):
        res = None

        # check that all chromosomes have the same tracknames_continuous
        for chromosome in self:
            if res is None:
                res = chromosome.tracknames_continuous
            else:
                assert res == chromosome.tracknames_continuous

        return res

    # XXX: should memoize these with an off-the-shelf decorator
    @property
    def mins(self):
        return self._accum_extrema("mins", partial(amin, axis=0))

    @property
    def maxs(self):
        return self._accum_extrema("maxs", partial(amax, axis=0))

    @property
    def sums(self):
        return self._accum_extrema("sums", add.reduce)

    @property
    def sums_squares(self):
        return self._accum_extrema("sums", add.reduce)

    @property
    def num_datapoints(self):
        return self._accum_extrema("sums", add.reduce)

class Chromosome(object):
    """
    implemented via an HDF5 File
    """
    # XXX: I need to handle the dirty case better, to allow "+" in mode
    def __init__(self, filename, mode="r", *args, **kwargs):
        h5file = openFile(filename, mode, *args, **kwargs)
        attrs = h5file.root._v_attrs

        # set or check file format version and dirty flag
        if "w" in mode:
            attrs.genomedata_format_version = FORMAT_VERSION
            attrs.dirty = True
        else:
            # XXX: the first version did not necessarily include this
            # attribute; after everything does this exception handling
            # should be removed
            try:
                assert attrs.genomedata_format_version == FORMAT_VERSION
            except AttributeError:
                if FORMAT_VERSION != 0:
                    raise

            assert not attrs.dirty

        self.filename = filename
        self.h5file = h5file

    def __iter__(self):
        h5file = self.h5file
        root = h5file.root

        for group in h5file.walkGroups():
            if group == root:
                continue

            yield Supercontig(group)

    def __getitem__(self, key):
        """
        returns the supercontig that contains this range if possible
        """
        if isinstance(key, slice):
            start = key.start
            end = key.stop
        else:
            # number
            start = key
            end = key + 1

        for supercontig in self:
            if start >= supercontig.start and end <= supercontig.end:
                return supercontig

        raise IndexError("no supercontigs matched %r" % key)

    def itercontinuous(self):
        for supercontig in self:
            try:
                yield supercontig, supercontig.continuous
            except NoSuchNodeError:
                continue

    def index_continuous(self, trackname):
        return self.tracknames_continuous.index(trackname)

    def close(self):
        return self.h5file.close()

    @property
    def name(self):
        return path(self.filename).name.rpartition(SUFFIX)[0]

    @property
    def attrs(self):
        return self.h5file.root._v_attrs

    @property
    def tracknames_continuous(self):
        return self.attrs.tracknames.tolist()

    @property
    def mins(self):
        return self.attrs.mins

    @property
    def maxs(self):
        return self.attrs.maxs

    @property
    def sums(self):
        return self.attrs.sums

    @property
    def sums_squares(self):
        return self.attrs.sums_squares

    @property
    def num_datapoints(self):
        return self.attrs.num_datapoints

    @property
    def end(self):
        return max(supercontig.end for supercontig in self)

class Supercontig(object):
    """
    implemented via a HDF5 Group
    """
    def __init__(self, h5group):
        self.h5group = h5group

    @property
    def continuous(self):
        return self.h5group.continuous

    @property
    def attrs(self):
        return self.h5group._v_attrs

    @property
    def name(self):
        return self.h5group._v_name

    @property
    def seq(self):
        return self.h5group.seq

    @property
    def start(self):
        return self.attrs.start

    @property
    def end(self):
        return self.attrs.end

    def project(self, pos):
        """
        project chromsomal coordinates to supercontig coordinates
        """
        return pos - self.start

def main(args=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())

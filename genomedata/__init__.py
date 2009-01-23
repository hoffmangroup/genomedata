#!/usr/bin/env python
from __future__ import division

__version__ = "$Revision$"

# Copyright 2009 Michael M. Hoffman <mmh1@washington.edu>

from os import extsep
import sys

from path import path
from tables import openFile, NoSuchNodeError

FORMAT_VERSION = 0

EXT = "h5" # XXX: change to .genomedata
SUFFIX = extsep + EXT

class Genome(object):
    """
    implemented via a file system directory
    """
    def __init__(self, dirname):
        self.dirpath = path(dirname)
        self.open_chromosomes = None

    def __iter__(self):
        for filepath in self.dirpath.files():
            chromosome = Chromosome(filepath)
            if self.open_chromosomes is not None:
                self.open_chromosomes.add(chromosome)

            yield chromosome

    def __getitem__(self):
        raise NotImplementedError

    def __enter__(self):
        self.open_chromosomes = set()

    def __exit__(self, exc_type, exc_value, exc_tb):
        for chromosome in self.open_chromosomes:
            chromosome.close()

class Chromosome(object):
    """
    implemented via an HDF5 File
    """
    def __init__(self, filename, mode="r", *args, **kwargs):
        h5file = openFile(filename, mode, *args, **kwargs)
        attrs = h5file.root._v_attrs

        # set or check file format version and dirty flag
        if mode == "w":
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
    def start(self):
        return self.attrs.start

    @property
    def end(self):
        return self.attrs.end

def main(args=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())

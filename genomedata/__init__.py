#!/usr/bin/env python
from __future__ import division

__version__ = "$Revision$"

# Copyright 2009 Michael M. Hoffman <mmh1@washington.edu>

from os import extsep
import sys

from path import path
from tables import openFile

FORMAT_VERSION = 0

EXT = "genomedata"
SUFFIX = extsep + EXT

class Genome(object):
    """
    implemented via a file system directory
    """
    def __init__(self, dirname):
        self.dirpath = path(dirname)

    @staticmethod
    def _filepath2name(filepath):
        return filepath.name.rpartition(SUFFIX)[0]

    def __iter__(self):
        """
        iterate through the chromosome names
        """
        for filepath in self.dirpath.files():
            yield self._filepath2name(filepath)

    def iterchromosomes(self):
        for filepath in self.dirpath.files():
            yield Chromosome(filepath)

    def iteritems(self):
        for filepath in self.dirpath.files():
            yield self._filepath2name(filepath), Chromosome(filepath)

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

        self.h5file = h5file

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

class Supercontig(object):
    """
    implemented via a HDF5 Group
    """
    pass

def main(args=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())

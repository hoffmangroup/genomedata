#!/usr/bin/env python
"""
Genomedata is a module to store and access large-scale functional
genomics data in a format which is both space-efficient and allows
efficient random-access.

Under the surface, genomedata is implemented as a collection of HDF5 files,
but genomedata provides a transparent interface to interact with your
underlying data without having to worry about the mess of repeatedly parsing
large data files or having to keep them in memory for random access.

.. Copyright 2009 Michael M. Hoffman <mmh1@washington.edu>

"""

from __future__ import division, with_statement

__version__ = "$Revision$"


from functools import partial
from os import extsep
import sys

from numpy import add, amin, amax, square
from path import path
from tables import openFile, NoSuchNodeError

FORMAT_VERSION = 0

EXT = "genomedata"
SUFFIX = extsep + EXT

class _InactiveDict(dict):
    """A fake dict that can't be added to."""
    def __setitem__(self, key, value):
        return

class Genome(object):
    """The root level of the genomedata object hierarchy.

    Implemented via a file system directory
    If you use this as a context manager, it will keep track of open
    Chromosomes and close them for you later when the context is left.
    
    """
    def __init__(self, dirname):
        """Create a Genome object from the genomedata objects in the directory.

        :param dirname: directory containing any chomosome files to include.
        :type dirname: string
        
        """
        self.dirpath = path(dirname)

        # used when the Genome instance is not used as a context
        # manager. replaced by __enter__()
        self.open_chromosomes = _InactiveDict()

        # a kind of refcounting for context managers
        self._context_count = 0

    def __iter__(self):
        """Return next chromosome, in sorted order, with memoization"""
        # sorted so that the order is always the same
        for filepath in sorted(self.dirpath.files("*" + SUFFIX)):

            # pass through __getitem__() to allow memoization
            yield self[filepath.namebase]

    def __getitem__(self, name):
        """Return a reference to a chromosome of the given name.

        :param name: name of the chromosome file (e.g. "chr1" if chr1.genomedata
                     is a file in the genomedata directory)
        :type name: string
        :rtype: Chromosome_
        
        """
        try:
            # memoization
            return self.open_chromosomes[name]
        except KeyError:
            pass

        res = Chromosome(self.dirpath / (name + SUFFIX))

        self.open_chromosomes[name] = res
        return res

    def __enter__(self):
        if self._context_count == 0:
            self.open_chromosomes = {}

        self._context_count += 1

        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        if self._context_count == 1:
            for name, chromosome in self.open_chromosomes.iteritems():
                chromosome.close()

            self.open_chromosomes = _InactiveDict()

        self._context_count -= 1

    def _accum_extrema(self, name, accumulator):
        res = None

        with self:
            self.tracknames_continuous # for assertion check

            for chromosome in self:
                new_extrema = getattr(chromosome, name)

                if res is None:
                    res = new_extrema
                else:
                    res = accumulator([res, new_extrema])

        return res

    @property
    def tracknames_continuous(self):
        """Return a list of the names of all data tracks stored."""
        res = None

        # check that all chromosomes have the same tracknames_continuous
        with self:
            for chromosome in self:
                if res is None:
                    res = chromosome.tracknames_continuous
                else:
                    assert res == chromosome.tracknames_continuous

        return res

    # XXX: should memoize these with an off-the-shelf decorator
    @property
    def mins(self):
        """Return the minimum value for each track.

        :rtype: numpy.array
        
        """
        return self._accum_extrema("mins", partial(amin, axis=0))

    @property
    def maxs(self):
        """Return a vector of the maximum value for each track.

        :rtype: numpy.array
        
        """
        return self._accum_extrema("maxs", partial(amax, axis=0))

    @property
    def sums(self):
        """Return a vector of the sum of the values for each track.

        :rtype: numpy.array
        
        """
        return self._accum_extrema("sums", add.reduce)

    @property
    def sums_squares(self):
        """Return a vector of the sum of squared values for each track's data.

        :rtype: numpy.array
        
        """
        return self._accum_extrema("sums_squares", add.reduce)

    @property
    def num_datapoints(self):
        """Return the number of datapoints in each track.

        :returns: a vector with an entry for each track.
        :rtype: numpy.array
        
        """
        return self._accum_extrema("num_datapoints", add.reduce)

    @property
    def means(self):
        """Return a vector of the mean value of each track.

        :rtype: numpy.array
        
        """
        with self:
            return self.sums / self.num_datapoints

    @property
    def vars(self):
        """Return a vector of the variance in the data for each track.

        :rtype: numpy.array
        
        """
        # this is an unstable way of calculating the variance,
        # but it should be good enough
        # Numerical Recipes in C, Eqn 14.1.7
        # XXX: best would be to switch to the pairwise parallel method
        # (see Wikipedia)
        with self:
            return (self.sums_squares / self.num_datapoints) - square(self.means)

class Chromosome(object):
    """The genomedata object corresponding to data for a given chromosome.
    
    Implemented via an HDF5 File
    
    """
    # XXX: I need to handle the dirty case better, to allow "+" in mode
    def __init__(self, filename, mode="r", *args, **kwargs):
        """
        :param filename: name of the chromosome file in the genomedata directory
        :param mode: mode of interaction with the chromosome file, with
                     ``r``: read, ``w``: write
        :type mode: string
        :param \*args: args passed on to openFile
        :param \*\*kwargs: keyword args passed on to openFile
        
        """
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
        """Return next supercontig in chromosome."""
        h5file = self.h5file
        root = h5file.root

        for group in h5file.walkGroups():
            if group == root:
                continue

            yield Supercontig(group)

    def __getitem__(self, key):
        """Return the supercontig that contains this range if possible.

        :param key: index or range of indices to find containing supercontig for
        :type key: slice or integer
        :rtype: Supercontig_
        
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
        """Return a generator over all superconig, continuous pairs."""
        for supercontig in self:
            try:
                yield supercontig, supercontig.continuous
            except NoSuchNodeError:
                continue

    def index_continuous(self, trackname):
        """Return the column index of the trackname in the continuous data.

        :param trackname: name of data track
        :type trackname: string
        :rtype: integer
        
        """
        return self.tracknames_continuous.index(trackname)

    def close(self):
        """Close the current chromosome file."""
        return self.h5file.close()

    @property
    def name(self):
        """Return the name of this chromosome."""
        return path(self.filename).name.rpartition(SUFFIX)[0]

    @property
    def attrs(self):
        """Return the attributes of this chromosome."""
        return self.h5file.root._v_attrs

    @property
    def tracknames_continuous(self):
        """Return a list of the data track names in this Chromosome."""
        return self.attrs.tracknames.tolist()

    @property
    def mins(self):
        """See Genome.mins."""
        return self.attrs.mins

    @property
    def maxs(self):
        """See Genome.maxs."""
        return self.attrs.maxs

    @property
    def sums(self):
        """See Genome.sums."""
        return self.attrs.sums

    @property
    def sums_squares(self):
        """See Genome.sums_squares."""
        return self.attrs.sums_squares

    @property
    def num_datapoints(self):
        """See Genome.num_datapoints."""
        return self.attrs.num_datapoints

    @property
    def end(self):
        """Return the position of the last base in this chromosome."""
        return max(supercontig.end for supercontig in self)

class Supercontig(object):
    """A container for a segment of data in one chromosome.
    
    Implemented via a HDF5 Group
    
    """
    def __init__(self, h5group):
        """
        :param h5group: group containing the Supercontig data
        :type h5group: HDF5 Group
        
        """
        self.h5group = h5group

    @property
    def continuous(self):
        """Return the underlying continuous data in this supercontig.

        :rtype: numpy.array
        
        """
        return self.h5group.continuous

    @property
    def attrs(self):
        """Return the attributes of this supercontig."""
        return self.h5group._v_attrs

    @property
    def name(self):
        """Return the name of this supercontig."""
        return self.h5group._v_name

    @property
    def seq(self):
        """Return the genomic sequence of this supercontig."""
        return self.h5group.seq

    @property
    def start(self):
        """Return the start position of this supercontig."""
        return self.attrs.start

    @property
    def end(self):
        """Return the end position of this supercontig."""
        return self.attrs.end

    def project(self, pos):
        """Project chromsomal coordinates to supercontig coordinates.

        :param pos: chromosome coordinate
        :type pos: integer
        :rtype: integer
        
        """
        return pos - self.start

def main(args=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())

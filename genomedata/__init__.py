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

from numpy import add, amin, amax, empty, NAN, square
from path import path
from tables import openFile, NoSuchNodeError
from warnings import warn

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

        :param name: name of the chromosome file (e.g. "chr1" if
                     chr1.genomedata is a file in the genomedata directory)
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

    @property
    def num_tracks_continuous(self):
        """Returns the number of continuous data tracks."""
        res = None

        # check that all chromosomes have the same tracknames_continuous
        with self:
            for chromosome in self:
                if res is None:
                    res = chromosome.num_tracks_continuous
                else:
                    assert res == chromosome.num_tracks_continuous

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
        self._seq = _ChromosomeSeqSlice(self)
        self._supercontigs = _Supercontigs(self)

    def __iter__(self):
        """Return next supercontig in chromosome."""
        h5file = self.h5file
        root = h5file.root

        for group in h5file.walkGroups():
            if group == root:
                continue

            yield Supercontig(group)

    def __getitem__(self, key):
        """Return the contiguous data corresponding to this bp slice

        If slice is taken over a supercontig boundary, missing data
        is filled in with NaN's automatically and a warning is printed.

        :param key: index or range of indices to get continuous data for
        :type key: slice or integer
        :rtype: numpy.array

        """
        if isinstance(key, tuple):
            key, cols = key
            if isinstance(cols, slice):
                if cols.start is None:
                    cols.start = 0
                if cols.stop is None:
                    cols.stop = self.num_tracks_continuous
                num_cols = cols.stop - cols.start
            else:
                num_cols = 1
            print str(key), str(cols), str(num_cols)
        else:
            cols = slice(None)
            num_cols = self.num_tracks_continuous
            print str(cols), str(num_cols)

        supercontigs = self.supercontigs[key]
        if len(supercontigs) == 0:
            warn("slice of chromosome data does not overlap any supercontig"
                 " (filling with 'NaN')")
        elif len(supercontigs) > 1:
            warn("slice of chromosome data spans more than one supercontig"
                 " (filling gaps with 'NaN')")

        start, end = _key_to_chrom_range(key, self)
        data = empty((end - start, num_cols),
                     dtype=self._continuous_dtype)
        data.fill(NAN)
        for supercontig in supercontigs:
            chr_start = max(start, supercontig.start)
            chr_end = min(end, supercontig.end)
            data[chr_start - start:chr_end - start, cols] = \
                supercontig.continuous[supercontig.project(chr_start):
                                           supercontig.project(chr_end),
                                       cols]

        return data

    def itercontinuous(self):
        """Return a generator over all supercontig, continuous pairs."""
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
    def _continuous_dtype(self):
        for supercontig in self:
            return supercontig._continuous_dtype

    @property
    def _seq_dtype(self):
        for supercontig in self:
            return supercontig._seq_dtype

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
    def num_tracks_continuous(self):
        """Returns the number of tracks in this chromosome"""
        return len(self.tracknames_continuous)

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
    def start(self):
        """Return the position of the first base in this chromosome."""
        return min(supercontig.start for supercontig in self)

    @property
    def end(self):
        """Return the position of the last base in this chromosome."""
        return max(supercontig.end for supercontig in self)

    @property
    def seq(self):
        return self._seq

    @property
    def supercontigs(self):
        """Return the supercontig that contains this range if possible.

        Get items with a slice or simple index
        :rtype: Supercontig_

        """
        return self._supercontigs

class _ChromosomeSeqSlice(object):
    def __init__(self, chromosome):
        self._chromosome = chromosome

    def __getitem__(self, key):
        """Get the underlying sequence that corresponds to this index (range).

        Insert "N"s if the index range spans no/multiple supercontigs.
        """

        supercontigs = self._chromosome.supercontigs[key]

        if len(supercontigs) == 0:
            warn("slice of chromosome sequence does not overlap any"
                 " supercontig (filling with 'N')")
        elif len(supercontigs) > 1:
            warn("slice of chromosome sequence spans more than one supercontig"
                 " (filling gaps with 'NaN')")

        start, end = _key_to_chrom_range(key, self._chromosome)
        seq = empty((end - start,), dtype=self._chromosome._seq_dtype)
        seq.fill(ord("N"))
        for supercontig in supercontigs:
            chr_start = max(start, supercontig.start)
            chr_end = min(end, supercontig.end)
            seq[chr_start - start:chr_end - start] = \
                supercontig.seq[supercontig.project(chr_start):
                                    supercontig.project(chr_end)]

        return seq

class _Supercontigs(object):
    def __init__(self, chromosome):
        self._chromosome = chromosome

    def __getitem__(self, key):
        """Get the supercontigs that contains this index range, if possible."""
        start, end = _key_to_chrom_range(key, self._chromosome)

        supercontigs = []
        for supercontig in self._chromosome:
            if start < supercontig.end and end > supercontig.start:
                supercontigs.append(supercontig)
                if start >= supercontig.start and end <= supercontig.end:
                    # Key entirely within one supercontig, so we're done
                    break
            # XXX: would be nice if we could count on supercontig ordering

        return supercontigs

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

    def project(self, pos):
        """Project chromsomal coordinates to supercontig coordinates.

        :param pos: chromosome coordinate
        :type pos: integer
        :rtype: integer

        """
        return pos - self.start

    @property
    def _seq_dtype(self):
        return self.seq.atom.dtype

    @property
    def _continuous_dtype(self):
        return self.continuous.atom.dtype

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

def _key_to_range(key):
    """Convert a slice or integer key to (start, stop) tuple."""
    if isinstance(key, slice):
        return (key.start, key.stop)
    else:
        start = int(key)
        return (start, start + 1)

def _key_to_chrom_range(key, chrom):
    """Key to range with enforcement within chromosome bounds."""
    start, end = _key_to_range(key)
    if start is None:
        start = chrom.start
    if end is None:
        end = chrom.end
    # Check chrom bounds
    if start < 0 or end > chrom.end:
        raise IndexError("Chromosome sequence index: %r out of chromosome"
                         " range: [0, %d)" % (key, chrom.end))
    return start, end


def main(args=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())

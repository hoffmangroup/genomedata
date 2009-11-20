#!/usr/bin/env python
"""
Genomedata is a module to store and access large-scale functional
genomics data in a format which is both space-efficient and allows
efficient random-access.

Under the surface, genomedata is implemented as a collection of HDF5 files,
but genomedata provides a transparent interface to interact with your
underlying data without having to worry about the mess of repeatedly parsing
large data files or having to keep them in memory for random access.

Copyright 2009 Michael M. Hoffman <mmh1@washington.edu>

"""

from __future__ import division, with_statement

__version__ = "$Revision$"

import sys

from functools import partial
from numpy import add, amin, amax, array, empty, NAN, square
from os import extsep
from path import path
from tables import openFile, NoSuchNodeError
from warnings import warn

FORMAT_VERSION = 0
DEFAULT_SEQ_DTYPE = "uint8"
DEFAULT_CONTINUOUS_DTYPE = "float32"

EXT = "genomedata"
SUFFIX = extsep + EXT

class _InactiveDict(dict):
    """A fake dict that can't be added to."""
    def __setitem__(self, key, value):
        return

class Genome(object):
    """The root level of the genomedata object hierarchy.

    Implemented via a file system directory.
    If you use this as a context manager, it will keep track of open
    Chromosomes and close them for you later when the context is left::

      with Genome("/path/to/genomedata") as genome:
        chromosome = genome["chr1"]
        [...]

    If not used as a context manager, you are responsible for closing
    chromosomes when you are done with them:

    >>> genome = Genome("/path/to/genomedata")
    >>> chromosome = genome["chr1"]
    [...]
    >>> chromosome.close()

    """
    def __init__(self, dirname):
        """Create a Genome object from the genomedata objects in the directory.

        :param dirname: directory containing any chomosome files to include
                        (usually just the genomedata directory).
        :type dirname: string

        Example:

        >>> genome = Genome("./genomedata.ctcf.pol2b/")
        >>> genome
        Genome("./genomedata.ctcf.pol2b/")

        """
        self.dirpath = path(dirname)

        # used when the Genome instance is not used as a context
        # manager. replaced by __enter__()
        self.open_chromosomes = _InactiveDict()

        # a kind of refcounting for context managers
        self._context_count = 0

    def __iter__(self):
        """Return next chromosome, in sorted order, with memoization.

        Example::

          for chromosome in genome:
            print chromosome.name
            for supercontig, continuous in chromosome.itercontinuous():
              [...]

        """
        # sorted so that the order is always the same
        for filepath in sorted(self.dirpath.files("*" + SUFFIX)):

            # pass through __getitem__() to allow memoization
            yield self[filepath.namebase]

    def __getitem__(self, name):
        """Return a reference to a chromosome of the given name.

        :param name: name of the chromosome file (e.g. "chr1" if
                     chr1.genomedata is a file in the genomedata directory)
        :type name: string
        :returns: :class:`Chromosome`

        Example:

        >>> genome["chrX"]
        Chromosome('/path/to/genomedata/chrX.genomedata')
        >>> genome["chrZ"]
        KeyError: 'Could not find chromosome: chrZ'

        """
        try:
            # memoization
            return self.open_chromosomes[name]
        except KeyError:
            pass

        try:
            res = Chromosome(self.dirpath / (name + SUFFIX))
        except IOError:
            raise KeyError("Could not find chromosome: %s" % name)

        self.open_chromosomes[name] = res
        return res

    def __enter__(self):
        if self._context_count == 0:
            self.open_chromosomes = {}

        self._context_count += 1

        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        # XXX: this and __enter__ have potential race conditions, if
        # _context_count is changed simultaneously by different
        # threads. should be synchronized

        if self._context_count == 1:
            for name, chromosome in self.open_chromosomes.iteritems():
                chromosome.close()

            self.open_chromosomes = _InactiveDict()

        self._context_count -= 1

    def __repr__(self):
        return "Genome('%s')" % self.dirpath

    def __str__(self):
        return repr(self)

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

        :returns: numpy.array

        """
        return self._accum_extrema("mins", partial(amin, axis=0))

    @property
    def maxs(self):
        """Return a vector of the maximum value for each track.

        :returns: numpy.array

        """
        return self._accum_extrema("maxs", partial(amax, axis=0))

    @property
    def sums(self):
        """Return a vector of the sum of the values for each track.

        :returns: numpy.array

        """
        return self._accum_extrema("sums", add.reduce)

    @property
    def sums_squares(self):
        """Return a vector of the sum of squared values for each track's data.

        :returns: numpy.array

        """
        return self._accum_extrema("sums_squares", add.reduce)

    @property
    def num_datapoints(self):
        """Return the number of datapoints in each track.

        :returns: a numpy.array vector with an entry for each track.

        """
        return self._accum_extrema("num_datapoints", add.reduce)

    @property
    def means(self):
        """Return a vector of the mean value of each track.

        :returns: numpy.array

        """
        with self:
            return self.sums / self.num_datapoints

    @property
    def vars(self):
        """Return a vector of the variance in the data for each track.

        :returns: numpy.array

        """
        # this is an unstable way of calculating the variance,
        # but it should be good enough
        # Numerical Recipes in C, Eqn 14.1.7
        # XXX: best would be to switch to the pairwise parallel method
        # (see Wikipedia)
        with self:
            return (self.sums_squares / self.num_datapoints) - \
                square(self.means)

class Chromosome(object):
    """The genomedata object corresponding to data for a given chromosome.

    Implemented via an HDF5 File

    Usually created by keying into a Genome object with the name of a
    chromosome, as in:

    >>> with Genome("/path/to/genomedata") as genome:
    ...     chromosome = genome["chrX"]
    ...     chromosome
    ...
    Chromosome('/path/to/genomedata/chrX.genomedata')

    """
    default_mode = "r"
    def __init__(self, filename, mode=default_mode, *args, **kwargs):
        """
        :param filename: name of the chromosome file in the
                         genomedata directory

        :param mode: mode of interaction with the chromosome file,
                     with ``r``: read, ``w``: write

        :type mode: string
        :param \*args: args passed on to openFile
        :param \*\*kwargs: keyword args passed on to openFile

        """

        # disabled possibilities:
        # , ``a``: append, ``r+``: append but force file to exist
        # already. See documentation for tables.openFile().

        h5file = openFile(filename, mode, *args, **kwargs)
        attrs = h5file.root._v_attrs

        # set or check file format version and dirty flag
        # XXX: need to handle the dirty case better, to allow "+" in mode
        assert mode in set("rw") # others not allowed yet

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
        self.mode = mode
        self.h5file = h5file
        self._seq = _ChromosomeSeqSlice(self)
        self._supercontigs = _Supercontigs(self)

    def __iter__(self):
        """Return next supercontig in chromosome.

        Seldom used in favor of the more direct:
        :meth:`Chromosome.itercontinuous`

        Example:

        >>> for supercontig in chromosome:
        ...     supercontig  # calls repr()
        ...
        <Supercontig('supercontig_0', 0:66115833)>
        <Supercontig('supercontig_1', 66375833:90587544)>
        <Supercontig('supercontig_2', 94987544:199501827)>

        """
        h5file = self.h5file
        root = h5file.root

        for group in h5file.walkGroups():
            if group == root:
                continue

            yield Supercontig(group)

    def __getitem__(self, key):
        """Return the continuous data corresponding to this bp slice

        :param key: key must index or slice bases, but can also index, slice,
                    or directly specify (string or list of strings) the data
                    tracks.

        :type key: <base_key>[, <track_key>]
        :returns: numpy.array

        If slice is taken over or outside a supercontig boundary,
        missing data is filled in with NaN's automatically and a
        warning is printed.

        Typical use:

        >>> chromosome = genome["chr4"]
        >>> chromosome[0:5]  # Get all data for the first five bases of chr4
        >>> chromosome[0, 0:2]  # Get data for first two tracks at chr4:0
        >>> chromosome[100, "ctcf"]  # Get "ctcf" track value at chr4:100

        """
        # XXX: Allow variable/negative steps, negative starts/stops, etc.

        # XXX: The problem with missing start or end indices is that
        # 1) it is unclear if 0 or self.start should be used for the
        #    start, and (preference for 0)
        # 2) the full length of the chromosome is not known, so self.end
        #    would need to be used, which might not be what is wanted.

        # Sanitize the input
        if isinstance(key, tuple):
            base_key, track_key = key
        else:
            base_key = key
            track_key = slice(None)  # All tracks

        # Treat direct indexing differently (just like numpy)
        base_index = False
        track_index = False

        if isinstance(base_key, int):
            base_index = True

        base_start, base_stop = _key_to_tuple(base_key)
        base_key = slice(base_start, base_stop)

        # First convert track_key toward slice
        if isinstance(track_key, basestring):
            track_key = self.index_continuous(track_key)

        if isinstance(track_key, int):
            track_key = slice(track_key, track_key + 1, 1)
            track_index = True

        if isinstance(track_key, slice):
            # Fix indices to number of tracks
            track_key = slice(*track_key.indices(self.num_tracks_continuous))
        else:
            raise TypeError("Unrecognized track indexing method: %s" %
                            track_key)

        nrows = base_key.stop - base_key.start
        ncols = len(xrange(track_key.start, track_key.stop, track_key.step))

        print track_key, ncols
        # Handle degenerate case
        dtype = self._continuous_dtype
        if nrows < 1 or ncols < 1:
            # Return empty array (matches numpy behavior)
            return array((), dtype=dtype)

        # At this point, base_key and track_key are guaranteed to be slices
        # with both start and end >= 0

        # Lookup appropriate data
        supercontigs = self.supercontigs[base_key]
        if len(supercontigs) == 0:
            warn("slice of chromosome data does not overlap any supercontig"
                 " (filling with 'NaN')")
        elif len(supercontigs) > 1:
            warn("slice of chromosome data spans more than one supercontig"
                 " (filling gaps with 'NaN')")

        data = empty((nrows, ncols), dtype=dtype)
        data.fill(NAN)

        for supercontig in supercontigs:
            assert base_key.start < supercontig.end and \
                base_key.stop > supercontig.start
            chr_start = max(base_key.start, supercontig.start)
            chr_end = min(base_key.stop, supercontig.end)
            data_slice = slice(chr_start - base_key.start,
                               chr_end - base_key.start)
            supercontig_slice = slice(supercontig.project(chr_start),
                                      supercontig.project(chr_end))
            # track_key must be a splice
            try:
                data[data_slice, :] = supercontig.continuous[supercontig_slice,
                                                          track_key]
            except NoSuchNodeError:
                # Allow the supercontig to not have a continuous dataset
                pass

        # Make output shape appropriate for indexing method (like numpy)
        if track_index:
            data = data[:, 0]
        if base_index:
            data = data[0]
        return data

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        if self.mode == self.default_mode:
            return "Chromosome('%s')" % (self.filename)
        else:
            return "Chromosome('%s', '%s')" % (self.filename, self.mode)

    def itercontinuous(self):
        """Return a generator over all supercontig, continuous pairs.

        This is the best way to efficiently iterate over the data since
        all specified data is in supercontigs.

        """
        for supercontig in self:
            try:
                yield supercontig, supercontig.continuous
            except NoSuchNodeError:
                continue

    def index_continuous(self, trackname):
        """Return the column index of the trackname in the continuous data.

        :param trackname: name of data track
        :type trackname: string
        :returns: integer

        This is used for efficient indexing into continuous data:

        >>> chromosome = genome["chr3"]
        >>> col_index = chromosome.index_continuous("sample_track")
        >>> data = chromosome[100:150, col_index]

        although for typical use, the track can be indexed directly:

        >>> data = chromosome[100:150, "sample_track"]

        """
        try:
            return self.tracknames_continuous.index(trackname)
        except ValueError:
            raise KeyError("Could not find continuous track: %s" % trackname)

    def close(self):
        """Close the current chromosome file.

        This only needs to be called when genomedata is not being used as
        a context manager. Using genomedata as a context manager makes
        life easy by memoizing chromosome access and guaranteeing the
        proper cleanup. See :class:`Genome`.

        """
        return self.h5file.close()

    @property
    def _continuous_dtype(self):
        for supercontig, continuous in self.itercontinuous():
            return supercontig._continuous_dtype
        return DEFAULT_CONTINUOUS_DTYPE

    @property
    def _seq_dtype(self):
        for supercontig in self:
            return supercontig._seq_dtype
        return DEFAULT_SEQ_DTYPE

    @property
    def name(self):
        """Return the name of this chromosome (same as __str__())."""
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
        """See :attr:`Genome.mins`"""
        return self.attrs.mins

    @property
    def maxs(self):
        """See :attr:`Genome.maxs`"""
        return self.attrs.maxs

    @property
    def sums(self):
        """See :attr:`Genome.sums`"""
        return self.attrs.sums

    @property
    def sums_squares(self):
        """See :attr:`Genome.sums_squares`"""
        return self.attrs.sums_squares

    @property
    def num_datapoints(self):
        """See :attr:`Genome.num_datapoints`"""
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

        :returns: :class:`Supercontig`

        Indexable with a slice or simple index:

        >>> chromosome.supercontigs[100]
        [<Supercontig('supercontig_0', 0:66115833)>]
        >>> chromosome.supercontigs[1:100000000]
        [<Supercontig('supercontig_0', 0:66115833)>, <Supercontig('supercontig_1', 66375833:90587544)>, <Supercontig('supercontig_2', 94987544:199501827)>]
        >>> chromosome.supercontigs[66115833:66375833]  # Between two supercontigs
        []

        """
        return self._supercontigs

class Supercontig(object):
    """A container for a segment of data in one chromosome.

    Implemented via a HDF5 Group

    """
    def __init__(self, h5group):
        """
        :param h5group: group containing the supercontig data
        :type h5group: HDF5 Group

        """
        self.h5group = h5group

    def __repr__(self):
        return "<Supercontig('%s', %d:%d)>" % (self.name, self.start, self.end)

    def __str__(self):
        return str(self.name)

    def project(self, pos):
        """Project chromsomal coordinates to supercontig coordinates.

        :param pos: chromosome coordinate
        :type pos: integer
        :returns: integer

        """
        return pos - self.start

    @property
    def _seq_dtype(self):
        try:
            return self.seq.atom.dtype
        except NoSuchNodeError:
            return DEFAULT_SEQ_DTYPE

    @property
    def _continuous_dtype(self):
        try:
            return self.continuous.atom.dtype
        except NoSuchNodeError:
            return DEFAULT_CONTINUOUS_DTYPE

    @property
    def continuous(self):
        """Return the underlying continuous data in this supercontig.

        :returns: numpy.array

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

        # If index was specific, don't return an array
        key_int = False
        if isinstance(key, int):
            key_int = True

        start, end = _key_to_tuple(key)
        length = end - start
        dtype = self._chromosome._seq_dtype
        if length <= 0:  # Handle degenerate case quickly
            return array((), dtype=dtype)

        seq = empty((length,), dtype=dtype)
        seq.fill(ord("N"))  # Assumes dtype is numeric type
        for supercontig in supercontigs:
            chr_start = max(start, supercontig.start)
            chr_end = min(end, supercontig.end)
            seq[chr_start - start:chr_end - start] = \
                supercontig.seq[supercontig.project(chr_start):
                                    supercontig.project(chr_end)]

        if key_int:
            seq = seq[0]
        return seq

class _Supercontigs(object):
    def __init__(self, chromosome):
        self._chromosome = chromosome

    def __getitem__(self, key):
        """Return list of supercontigs containing any of this index range"""
        start, end = _key_to_tuple(key)
        if start < self._chromosome.start:
            start = self._chromosome.start
        if end > self._chromosome.end:
            end = self._chromosome.end

        supercontigs = []
        for supercontig in self._chromosome:
            if start < supercontig.end and end > supercontig.start:
                supercontigs.append(supercontig)
                if start >= supercontig.start and end <= supercontig.end:
                    # Key entirely within one supercontig, so we're done
                    break
            # XXX: would be nice if we could count on supercontig ordering

        return supercontigs

def _key_to_tuple(key):
    """Key to (start, stop)"""
    if isinstance(key, int):
        start = key
        end = key + 1
    elif isinstance(key, slice):
        if key.start is None or key.stop is None:
            raise NotImplementedError("Both start and end must be specified in"
                                      " chromosomal slices")
        elif key.step is not None and key.step != 0:
            raise NotImplementedError("Chromosome slicing does not support"
                                      " non-contiguous retrieval")
        else:
            start = key.start
            end = key.stop
    else:
        raise NotImplementedError("Unsupported index found: %s" % key)

    if start < 0 or end < 0:
        raise NotImplementedError("Negative (wrapping) indices not supported")
    elif start > end:
        raise IndexError("Start index can be at most the end index")

    return start, end


def main(args=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())

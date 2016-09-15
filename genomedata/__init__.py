#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
Genomedata is a module to store and access large-scale functional
genomics data in a format which is both space-efficient and allows
efficient random-access.

Under the surface, genomedata is implemented as a collection of HDF5 files,
but genomedata provides a transparent interface to interact with your
underlying data without having to worry about the mess of repeatedly parsing
large data files or having to keep them in memory for random access.

Copyright 2009-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

"""

__version__ = "1.3.6"


import sys

import tables
from functools import partial
from numpy import (add, amin, amax, append, array, empty, float32, inf,
                   nan, ndarray, square, uint8)
from os import extsep
from path import path
from tables import Float32Atom, NoSuchNodeError, open_file, UInt8Atom
from warnings import warn

FORMAT_VERSION = 1
SEQ_DTYPE = uint8
SEQ_ATOM = UInt8Atom()

CONTINUOUS_DTYPE = float32
CONTINUOUS_ATOM = Float32Atom(dflt=nan)
CONTINUOUS_CHUNK_SHAPE = (10000, 1)

EXT = "genomedata"
SUFFIX = extsep + EXT

# If there are fewer than this many chomosomes, the default behavior
# is to implement the Genomedata archive as a directory. If there are
# more than this many, it will be a single file by default.
FILE_MODE_CHROMS = 100

class _InactiveDict(dict):
    """A fake dict that can't be added to."""
    def __setitem__(self, key, value):
        return


def _open_file(filename, *args, **kwargs):
    # From pytables 3 docs:
    # [open_file] recognizes the (lowercase) names of parameters present in tables/parameters.py
    if not "buffer_times" in kwargs:
        # eliminate spurious PerformanceWarning
        kwargs["buffer_times"] = inf

    return open_file(str(filename), *args, **kwargs)


class Genome(object):
    """The root level of the genomedata object hierarchy.

    If you use this as a context manager, it will keep track of any open
    Chromosomes and close them (and the Genome object) for you later when
    the context is left::

      with Genome("/path/to/genomedata") as genome:
        chromosome = genome["chr1"]
        [...]

    If not used as a context manager, you are responsible for closing
    the Genomedata archive once you are done:

    >>> genome = Genome("/path/to/genomedata")
    >>> chromosome = genome["chr1"]
    [...]
    >>> genome.close()

    """
    def __init__(self, filename, *args, **kwargs):
        """Create a Genome object from a genomdata archive.

        :param filename: the root of the Genomedata object
                         hierarchy. This can either be a .genomedata
                         file that contains the entire genome or a
                         directory containing multiple chromosome files.
        :type filename: string
        :param \*args: args passed on to open_file if single file or to
                       Chromosome if directory
        :param \*\*kwargs: keyword args passed on to open_file if single file
                           or to Chromosome if directory

        Example:

        >>> genome = Genome("./genomedata.ctcf.pol2b/")
        >>> genome
        Genome("./genomedata.ctcf.pol2b/")
            [...]
        >>> genome.close()
        >>> genome = Genome("./cat_chipseq.genomedata", mode="r")
            [...]
        >>> genome.close()

        """
        # so that Genome.__del__() won't throw an exception if there
        # is an error during __init__()
        self._isopen = False
        self.filename = filename
        self.args = args
        self.kwargs = kwargs

        # Process path for internal use, following symbolic links
        # until we get to the eventual file or directory
        filepath = path(filename).expand()
        while filepath.islink():
            filepath = filepath.readlinkabs()

        if not filepath.exists():
            raise IOError("Could not find Genomedata archive: %s" % filepath)

        if filepath.isfile():
            # Open the Genomedata file
            isfile = True
            self.h5file = _open_file(filepath, *args, **kwargs)
            self._file_attrs = self.h5file.root._v_attrs
        elif filepath.isdir():
            # Genomedata directory
            isfile = False
        else:
            raise ValueError("Genomedata archive must be file or directory: %s"
                             % filepath)

        self._path = filepath
        self._isfile = isfile
        # Keep track of open chromosomes
        self.open_chromosomes = {}
        # a kind of refcounting for context managers
        self._context_count = 0
        self._isopen = True

        format_version = self.format_version
        if format_version is not None and format_version > FORMAT_VERSION:
            raise NotImplementedError("This archive has format version %s,"
                                      " but the installed Genomedata software"
                                      " unly supports format version %d"
                                      % (format_version, FORMAT_VERSION))

    def __iter__(self):
        """Return next chromosome, in sorted order, with memoization.

        Example::

          for chromosome in genome:
            print chromosome.name
            for supercontig, continuous in chromosome.itercontinuous():
              [...]

        """
        assert self.isopen
        if self._isfile:  # Chromosomes are groups
            # Iterate over child group of root
            for group in self.h5file.iter_nodes("/", classname="Group"):
                groupname = group._v_name
                yield self[groupname]
        else:  # Chromosomes are files
            # sorted so that the order is always the same
            for filepath in sorted(self._path.files("*" + SUFFIX)):
                # pass through __getitem__() to allow memoization
                yield self[filepath.namebase]

    def __getitem__(self, name):
        """Return a reference to a chromosome of the given name.

        :param name: name of the chromosome (e.g. "chr1" if
                     chr1.genomedata is a file in the Genomedata archive
                     or chr1 is a top-level group in the single-file
                     Genomedata archive)
        :type name: string
        :returns: :class:`Chromosome`

        Example:

        >>> genome["chrX"]
        <Chromosome 'chrX', file='/path/to/genomedata/chrX.genomedata'>
        >>> genome["chrZ"]
        KeyError: 'Could not find chromosome: chrZ'

        """
        assert self.isopen
        try:
            # memoization
            return self.open_chromosomes[name]
        except KeyError:
            pass

        try:
            if self._isfile:
                res = Chromosome(self.h5file, where="/" + name)
            else:
                res = Chromosome._fromfilename(
                    self._path.joinpath(name + SUFFIX),
                    *self.args, **self.kwargs)
        except (IOError, NoSuchNodeError):
            raise KeyError("Could not find chromosome: %s" % name)

        self.open_chromosomes[name] = res
        return res

    def __contains__(self, name):
        """Return if there is a chromosome of the given name

        :param name: name of the chromosome (e.g. "chr1" if
                     chr1.genomedata is a file in the Genomedata archive
                     or chr1 is a top-level group in the single-file
                     Genomedata archive)
        :type name: string
        :returns: boolean

        Example:

        >>> "chrX" in Genome
        True
        >>> "chrZ" in Genome
        False

        """
        try:
            self[name]
        except KeyError:
            return False  # Couldn't find chromosome
        else:
            return True  # No errors opening chromosome

    def __enter__(self):
        assert self.isopen
        self._context_count += 1
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        # XXX: this and __enter__ have potential race conditions, if
        # _context_count is changed simultaneously by different
        # threads. should be synchronized
        self._context_count -= 1
        if self._context_count == 0:
            self.close()

    def __del__(self):
        if self.isopen:
            self.close()

    def close(self):
        """Close this Genomedata archive and any open chromosomes

        If the Genomedata archive is a directory, this closes all open
        chromosomes. If it is a single file, this closes that file.
        This should only be used if Genome is not a context manager
        (see :class:`Genome`). The behavior is undefined if this is
        called while Genome is being used as a context manager.

        """
        assert self.isopen

        # Whether a single file or a directory, close all the chromosomes
        # so they know they shouldn't be read. Do this before closing
        # Genome.h5file in case the chromosomes need access to it in closing.
        for name, chromosome in self.open_chromosomes.iteritems():
            # Only close those not closed manually by the user
            if chromosome.isopen:
                chromosome.close()

        if self._isfile:
            self.h5file.close()

        self.open_chromosomes = {}
        self._isopen = False

    def __repr__(self):
        items = ["'%s'" % self.filename]
        if self.args:
            items.append("*%r" % self.args)
        if self.kwargs:
            items.append("**%r" % self.kwargs)
        return "Genome(%s)" % ", ".join(items)

    def __str__(self):
        return repr(self)

    def _accum_extrema(self, name, accumulator):
        self.tracknames_continuous  # for assertion check

        extrema = [getattr(chromosome, name) for chromosome in self]
        return accumulator(extrema)

    def erase_data(self, trackname):
        """Erase all data for the given track across all chromosomes

        The Genome object must have been created with
        :param mode:="r+". Behavior is undefined if this is not the case.

        Currently sets the dirty bit, which can only be erased with
        genomedata-close-data

        """
        assert self.isopen
        for chromosome in self:
            chromosome._erase_data(trackname)

    def add_track_continuous(self, trackname):
        """Add a new track

        The Genome object must have been created with
        :param mode:="r+". Behavior is undefined if this is not the case.

        Currently sets the dirty bit, which can only be erased with
        genomedata-close-data

        """
        assert self.isopen
        if self.format_version < 1:
            raise NotImplementedError("""Adding tracks is only supported \
for archives created with Genomedata version 1.2.0 or later.""")

        if self._isfile:
            # Update tracknames attribute on file
            attrs = self._file_attrs
            if "tracknames" in attrs:
                tracknames = attrs.tracknames
                if trackname in tracknames:
                    raise ValueError("%s already has a track of name: %s"
                                     % (self.filename, trackname))
            else:
                tracknames = array([])

            attrs.tracknames = append(tracknames, trackname)

        # Let the chromosomes handle the rest
        for chromosome in self:
            chromosome._add_track_continuous(trackname)

    @property
    def isopen(self):
        """Return a boolean indicating if the Genome is still open"""
        return self._isopen

    @property
    def tracknames_continuous(self):
        """Return a list of the names of all data tracks stored."""
        assert self.isopen
        if self._isfile:
            # Tracknames are stored at the root of each file, so we can
            # access them directly in this case
            return self._file_attrs.tracknames.tolist()
        else:
            # check that all chromosomes have the same tracknames_continuous
            res = None
            for chromosome in self:
                if res is None:
                    res = chromosome.tracknames_continuous
                else:
                    assert res == chromosome.tracknames_continuous

        return res

    def index_continuous(self, trackname):
        """Return the column index of the trackname in the continuous data.

        :param trackname: name of data track
        :type trackname: string
        :returns: integer

        This is used for efficient indexing into continuous data:

        >>> col_index = genome.index_continuous("sample_track")
        >>> data = genome["chr3"][100:150, col_index]

        although for typical use, the track can be indexed directly:

        >>> data = genome["chr3"][100:150, "sample_track"]

        """
        try:
            return self.tracknames_continuous.index(trackname)
        except ValueError:
            raise KeyError("Could not find continuous track: %s" % trackname)

    @property
    def num_tracks_continuous(self):
        """Returns the number of continuous data tracks."""
        try:
            return len(self.tracknames_continuous)
        except AttributeError:
            return 0

    @property
    def format_version(self):
        """Genomedata format version

        None means there are no chromosomes in it already.
        """
        assert self.isopen
        if self._isfile:
            try:
                return self._file_attrs.genomedata_format_version
            except AttributeError:
                pass

        # else: self is a directory
        chromosomes = iter(self)
        try:
            first_chromosome = next(chromosomes)
        except StopIteration:
            return None

        res = first_chromosome._format_version

        assert all(res == chromosome._format_version
                   for chromosome in chromosomes)

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
        return (self.sums_squares / self.num_datapoints) - \
            square(self.means)


class Chromosome(object):
    """The Genomedata object corresponding to data for a given chromosome.

    Usually created by keying into a Genome object with the name of a
    chromosome, as in:

    >>> with Genome("/path/to/genomedata") as genome:
    ...     chromosome = genome["chrX"]
    ...     chromosome
    ...
    <Chromosome 'chrX', file='/path/to/genomedata/chrX.genomedata'>

    """
    class ChromosomeDirtyError(Exception):
        pass

    default_where = "/"
    default_mode = "r"

    def __init__(self, h5file, where=default_where, name=None):
        """
        :param h5file: tables.File object for the h5file which contains
                       the chromosome to be opened. If the Genomedata archive
                       is a single file, :param where: should specify the
                       path to the chromosome group within the file.
        :param where: path or Node to the root of the chromosome within
                      the Genomedata file.
        :param name: name of the Chromosome. If None, the name will try
                     to be parsed from :param where:.
        :type file: tables.File
        :type where: string or tables.Node
        :type name: string or None

        """
        # If file is a string, open the h5 file
        if isinstance(h5file, tables.File):
            if name is None:
                name = where.rpartition("/")[2]
        else:
            raise NotImplementedError("Chromosome file of unsupported"
                                      " type: %r" % h5file)

        # Now, open the group that is the root of the chromosome
        h5group = h5file.get_node(where, classname="Group")

        # XXX: even though each chromosome has its own dirty bit, and
        # the metadata only needs to be recalculated on those where it is
        # set, right now we can't guarantee that the user hasn't changed
        # the values directly, so we have to set the dirty bit on every
        # opened chromosome group. This can be improved by tracking what
        # the user changes.
        attrs = h5group._v_attrs
        if h5file.mode in set(["w", "r+", "a"]):
            # Make sure there is a genomedata_format_version
            file_attrs = h5file.root._v_attrs
            if "genomedata_format_version" not in file_attrs:
                # Set as first version (before it was standard)
                file_attrs.genomedata_format_version = 0

            attrs.dirty = True
        else:
            if attrs.dirty:
                raise self.ChromosomeDirtyError("""
Chromosome has been modified (or loaded with a mode of "w", "r+", or "a")
since being closed with genomedata-close-data.""")

        self.filename = h5file.filename
        self._name = name
        self.h5file = h5file
        self.h5group = h5group
        self._isfile = (where == self.default_where)
        self._seq = _ChromosomeSeqSlice(self)
        self._supercontigs = _Supercontigs(self)
        self._isopen = True

    @classmethod
    def _fromfilename(cls, filename, mode=default_mode, *args, **kwargs):
        """
        :param filename: name of the chromosome (.genomedata) file to access

        :param mode: mode of interaction with the chromosome file,
                     with ``r``: read, ``w``: write, ``a``: append,
                     ``r+``: append but force file to exist (see documentation
                     for tables.open_file().)

        :type mode: string
        :param \*args: args passed on to open_file
        :param \*\*kwargs: keyword args passed on to open_file

        """
        filepath = path(filename).expand()
        try:
            h5file = _open_file(filepath, mode=mode, *args, **kwargs)
        except IOError:
            raise IOError("Could not find file: %r" % filename)

        name = filepath.name.rpartition(SUFFIX)[0]

        return cls(h5file, name=name)

    def __iter__(self):
        """Return next supercontig in chromosome.

        .. versionadded:: 1.2
           Supercontigs are ordered by start index

        Seldom used in favor of the more direct:
        :meth:`Chromosome.itercontinuous`

        Example:

        >>> for supercontig in chromosome:
        ...     supercontig  # calls repr()
        ...
        <Supercontig 'supercontig_0', [0:66115833]>
        <Supercontig 'supercontig_1', [66375833:90587544]>
        <Supercontig 'supercontig_2', [94987544:199501827]>

        """
        assert self.isopen
        supercontigs = []
        for group in self.h5group:
            supercontig = Supercontig(group)
            supercontigs.append((supercontig.start, supercontig))

        supercontigs.sort()
        for start, supercontig in supercontigs:
            yield supercontig

    def __getitem__(self, key):
        """Return the continuous data corresponding to this bp slice

        :param key: base_key must index or slice bases
                    track_key specify data tracks with index, slice, string,
                    list of strings, list of indexes, or array of indexes

        but can also index, slice,
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
        assert self.isopen
        # XXX: Allow variable/negative steps, negative starts/stops, etc.

        # XXX: The problem with missing start or end indices is that
        # 1) it is unclear if 0 or self.start should be used for the
        #    start, and (preference for 0)
        # 2) the full length of the chromosome is not known, so self.end
        #    would need to be used, which might not be what is wanted. -OJB
        # XXX: I think this should no longer be a problem now that
        # self.end is the full length of the chromsome? need to check this -MMH

        # Sanitize the input
        if isinstance(key, tuple):
            base_key, track_key = key
        else:
            base_key = key
            track_key = slice(None)  # All tracks

        # just like NumPy, direct indexing results in output shape
        # change (at end of method)
        base_direct_index = isinstance(base_key, int)
        track_direct_index = isinstance(track_key, (str, int))

        # convert base_key
        base_key = slice(*_key_to_tuple(base_key))

        # First convert track_key toward slice
        if isinstance(track_key, (list, ndarray)):
            track_indexes = array([self._index_continuous(item)
                                   for item in track_key])
            track_min = track_indexes.min()
            track_key = slice(track_min, track_indexes.max() + 1, 1)
            track_subset_indexes = track_indexes - track_min

        else:
            track_subset_indexes = slice(None)  # everything
            if isinstance(track_key, str):
                track_key = self.index_continuous(track_key)
            if isinstance(track_key, int):
                track_key = slice(track_key, track_key + 1, 1)

        if isinstance(track_key, slice):
            # Fix indices to number of tracks
            track_key = slice(*track_key.indices(self.num_tracks_continuous))
        else:
            raise TypeError("Unrecognized track indexing method: %s" %
                            track_key)

        nrows = base_key.stop - base_key.start
        ncols = len(xrange(track_key.start, track_key.stop, track_key.step))
        dtype = self._continuous_dtype

        # Handle degenerate case
        if nrows < 1 or ncols < 1:
            # Return empty array (matches numpy behavior)
            return array([], dtype=dtype)

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
        data.fill(nan)

        for supercontig in supercontigs:
            assert (base_key.start < supercontig.end and
                    base_key.stop > supercontig.start)
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

        # get a subset of tracks
        data = data[:, track_subset_indexes]

        # Make output shape appropriate for indexing method (like numpy)
        if track_direct_index:
            data = data[:, 0]
        if base_direct_index:
            data = data[0]
        return data

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        return "<Chromosome '%s', file='%s'>" % (self.name, self.filename)

    def itercontinuous(self):
        """Return a generator over all supercontig, continuous pairs.

        .. versionadded:: 1.2
           Supercontigs are ordered by increasing supercontig.start.

        This is the best way to efficiently iterate over the data since
        all specified data is in supercontigs::

            for supercontig, continuous in chromosome.itercontinuous():
                print supercontig, supercontig.start, supercontig.end
                [...]

        """
        assert self.isopen
        for supercontig in self:
            try:
                yield supercontig, supercontig.continuous
            except NoSuchNodeError:
                continue

    def _index_continuous(self, track_key):
        """
        Convert track_key to index only when it is a basestring.
        Otherwise return track_key unchanged.
        """
        if isinstance(track_key, str):
            return self.index_continuous(track_key)

        return track_key

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
        assert self.isopen
        try:
            return self.tracknames_continuous.index(trackname)
        except ValueError:
            raise KeyError("Could not find continuous track: %s" % trackname)

    def close(self):
        """Close the current chromosome file.

        This only needs to be called when Genomedata files are manually
        opened as Chromosomes. Otherwise, :meth:`Genome.close`
        should be called to close any open chromosomes or Genomedata files.
        The behavior is undefined if this is called on a Chromosome accessed
        through a Genome object.
        Using Genomedata as a context manager makes
        life easy by memoizing chromosome access and guaranteeing the
        proper cleanup. See :class:`Genome`.

        """
        assert self.isopen

        if self.attrs.dirty:
            warn("Closing Chromosome with modified data. Metadata needs to"
                 " be recalculated by calling genomedata-close-data on the"
                 " Genomedata archive before re-accessing it")

        if self._isfile:
            self.h5file.close()

        self._isopen = False

    def _erase_data(self, trackname):
        """Erase all data for the given track

        The Genome object or this Chromosome must have been created with
        :param mode:="r+". Behavior is undefined if this is not the case.

        Currently sets the dirty bit, which can only be erased with
        genomedata-close-data GENOMEDATA

        """
        assert self.isopen
        col_index = self.index_continuous(trackname)
        self.attrs.dirty = True
        for supercontig, continuous in self.itercontinuous():
            continuous[:, col_index] = nan

    def _add_track_continuous(self, trackname):
        """Add a new track

        The Genome object must have been created with
        :param mode:="r+". Behavior is undefined if this is not the case.

        Currently sets the dirty bit, which can only be erased with
        genomedata-close-data

        """
        assert self.isopen
        if self._isfile:
            # Update tracknames attribute with new trackname
            file_attrs = self._file_attrs
            if "tracknames" in file_attrs:
                tracknames = file_attrs.tracknames
                if trackname in tracknames:
                    raise ValueError("%s already has a track of name: %s"
                                     % (self.filename, trackname))
            else:
                tracknames = array([])

            file_attrs.tracknames = append(tracknames, trackname)
        # else: hope the Genome object updated its own tracknames

        self.attrs.dirty = True  # dirty specific to chromosome

        # Extend supercontigs by a column (or create them)
        for supercontig in self:
            supercontig_length = supercontig.end - supercontig.start
            try:
                continuous = supercontig.continuous
            except NoSuchNodeError:
                # Define an extendible array in the second dimension (0)
                supercontig_shape = (supercontig_length, 0)
                self.h5file.create_earray(supercontig.h5group, "continuous",
                                          CONTINUOUS_ATOM, supercontig_shape,
                                          chunkshape=CONTINUOUS_CHUNK_SHAPE)
                continuous = supercontig.continuous

            # Add column to supercontig continuous array
            # "truncate" also extends with default values
            continuous.truncate(continuous.nrows + 1)

    @property
    def isopen(self):
        """Return a boolean indicating if the Chromosome is still open"""
        return self._isopen

    @property
    def _continuous_dtype(self):
        for supercontig, continuous in self.itercontinuous():
            return supercontig._continuous_dtype
        return CONTINUOUS_DTYPE

    @property
    def _seq_dtype(self):
        for supercontig in self:
            return supercontig._seq_dtype
        return SEQ_DTYPE

    @property
    def name(self):
        """Return the name of this chromosome (same as __str__())."""
        return self._name

    @property
    def attrs(self):
        """Return the attributes for this Chromosome.

        This may also include Genome-wide attributes if the archive
        is implemented as a directory.

        """
        assert self.isopen
        return self.h5group._v_attrs

    @property
    def _file_attrs(self):
        assert self.isopen
        return self.h5file.root._v_attrs

    @property
    def tracknames_continuous(self):
        """Return a list of the data track names in this Chromosome."""
        assert self.isopen
        return self._file_attrs.tracknames.tolist()

    @property
    def num_tracks_continuous(self):
        """Return the number of tracks in this chromosome"""
        try:
            return len(self.tracknames_continuous)
        except AttributeError:
            return 0

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
    def _format_version(self):
        """See :attr:`Genome.format_version`"""
        try:
            return self._file_attrs.genomedata_format_version
        except AttributeError:
            try:
                return self.attrs.genomedata_format_version
            except AttributeError:
                # original version did not have
                # genomedata_format_version attribute
                return 0

    @property
    def start(self):
        """Return the index of the first base in this chromosome.

        For :attr:`Genome.format_version` > 0, this will always be 0.
        For == 0, this will be the start of the first supercontig.

        """
        if self._format_version == 0:
            return min(supercontig.start for supercontig in self)
        else:
            return self.attrs.start

    @property
    def end(self):
        """Return the index past the last base in this chromosome.

        For :attr:`Genome.format_version` > 0, this will be
        the number of bases of sequence in the chromosome. For == 0,
        this will be the end of the last supercontig.

        This is the end in half-open coordinates, making slicing simple:

        >>> chromosome.seq[chromosome.start:chromosome.end]

        """
        if self._format_version == 0:
            return max(supercontig.end for supercontig in self)
        else:
            return self.attrs.end

    @property
    def seq(self):
        """Return the genomic sequence of this chromosome.

        If the index or slice spans a non-supercontig range, N's are
        inserted in place of the missing data and a warning is issued.

        Example:

        >>> chromosome = genome["chr1"]
        >>> for supercontig in chromosome:
        ...     print repr(supercontig)
        ...
        <Supercontig 'supercontig_0', [0:121186957]>
        <Supercontig 'supercontig_1', [141476957:143422081]>
        <Supercontig 'supercontig_2', [143522081:247249719]>
        >>> chromosome.seq[0:10].tostring()  # Inside supercontig
        'taaccctaac'
        >>> chromosome.seq[121186950:121186970].tostring()  \
# supercontig boundary
        'agAATTCNNNNNNNNNNNNN'
        >>> chromosome.seq[121186957:121186960].tostring()  \
# not in supercontig
        UserWarning: slice of chromosome sequence does not overlap any \
supercontig (filling with 'N')
        'NNN'

        The entire sequence for a chromosome can be retrieved with:

        >>> chromosome.seq[chromosome.start:chromosome.end]

        """
        return self._seq

    @property
    def supercontigs(self):
        """Return the supercontig that contains this range if possible.

        :returns: :class:`Supercontig`

        Indexable with a slice or simple index:

        >>> chromosome.supercontigs[100]
        [<Supercontig 'supercontig_0', [0:66115833]>]
        >>> chromosome.supercontigs[1:100000000]
        [<Supercontig 'supercontig_0', [0:66115833]>, \
<Supercontig 'supercontig_1', [66375833:90587544]>, \
<Supercontig 'supercontig_2', [94987544:199501827]>]
        >>> chromosome.supercontigs[66115833:66375833]  \
# Between two supercontigs
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
        return "<Supercontig '%s', [%d:%d]>" % (self.name, self.start,
                                                self.end)

    def __str__(self):
        return str(self.name)

    def project(self, pos, bound=False):
        """Project chromsomal coordinates to supercontig coordinates.

        :param pos: chromosome coordinate
        :param bound: bound result to valid supercontig coordinates
        :type pos: integer
        :type bound: boolean
        :returns: integer

        """
        if bound:
            pos = max(pos, self.start)
            pos = min(pos, self.end)

        return int(pos - self.start)

    @property
    def _seq_dtype(self):
        try:
            return self.seq.atom.dtype
        except NoSuchNodeError:
            return SEQ_DTYPE

    @property
    def _continuous_dtype(self):
        try:
            return self.continuous.atom.dtype
        except NoSuchNodeError:
            return CONTINUOUS_DTYPE

    @property
    def continuous(self):
        """Return the underlying continuous data in this supercontig.
        To read the whole dataset into memory as a `numpy.array`, use
        continuous.read()

        :returns: `tables.EArray`

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
        """See :attr:`Chromosome.seq`."""
        return self.h5group.seq

    @property
    def start(self):
        """Return the index of the first base in this supercontig.

        The first base is index 0.

        """
        return int(self.attrs.start)

    @property
    def end(self):
        """Return the index past the last base in this supercontig.

        This is the end in half-open coordinates, making slicing simpler:

        >>> supercontig.seq[supercontig.start:supercontig:end]

        """
        return int(self.attrs.end)


class _ChromosomeSeqSlice(object):
    def __init__(self, chromosome):
        assert isinstance(chromosome, Chromosome)
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
            return array([], dtype=dtype)

        seq = empty((length,), dtype=dtype)
        seq.fill(ord("N"))  # Assumes dtype is numeric type
        for supercontig in supercontigs:
            chr_start = max(start, supercontig.start)
            chr_end = min(end, supercontig.end)
            dest_start = chr_start - start
            dest_end = chr_end - start
            sc_start = supercontig.project(chr_start)
            sc_end = supercontig.project(chr_end)
            seq[dest_start:dest_end] = supercontig.seq[sc_start:sc_end]

        if key_int:
            seq = seq[0]
        return seq


class _Supercontigs(object):
    def __init__(self, chromosome):
        assert isinstance(chromosome, Chromosome)
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


def main(argv=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())

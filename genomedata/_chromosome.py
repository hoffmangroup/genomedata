from functools import partial
from operator import attrgetter
from warnings import warn

from numpy import (add, amin, amax, array, empty, float32, nan, ndarray,
                   uint8)
from six.moves import range

from ._util import OverlapWarning

CONTINUOUS_DTYPE = float32

SEQ_DTYPE = uint8


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

    default_mode = "r"

    def __init__(self, name=None):
        self._name = name
        self._seq = _ChromosomeSeqSlice(self)
        self._supercontigs = _Supercontigs(self)
        self._isopen = True

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
        pass  # NB: implemented by subclass

    def _check_region_inside_supercontigs(self, supercontigs, start, end):
        """
        Checks if the given region given by start and end overlap with
        regions not covered by the given supercontigs
        """

        num_supercontigs = len(supercontigs)
        # If there are no supercontigs
        if num_supercontigs == 0:
            raise ValueError("{} {} {} sequence does not overlap any "
                             "supercontig (nothing written)".format(
                                 self._name, start, end))

        # If there is at least 1 supercontig
        elif num_supercontigs > 0:
            # If there is more than 1 supercontig
            if num_supercontigs > 1:
                # If there is a gap between any given supercontig
                sorted_contigs = sorted(supercontigs, key=attrgetter("start"))
                if any(((sorted_contigs[i].start - sorted_contigs[i-1].end) > 0
                        for i in range(1, num_supercontigs))):
                    # Do not support writing between gaps in supercontigs
                    raise ValueError("{} {} {} sequence overlaps gaps in"
                                     "supercontigs (nothing written)".format(
                                         self._name,
                                         start,
                                         end))

            # If the specified chromosomal base start index is before the
            # first supercontig coordinate
            first_supercontig = min(supercontigs, key=attrgetter("start"))
            if start < first_supercontig.start:
                # Do not support writing before the first supercontig
                raise ValueError("{} {} {} sequence overlaps before beginning "
                                 " of earliest supercontig (nothing "
                                 "written)".format(
                                     self._name,
                                     start,
                                     end))

            # If the specified chromosomal end index is after the furthest
            # supercontig coordinate
            last_supercontig = max(supercontigs, key=attrgetter("end"))
            if end > last_supercontig.end:
                # Do not support writing after the last supercontig
                raise ValueError("{} {} {} sequence overlaps after end of "
                                 "furthest supercontig (nothing "
                                 "written)".format(
                                     self._name,
                                     start,
                                     end))

    def _get_base_and_track_key(self, key):
        """
        Takes a given key indexed into a chromsome and splits the key into
        a base and track key if necessary
        """
        # Sanitize the input
        if isinstance(key, tuple):
            base_key, track_key = key
        else:
            base_key = key
            track_key = slice(None)  # All tracks

        return base_key, track_key

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

        # Get the chromosomal base key and track key
        base_key, track_key = self._get_base_and_track_key(key)

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
        ncols = len(range(track_key.start, track_key.stop, track_key.step))
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
                 " (filling with 'NaN')", category=OverlapWarning)
        elif len(supercontigs) > 1:
            warn("slice of chromosome data spans more than one supercontig"
                 " (filling gaps with 'NaN')", category=OverlapWarning)

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
            except MissingContinuousData:
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

    def __setitem__(self, key, value):
        if (not self.isopen or
           self.mode != "r+"):  # r+ is the only open mode for writing
            raise IOError("Genomedata archive not opened for writing")

        # Split the given key to its chromosomal and track key
        base_key, track_key = self._get_base_and_track_key(key)

        # Convert base key to slice for indexing
        base_key = slice(*_key_to_tuple(base_key))

        # Convert track key to numbered list if necessary
        if isinstance(track_key, (list, ndarray)):
            track_key = array([self._index_continuous(item)
                               for item in track_key])

        # Get a list of supercontigs from our base key
        supercontigs = self.supercontigs[base_key]

        # Check if the given region overlaps with any gap in the assembly
        self._check_region_inside_supercontigs(supercontigs,
                                               base_key.start,
                                               base_key.stop)

        # For each supercontig
        for supercontig in supercontigs:

            # Ensure the chromosomal base indicies overlap with the
            # supercontig region
            # XXX: This is the same as in read and this should be
            # guaranteed at this point
            assert (base_key.start < supercontig.end and
                    base_key.stop > supercontig.start)

            # Cap the chromosomal coordinates to the supercontig start and
            # end if necessary
            chr_start = max(base_key.start, supercontig.start)
            chr_end = min(base_key.stop, supercontig.end)

            # Get the indexes for the underlying continuous data
            supercontig_slice = slice(supercontig.project(chr_start),
                                      supercontig.project(chr_end))

            # Write values to supercontig
            supercontig.continuous[supercontig_slice, track_key] = value

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
            except MissingContinuousData:
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
        for supercontig, continuous in self.itercontinuous():
            continuous[:, col_index] = nan

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
    def mode(self):
        """Return the read/write properties for this chromosome as string"""
        return Chromosome.default_mode

    @property
    def num_tracks_continuous(self):
        """Return the number of tracks in this chromosome"""
        try:
            return len(self.tracknames_continuous)
        except AttributeError:
            return 0

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
        OverlapWarning: slice of chromosome sequence does not overlap any \
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


class _ChromosomeList(object):
    def __init__(self):
        self.open_chromosomes = {}  # NB: maintained by subclass

    def __getitem__(self, name):
        """Returns an existing opened chromosome if it exists. Subclasses
        assumed to implement this method and call this for caching."""
        # If we have opened the chromosome already
        if (name in self.open_chromosomes and
           self.open_chromosomes[name].isopen):  # NB: manual closure check
            # Return cached object
            return self.open_chromosomes[name]

        return None

    def _format_version(self):
        """Get version information metadata. None by default. Subclasses to
        implement proper versioning if possible"""
        return None

    def _accum_extrema(self, name, accumulator):
        self.tracknames_continuous  # for assertion check

        extrema = [getattr(chromosome, name) for chromosome in self]
        return accumulator(extrema)

    def tracknames_continuous(self):
        raise NotImplementedError

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

    def close(self):
        self.open_chromosomes = {}


class Supercontig(object):
    """A container for a segment of data in one chromosome.

    """
    def __repr__(self):
        return "<Supercontig '%s', [%d:%d]>" % (self.name, self.start,
                                                self.end)

    def __str__(self):
        return str(self.name)

    def project(self, pos, bound=False):
        """Project chromosomal coordinates to supercontig coordinates.

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
        raise NotImplementedError

    @property
    def _continuous_dtype(self):
        raise NotImplementedError

    @property
    def continuous(self):
        """Return the underlying continuous data in this supercontig.
        To read the whole dataset into memory as a `numpy.array`, use
        continuous.read()

        """
        raise NotImplementedError

    @property
    def attrs(self):
        """Return the attributes of this supercontig."""
        raise NotImplementedError

    @property
    def name(self):
        """Return the name of this supercontig."""
        raise NotImplementedError

    @property
    def seq(self):
        """See :attr:`Chromosome.seq`."""
        raise NotImplementedError

    @property
    def start(self):
        """Return the index of the first base in this supercontig.

        The first base is index 0.

        """
        raise NotImplementedError

    @property
    def end(self):
        """Return the index past the last base in this supercontig.

        This is the end in half-open coordinates, making slicing simpler:

        >>> supercontig.seq[supercontig.start:supercontig:end]

        """
        raise NotImplementedError


class MissingContinuousData(Exception):
    """Raised when a continuous section of a supercontig does not exist"""
    pass


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
                 " supercontig (filling with 'N')", OverlapWarning)
        elif len(supercontigs) > 1:
            warn("slice of chromosome sequence spans more than one supercontig"
                 " (filling gaps with 'NaN')", OverlapWarning)

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

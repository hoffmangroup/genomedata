from warnings import warn

from numpy import array, append, inf, nan
from path import Path
import tables

from ._chromosome import (CONTINUOUS_DTYPE, Chromosome, _ChromosomeList,
                          MissingContinuousData, Supercontig, SEQ_DTYPE)
from ._util import (decode_trackname, GenomedataDirtyWarning,
                    GENOMEDATA_ENCODING, SUFFIX)

CONTINUOUS_ATOM = tables.Float32Atom(dflt=nan)
CONTINUOUS_CHUNK_SHAPE = (10000, 1)

FILTERS_GZIP = tables.Filters(complevel=1)


def _open_file(filename, *args, **kwargs):
    # From pytables 3 docs:
    # [open_file] recognizes the (lowercase) names of parameters present in
    # tables/parameters.py
    if "buffer_times" not in kwargs:
        # eliminate spurious PerformanceWarning
        kwargs["buffer_times"] = inf

    return tables.open_file(str(filename), *args, **kwargs)


class _HDF5FileChromosome(Chromosome):

    default_where = "/"

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

        self.h5file = h5file
        self.filename = h5file.filename
        self.h5group = h5group

        super().__init__(name)

    def __iter__(self):
        assert self.isopen
        supercontigs = []
        for group in self.h5group:
            supercontig = _HDF5Supercontig(group)
            supercontigs.append((supercontig.start, supercontig))

        supercontigs.sort()
        for _start, supercontig in supercontigs:
            yield supercontig

    def close(self):
        if self.attrs.dirty:
            warn("Closing Chromosome with modified data. Metadata needs to"
                 " be recalculated by calling genomedata-close-data on the"
                 " Genomedata archive before re-accessing it",
                 category=GenomedataDirtyWarning)

        super().close()

    def _erase_data(self, trackname):
        """
        Currently sets the dirty bit, which can only be erased with
        genomedata-close-data GENOMEDATA
        """
        self.attrs.dirty = True
        return super()._erase_data(trackname)

    def add_trackname(self, trackname):
        _hdf5_add_trackname(self.h5file, trackname)

    def _add_track_continuous(self):
        """Add a new track

        The Genome object must have been created with
        :param mode:="r+". Behavior is undefined if this is not the case.

        Currently sets the dirty bit, which can only be erased with
        genomedata-close-data

        """
        self.attrs.dirty = True  # dirty specific to chromosome

        # Extend supercontigs by a column (or create them)
        for supercontig in self:
            supercontig_length = supercontig.end - supercontig.start
            try:
                continuous = supercontig.continuous
            except MissingContinuousData:
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
    def tracknames_continuous(self):
        """Return a list of the data track names in this Chromosome."""
        assert self.isopen
        # Tracknames are always in the root group of the h5file
        return [decode_trackname(trackname)
                for trackname in self._file_attrs.tracknames]

    @property
    def _file_attrs(self):
        assert self.isopen
        return self.h5file.root._v_attrs

    @property
    def mode(self):
        """Return the read/write properties for this chromosome as string"""
        return self.h5file.mode

    @property
    def attrs(self):
        """Return the attributes for this Chromosome.

        This may also include Genome-wide attributes if the archive
        is implemented as a directory.

        """
        assert self.isopen
        return self.h5group._v_attrs

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


class _HDF5DirectoryChromosome(_HDF5FileChromosome):
    def __init__(self, filename, mode=Chromosome.default_mode, *args,
                 **kwargs):
        r"""
        :param filename: name of the chromosome (.genomedata) file to access

        :param mode: mode of interaction with the chromosome file,
                     with ``r``: read, ``w``: write, ``a``: append,
                     ``r+``: append but force file to exist (see documentation
                     for tables.open_file().)

        :type mode: string
        :param *args: args passed on to open_file
        :param **kwargs: keyword args passed on to open_file

        """
        filepath = Path(filename).expand()
        try:
            h5file = _open_file(filepath, mode=mode, *args, **kwargs)
        except IOError:
            raise IOError("Could not find file: %r" % filename)

        name = filepath.name.rpartition(SUFFIX)[0]

        super().__init__(h5file, name=name)

    def close(self):
        super().close()
        self.h5file.close()


class _HDF5DirectoryChromosomeList(_ChromosomeList):
    def __init__(self, filepath, *args, **kwargs):
        # Folder path containing chromosome genomedata files
        self.filepath = filepath
        # Mode for accessing chromosome genomedata files
        self.args = args
        self.kwargs = kwargs

        super().__init__()

    def __iter__(self):
        # sorted so that the order is always the same
        for filepath in sorted(self.filepath.files("*" + SUFFIX)):
            # pass through __getitem__() to allow memoization
            yield self[filepath.stem]

    def __getitem__(self, name):
        # If we have opened the chromosome already
        res = super().__getitem__(name)
        # Otherwise
        if not res:
            # Return a new chromosome from from hdf5 file
            try:
                res = _HDF5DirectoryChromosome(
                    self.filepath.joinpath(name + SUFFIX),
                    *self.args, **self.kwargs)
            except IOError:
                raise KeyError("Could not find chromosome: %s" % name)

            self.open_chromosomes[name] = res

        return res

    def _format_version(self):
        res = None
        chromosomes = iter(self)

        try:
            first_chromosome = next(chromosomes)
        except StopIteration:
            return res

        # Try to get the format version from a chromsome file
        try:
            res = first_chromosome._format_version

            assert all(res == chromosome._format_version
                       for chromosome in chromosomes)
        # Otherwise assume we have no format version information
        except AttributeError:
            return res

        return res

    def tracknames_continuous(self):
        # check that all chromosomes have the same tracknames_continuous
        res = None
        for chromosome in self:
            if res is None:
                res = chromosome.tracknames_continuous
            else:
                assert res == chromosome.tracknames_continuous

        return res

    def add_trackname(self, trackname):
        for chromosome in self:
            chromosome.add_trackname(trackname)

    def create(self, name):
        res = self[name]
        res.attrs.dirty = True
        return res

    def close(self):
        # Close all the chromosomes so they know they shouldn't be read
        # NB: calls __getitem__ which checks for valid open chromosomes
        for chrom in self:
            chrom.close()

        super().close()


class _HDF5SingleFileChromosomeList(_ChromosomeList):
    def __init__(self, filepath, *args, **kwargs):
        self.h5file = _open_file(filepath, *args, **kwargs)
        self._file_attrs = self.h5file.root._v_attrs

        super().__init__()

    def __iter__(self):
        # Iterate over child group of root
        for group in self.h5file.iter_nodes("/", classname="Group"):
            groupname = group._v_name
            # Pass through to __getitem__()
            yield self[groupname]

    def __getitem__(self, name):
        # If we have opened the chromosome already
        res = super().__getitem__(name)
        # Otherwise
        if not res:
            # Return a new chromosome from from hdf5 file
            try:
                res = _HDF5FileChromosome(self.h5file, where="/" + name)
            except tables.NoSuchNodeError:
                raise KeyError("Could not find chromosome: %s" % name)
            self.open_chromosomes[name] = res

        return res

    def _format_version(self):
        # Get the Genomedata format version from the HDF5 file configurations
        try:
            return self._file_attrs.genomedata_format_version
        except AttributeError:
            return None

    def add_trackname(self, trackname):
        _hdf5_add_trackname(self.h5file, trackname)

    def tracknames_continuous(self):
        # Tracknames are stored at the root of each file, so we can
        # access them directly in this case
        return [decode_trackname(trackname)
                for trackname in self._file_attrs.tracknames]

    def create(self, name):
        self.h5file.create_group("/", name, filters=FILTERS_GZIP)

        res = self[name]
        res.attrs.dirty = True
        return res

    def close(self):
        # Close all the chromosomes so they know they shouldn't be read. Do
        # this before closing h5file in case the chromosomes need access
        # to it in closing.
        # NB: calls __getitem__ which checks for valid open chromosomes
        for chrom in self:
            chrom.close()

        self.h5file.close()
        self.isopen = False

        super().close()


def _hdf5_add_trackname(h5file, trackname):
    h5file_attr = h5file.root._v_attrs
    if "tracknames" in h5file_attr:
        tracknames = h5file_attr.tracknames
        if trackname in tracknames:
            raise ValueError("%s already has a track of name: %s" %
                             (h5file.filename, trackname))
    else:
        tracknames = array([])

    h5file_attr.tracknames = append(tracknames,
                                    trackname.encode(GENOMEDATA_ENCODING))


class _HDF5Supercontig(Supercontig):
    """A container for a segment of data in one chromosome.

    Implemented via a HDF5 Group

    """
    def __init__(self, h5group):
        """
        :param h5group: group containing the supercontig data
        :type h5group: HDF5 Group

        """
        self.h5group = h5group

    @property
    def _seq_dtype(self):
        try:
            return self.seq.atom.dtype
        except tables.NoSuchNodeError:
            return SEQ_DTYPE

    @property
    def _continuous_dtype(self):
        try:
            return self.continuous.atom.dtype
        except MissingContinuousData:
            return CONTINUOUS_DTYPE

    @property
    def continuous(self):
        """Return the underlying continuous data in this supercontig.
        To read the whole dataset into memory as a `numpy.array`, use
        continuous.read()

        :returns: `tables.EArray`

        """
        try:
            res = self.h5group.continuous
        except tables.NoSuchNodeError:
            raise MissingContinuousData

        return res

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

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


from importlib.metadata import version
import sys

from numpy import square
from path import Path

from ._hdf5 import _HDF5DirectoryChromosomeList, _HDF5SingleFileChromosomeList
from ._bigwig import _BigWigChromosomeList, is_big_wig

# Allow raising a PackageNotFoundError if somehow genomedata was not
# installed
__version__ = version("genomedata")

FORMAT_VERSION = 1


# If there are fewer than this many chomosomes, the default behavior
# is to implement the Genomedata archive as a directory. If there are
# more than this many, it will be a single file by default.
FILE_MODE_CHROMS = 100


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
        r"""Create a Genome object from a genomedata archive.

        :param filename: the root of the Genomedata object
                         hierarchy. This can either be a .genomedata
                         file that contains the entire genome or a
                         directory containing multiple chromosome files.
        :type filename: string
        :param *args: args passed on to open_file if single file or to
                      Chromosome if directory
        :param **kwargs: keyword args passed on to open_file if single file
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
        filepath = Path(filename).expand()
        while filepath.islink():
            filepath = filepath.readlinkabs()

        if not filepath.exists():
            raise IOError("Could not find Genomedata archive: %s" % filepath)

        # If it's a file we are opening
        if filepath.is_file():
            # Check if the file type is bigWig
            # NB: Could consider checking by filename extension only
            if is_big_wig(filepath):
                self._chromosomes = _BigWigChromosomeList(filepath)
            # Otherwise assume and attempt to open the Genomedata file
            else:
                self._chromosomes = _HDF5SingleFileChromosomeList(
                    filepath, *args, **kwargs)
        elif filepath.is_dir():
            # Genomedata directory
            self._chromosomes = _HDF5DirectoryChromosomeList(filepath, *args,
                                                             **kwargs)
        else:
            raise ValueError("Genomedata archive must be file or directory: %s"
                             % filepath)

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
        for chromosome in self._chromosomes:
            yield chromosome

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

        return self._chromosomes[name]

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

    def _create_chromosome(self, name, mode):
        name = "_".join(name.split())  # Remove any whitespace

        res = self._chromosomes.create(name)

        return res

    def close(self):
        """Close this Genomedata archive and any open chromosomes

        If the Genomedata archive is a directory, this closes all open
        chromosomes. If it is a single file, this closes that file.
        This should only be used if Genome is not a context manager
        (see :class:`Genome`). The behavior is undefined if this is
        called while Genome is being used as a context manager.

        """
        assert self.isopen

        self._chromosomes.close()

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
        # TODO: Add sensible error message for non HDF5 formats
        if self.format_version < 1:
            raise NotImplementedError("""Adding tracks is only supported \
for archives created with Genomedata version 1.2.0 or later.""")

        self._chromosomes.add_trackname(trackname)

        # Open continuous data regions per chromosome
        for chromosome in self:
            chromosome._add_track_continuous()

    @property
    def isopen(self):
        """Return a boolean indicating if the Genome is still open"""
        return self._isopen

    @property
    def tracknames_continuous(self):
        """Return a list of the names of all data tracks stored."""
        assert self.isopen

        return self._chromosomes.tracknames_continuous()

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

        None means there are no chromosomes in it already or there is no
        information available.
        """
        assert self.isopen

        return self._chromosomes._format_version()

    # XXX: should memoize these with an off-the-shelf decorator
    @property
    def mins(self):
        """Return the minimum value for each track.

        :returns: numpy.array

        """
        return self._chromosomes.mins

    @property
    def maxs(self):
        """Return a vector of the maximum value for each track.

        :returns: numpy.array

        """
        return self._chromosomes.maxs

    @property
    def sums(self):
        """Return a vector of the sum of the values for each track.

        :returns: numpy.array

        """
        return self._chromosomes.sums

    @property
    def sums_squares(self):
        """Return a vector of the sum of squared values for each track's data.

        :returns: numpy.array

        """
        return self._chromosomes.sums_squares

    @property
    def num_datapoints(self):
        """Return the number of datapoints in each track.

        :returns: a numpy.array vector with an entry for each track.

        """
        return self._chromosomes.num_datapoints

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


def main(argv=sys.argv[1:]):
    pass


if __name__ == "__main__":
    sys.exit(main())

import struct

import pyBigWig
from numpy import array

from ._chromosome import (Chromosome, CONTINUOUS_DTYPE, _ChromosomeList,
                          _key_to_tuple, Supercontig)


BIG_WIG_SIGNATURE = 0x888FFC26
BIG_WIG_SIGNATURE_BYTE_SIZE = 4


def is_big_wig(filename):
    """ Checks that the given filename refers to a valid bigWig file """
    with open(filename, "rb") as big_wig_file:
        signature_string = big_wig_file.read(BIG_WIG_SIGNATURE_BYTE_SIZE)

    # If we have successfully read in the number of signature bytes
    if len(signature_string) == BIG_WIG_SIGNATURE_BYTE_SIZE:
        # Check the bigWig signature

        # unpack returns a tuple regardless of length
        # the kent reference checks both little endian and big endian packing
        # of the 4 byte signature
        little_endian_signature = struct.unpack("<L", signature_string)[0]
        big_endian_signature = struct.unpack(">L", signature_string)[0]

        if (little_endian_signature == BIG_WIG_SIGNATURE or
           big_endian_signature == BIG_WIG_SIGNATURE):
            return True

    # Otherwise no signature found
    return False


# def merge_intervals(start_end_list):
#     """Takes a list of tuples consisting of start and end coordinates and
#     returns reduced list of tuples based on adjacent intervals. Assumes that
#     there are no overlapping intervals and the coordinates are sorted."""

#     first_interval = start_end_list[0]
#     merged_start = first_interval[0]
#     merged_end = first_interval[1]

#     merged_intervals = []

#     # For each start and end coordinate past the first interval
#     for start, end in start_end_list[1:]:
#         # If the start of the new interval matches the end of the current
#         # merged interval
#         if start == merged_end:
#             # Merge the interval by using the new end coord of current interval
#             merged_end = end
#         # Otherwise no further merging can take place
#         else:
#             # Put current merged coordinates into list to return
#             merged_intervals.append((merged_start, merged_end))
#             # Mark new coordinates for merging
#             merged_start = start
#             merged_end = end

#     # Append remaining merged coordinates to list
#     merged_intervals.append((merged_start, merged_end))

#     # Return the list of merged intervals
#     return merged_intervals


class _BigWigChromosomeList(_ChromosomeList):

    def __init__(self, filepath):
        # File path for bigWig file
        self.filepath = filepath
        # NB: Assumes a valid bigWig file at this point
        self.bw_file = pyBigWig.open(filepath)
        self.bw_file_header = self.bw_file.header()

        super().__init__()

    def __iter__(self):
        for chromosome_name in sorted(self.bw_file.chroms().keys()):
            yield self[chromosome_name]

    def __getitem__(self, name):
        # If we have opened the chromosome already
        res = super().__getitem__(name)
        # Otherwise
        if not res:
            # Return a new chromosome from the bigWig
            try:
                res = _BigWigChromosome(name, self.bw_file, self.filepath)
                # Cache the result
                self.open_chromosomes[name] = res
            # Otherwise raise an error
            except ValueError:
                raise KeyError("Could not find chromosome: %s" % name)

        return res

    # NB: One implicit trackname across chromosomes for stat retrieval
    @property
    def mins(self):
        return array(self.bw_file_header['minVal'])

    @property
    def maxs(self):
        return array(self.bw_file_header['maxVal'])

    @property
    def sums(self):
        return array(self.bw_file_header['sumData'])

    @property
    def sums_squares(self):
        return array(self.bw_file_header['sumSquared'])

    @property
    def num_datapoints(self):
        return array(self.bw_file_header['nBasesCovered'])

    def tracknames_continuous(self):
        # Return filepath to bigWig as implicit trackname
        return [self.filepath]

    def add_trackname(self, trackname):
        raise NotImplementedError

    def create(self, _name):
        raise NotImplementedError

    def close(self):
        # Close all the chromosomes so they know they shouldn't be read
        self.bw_file.close()

        super().close()


class _BigWigChromosome(Chromosome):

    def __init__(self, name, big_wig_file, filepath):
        self.bw_file = big_wig_file
        self.filepath = filepath

        # Check if the chromosome exists
        if name not in self.bw_file.chroms():
            raise ValueError

        # Get start and end positions from the chromosome
        self._name = name

        intervals = self.bw_file.intervals(name)  # NB: sorted list of tuples
        self._start = intervals[0][0]  # tuple: (start, end, value)
        self._end = intervals[-1][1]

        super().__init__(name)

    def __iter__(self):
        # Create a single supercontig for the entire chromosome in absense of
        # assembly information

        # NB: Since there is one supercontig per chromosome, nearly all
        # properties are shared between chromosome/supercontig. We keep the
        # underlying supercontig for interface consistency
        supercontig = _BigWigSupercontig(self)

        yield supercontig

    def __getitem__(self, key):
        """Same functionality as parent Chromosome without second key
        specifying a track or track range"""

        # Convert key to valid chromosomal indexing coordinates
        start, end = _key_to_tuple(key)

        range_length = end - start

        # If there's no length to the slice
        if range_length < 1:
            # Return an empty numpy array
            return array([], dtype=CONTINUOUS_DTYPE)

        # Otherwise return the data in range
        return self.bw_file.values(self.name, start, end, numpy=True)

    @property
    def seq(self):
        # Raise an error for a bigWig chromsome type
        raise NotImplementedError

    def close(self):
        super().close()

    # NB: Called from the public erase_data interface
    def _erase_data(self, _trackname):
        raise NotImplementedError

    def add_trackname(self, trackname):
        raise NotImplementedError

    def _add_track_continuous(self):
        raise NotImplementedError

    @property
    def tracknames_continuous(self):
        # Return filepath to bigWig as implicit trackname
        return [self.filepath]

    # NB: Undocumented public interfaces
    @property
    def mode(self):
        """Return the read/write properties for this chromosome as string"""
        pass 

    @property
    def attrs(self):
        pass

    @property
    def _format_version(self):
        return None

    @property
    def start(self):
        """Return the index of the first base in this chromosome.

        For :attr:`Genome.format_version` > 0, this will always be 0.
        For == 0, this will be the start of the first supercontig.

        """
        # if self._format_version == 0:
        #     return min(supercontig.start for supercontig in self)
        # else:
        #     return self.attrs.start
        return self._start

    @property
    def end(self):
        """Return the index past the last base in this chromosome.

        For :attr:`Genome.format_version` > 0, this will be
        the number of bases of sequence in the chromosome. For == 0,
        this will be the end of the last supercontig.

        This is the end in half-open coordinates, making slicing simple:

        >>> chromosome.seq[chromosome.start:chromosome.end]

        """
        return self._end


class _BigWigSupercontig(Supercontig):
    """A container for a segment of data in one chromosome."""

    class _Continuous():
        """An indexable-type for the continuous region.
        This type also should implement the .read() method which returns the
        whole region as a numpy array """
        def __init__(self, parent_chromosome):
            self.parent_chromosome = parent_chromosome

        def read(self):
            return self.parent_chromsome[self.parent_chromosome.start:
                                         self.parent_chromosome.end]

        def __getitem__(self, key):
            # Convert key to valid chromosomal indexing coordinates
            start, end = _key_to_tuple(key)

            # Convert key to supercontig indexing coordinates
            start += self.parent_chromosome.start
            end += self.parent_chromosome.start
            return self.parent_chromosome[start:end]

    def __init__(self, parent_chromosome):
        self.parent_chromosome = parent_chromosome
        self._continuous = self._Continuous(parent_chromosome)

    @property
    def _seq_dtype(self):
        # There should be no sequence data or type information
        raise NotImplementedError

    @property
    def _continuous_dtype(self):
        return CONTINUOUS_DTYPE

    @property
    def continuous(self):
        return self._continuous

    @property
    def attrs(self):
        raise NotImplementedError

    @property
    def name(self):
        raise NotImplementedError

    @property
    def seq(self):
        raise NotImplementedError

    @property
    def start(self):
        """Return the index of the first base in this supercontig.

        The first base is index 0.

        """
        return int(self.parent_chromosome.start)

    @property
    def end(self):
        """Return the index past the last base in this supercontig.

        This is the end in half-open coordinates, making slicing simpler:

        >>> supercontig.seq[supercontig.start:supercontig:end]

        """
        return int(self.parent_chromosome.end)
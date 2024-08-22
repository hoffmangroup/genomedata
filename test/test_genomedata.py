#!/usr/bin/env python

from __future__ import absolute_import, division, print_function
from six.moves import zip

"""
test_genomedata: DESCRIPTION
"""

__version__ = "$Revision$"

# Copyright 2010-2013 Michael M. Hoffman <mmh1@uw.edu>
import os
from tempfile import mkdtemp, mkstemp
import unittest
import warnings

from numpy import array, isnan, logical_and, logical_not, nan
from path import Path

from genomedata import Genome
from genomedata.load_genomedata import load_genomedata
from genomedata._load_data import load_data
from genomedata._close_data import close_data
from genomedata._erase_data import erase_data
from genomedata._hardmask import hardmask_data
from genomedata._open_data import open_data
from genomedata._util import (GENOMEDATA_ENCODING, GenomedataDirtyWarning,
                              OverlapWarning)

DEFAULT_TRACK_FILTER_THRESHOLD = 0.5
UNFILTERED_TRACKNAME = "zunfiltered"
UNFILTERED_TRACK_FILENAME = "unfiltered.bed"


def test_filename(filename):
    return os.path.join("data", filename)


def seq2str(seq):
    return seq.tobytes().decode(GENOMEDATA_ENCODING).lower()


def make_temp_dir():
    return mkdtemp(prefix="genomedata.test.")


class GenomedataTesterBase(unittest.TestCase):
    def setUp(self):
        # Defaults
        # XXX: adding verbosity when unittest is run with verbosity would be
        # useful
        self.verbose = False
        self.write = True
        self.mode = "dir"
        self.tracks = {"vertebrate":
                       "chr1.phyloP44way.vertebrate.short.wigFix",
                       "placental": "chr1.phyloP44way.placental.short.wigFix",
                       "primate": "chr1.phyloP44way.primate.short.wigFix"}
        # UNFILTERED_TRACKNAME: "unfiltered.bed"} # Set as last track
        # Track to be added by test_add_track
        self.new_track = ("dnase", "chr1.wgEncodeDukeDNaseSeqBase"
                          "OverlapSignalK562V2.wig")
        # File to be used to filter out the added track (default none)
        self.filter = "filter.bed"
        self.filter_threshold = DEFAULT_TRACK_FILTER_THRESHOLD

        # Potentially override defaults
        self.init()  # Call to sub-classed method

    def assertArraysEqual(self, observed, expected):
        expected = array(expected, dtype=observed.dtype)
        not_equal = (observed != expected)
        # Do some special stuff to allow testing equality between NaN's
        both_nan = logical_and(isnan(observed), isnan(observed))
        if logical_and(not_equal, logical_not(both_nan)).any():
            self.fail("%r != %r" % (observed, expected))

    def test_interface(self):
        original_num_datapoints = 0
        if self.write:
            mode = "r+"
        else:
            mode = "r"
        # catch_warnings acts as a context manager storing the original warning
        # filter and resetting it at the end. All non user warnings should
        # still be displayed
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", GenomedataDirtyWarning)
            warnings.simplefilter("ignore", OverlapWarning)
            with Genome(self.gdfilepath, mode=mode) as genome:
                original_num_datapoints = genome.num_datapoints

                self.assertTrue("chr1" in genome)
                self.assertFalse("chrZ" in genome)

                chromosome = genome["chr1"]

                # Test tracknames are as expected
                self.assertEqual(sorted(chromosome.tracknames_continuous),
                                 sorted(self.tracknames))

                # Test tracknames are consistent
                self.assertEqual(sorted(genome.tracknames_continuous),
                                 sorted(chromosome.tracknames_continuous))

                # Test chromosome attributes
                self.assertEqual(chromosome.start, 0)
                self.assertEqual(chromosome.end, 24950)

                # Test sequence inside of data range
                self.assertEqual(seq2str(chromosome.seq[0:20]),
                                 "taaccctaaccctaacccta")

                # Test sequence outside of data range
                self.assertEqual(seq2str(chromosome.seq[30000]), "n")

                # Track ordering should be: placental, primate, vertebrate
                self.assertEqual(chromosome.tracknames_continuous,
                                 self.tracknames)

                # Given track ordering, check multi-track data retrieval
                self.assertArraysEqual(chromosome[290, 0:3],
                                       [-2.297, -2.327, -2.320])

                # test multi-track data retrieval by list
                self.assertArraysEqual(chromosome[290, ["placental", "primate",
                                                        "vertebrate"]],
                                       chromosome[290, 0:3])
                self.assertArraysEqual(chromosome[290, ["placental",
                                                        "vertebrate"]],
                                       [-2.297, -2.320])
                self.assertArraysEqual(chromosome[290, [0, 2]],
                                       [-2.297, -2.320])
                self.assertArraysEqual(chromosome[290, [2, 0]],
                                       [-2.320, -2.297])

                self.assertArraysEqual(chromosome[290, array([1, 0])],
                                       [-2.327, -2.297])

                # Test filling of unassigned continuous segments
                chromosome = genome["chrY"]
                # Get first supercontig
                for supercontig in chromosome:
                    break
                self.assertArraysEqual(supercontig.continuous[0, 2], nan)

                # If we are testing writing to archives
                if self.write:
                    # Test writing scalars to various indexing methods
                    chromosome = genome["chr1"]
                    # Test writing scalar to multiple tracks
                    chromosome[290] = 100.0
                    # Test writing scalar to tracks by named list
                    chromosome[291,
                               ["placental", "primate", "vertebrate"]] = 101.0
                    # Test writing scalar to select tracks by named list
                    chromosome[292, ["placental", "vertebrate"]] = 102.0
                    # Test writing scalar to tracks by index
                    chromosome[293, [0, 2]] = 103.0
                    chromosome[294, [2, 0]] = 104.0

                    # Test writing arrays to various indexing methods
                    # Test writing an array to a single index
                    chromosome[295] = [105.0, 106.0, 107.0]
                    # Test writing a subarray to a index subset
                    chromosome[296,
                               ["placental", "vertebrate"]] = [108.0, 109.0]

                    # Test removing datapoints by writing NaN
                    chromosome[297, ["primate"]] = nan

                    # Test writing around supercontig boundaries
                    # <Supercontig 'supercontig_0', [0:24950]>
                    # Test writing outside a supercontig
                    try:
                        chromosome[300000] = 110.0
                    except ValueError:
                        pass  # we expect a value error here

                    # Test writing overlap across supercontig to no supercontig
                    try:
                        chromosome[24900:30000] = 111.0
                    except ValueError:
                        pass  # we expect a value error here

        # Check write output after closing if testing writes
        if self.write:
            # Close with newly written data
            close_data(self.gdfilepath, verbose=self.verbose)
            # Read data and verify new data and parameters
            with Genome(self.gdfilepath) as genome:
                chromosome = genome["chr1"]

                self.assertArraysEqual(chromosome[290],
                                       [100.0, 100.0, 100.0])
                self.assertArraysEqual(chromosome[291],
                                       [101.0, 101.0, 101.0])
                # L14 in primate wigFix
                self.assertArraysEqual(chromosome[292],
                                       [102.0, 0.371, 102.0])
                self.assertArraysEqual(chromosome[293],
                                       [103.0, 0.372, 103.0])
                self.assertArraysEqual(chromosome[294],
                                       [104.0, 0.372, 104.0])
                self.assertArraysEqual(chromosome[295],
                                       [105.0, 106.0, 107.0])
                self.assertArraysEqual(chromosome[296],
                                       [108.0, -2.327, 109.0])
                # Check if one datapoint was successfully removed
                self.assertArraysEqual(original_num_datapoints,
                                       [genome.num_datapoints[0],
                                        genome.num_datapoints[1] + 1,
                                        genome.num_datapoints[2]])

    def test_repr_str(self):
        genome = Genome(self.gdfilepath, mode="r")
        self.assertEqual(repr(genome), "Genome('%s', **{'mode': 'r'})" %
                         self.gdfilepath)
        chr = genome["chr1"]
        if self.mode == "dir":
            self.assertEqual(repr(chr),
                             "<Chromosome 'chr1', file='%s/chr1.genomedata'>" %
                             self.gdfilepath)
            self.assertEqual(str(chr), "chr1")
        elif self.mode == "file":
            self.assertEqual(repr(chr),
                             "<Chromosome 'chr1', file='%s'>" %
                             self.gdfilepath)
            self.assertEqual(str(chr), "chr1")

        genome.close()

    def test_no_context(self):
        genome = Genome(self.gdfilepath)
        chr1 = genome["chr1"]
        genome.tracknames_continuous  # test access
        chr1[100:1000]  # test access: at one point segfaulted
        chr2 = genome["chrY"]
        chr2.close()  # Make sure manual close doesn't break it
        self.assertTrue(chr1.isopen)
        self.assertFalse(chr2.isopen)
        genome.close()
        self.assertFalse(chr1.isopen)
        self.assertRaises(Exception, next, iter(chr1))

    def test_open_chromosomes(self):
        genome = Genome(self.gdfilepath)
        with genome:
            chr1 = genome["chr1"]
            chr2 = genome["chr1"]  # Memoized
            self.assertEqual(chr1, chr2)
            genome["chrY"]
            self.assertEqual(len(genome._chromosomes.open_chromosomes), 2)

        self.assertEqual(genome._chromosomes.open_chromosomes, {})

    def set_gdfilepath(self, filename):
        self.gdfilepath = Path(filename).expand()


class GenomedataTester(GenomedataTesterBase):
    def setUp(self):
        GenomedataTesterBase.setUp(self)

        # Create Genomedata collection from test files
        seqs = ["chr1.short.fa", "chrY.short.fa.gz"]
        # Placental includes data for chr1 and chrY
        if self.mode == "dir":
            gdfilename = make_temp_dir()

        elif self.mode == "file":
            tempfile, gdfilename = mkstemp(prefix="genomedata", text=True)
            os.close(tempfile)
            os.remove(gdfilename)  # Allow load_genomedata to create it
        else:
            self.fail("Unrecognized mode: %s" % self.mode)

        self.set_gdfilepath(gdfilename)
        self.tracknames = sorted(self.tracks.keys())

        # Get resource paths instead of filenames
        seqfiles = [test_filename(file) for file in seqs]
        trackfiles = [test_filename(self.tracks[trackname])
                      for trackname in self.tracknames]
        self.seqfiles = seqfiles
        self.trackfiles = trackfiles

        tracktuples = [(trackname, trackfile) for trackname, trackfile in
                       zip(self.tracknames, self.trackfiles)]

        load_genomedata(self.gdfilepath, tracktuples, seqfiles,
                        verbose=self.verbose, mode=self.mode)

        if self.mode == "dir":
            self.chroms = [val.split(".")[0] for val in seqs]
            for chrom in self.chroms:
                filename = os.extsep.join([chrom, "genomedata"])
                filepath = self.gdfilepath.joinpath(filename)
                self.assertTrue(filepath.is_file(),
                                "Chromosome file was not found: %s" % filepath)
        elif self.mode == "file":
            self.assertTrue(self.gdfilepath.is_file(),
                            "Genomedata archive was not created: %r" %
                            self.gdfilepath)
        else:
            self.fail("Unrecognized mode: %s" % self.mode)

    def tearDown(self):
        return
        if self.mode == "dir":
            self.gdfilepath.rmtree()
        elif self.mode == "file":
            self.gdfilepath.remove()
        else:
            self.fail("Unrecognized mode: %s" % self.mode)

    def test_add_track(self):
        new_track_name, new_track_file = self.new_track

        # Open new track
        genome = Genome(self.gdfilepath, mode="r+")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", GenomedataDirtyWarning)
            with genome:
                genome.add_track_continuous(new_track_name)

        # Load data for new track
        load_data(self.gdfilepath, new_track_name,
                  test_filename(new_track_file), verbose=self.verbose)

        # Close data with new track
        close_data(self.gdfilepath, verbose=self.verbose)

        # Make sure addition was successful
        genome = Genome(self.gdfilepath)
        with genome:
            # Track ordering should now end with dnase
            self.assertEqual(genome.tracknames_continuous,
                             self.tracknames + [new_track_name])

            # Given track ordering, check single track data retrieval
            self.assertArraysEqual(genome["chr1"][155:168, new_track_name],
                                   [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0])

    def test_filter_track(self):
        # Add filter track
        open_data(self.gdfilepath, [UNFILTERED_TRACKNAME],
                  verbose=self.verbose)
        load_data(self.gdfilepath, UNFILTERED_TRACKNAME,
                  test_filename(UNFILTERED_TRACK_FILENAME),
                  verbose=self.verbose)
        close_data(self.gdfilepath, verbose=self.verbose)

        # Perform filtering on data
        hardmask_data(self.gdfilepath, test_filename(self.filter),
                      [UNFILTERED_TRACKNAME],
                      lambda x: x < self.filter_threshold,
                      verbose=self.verbose)

        # Make sure filtering was successful
        genome = Genome(self.gdfilepath)
        with genome:
            self.assertArraysEqual(genome["chr1"][0:4,
                                   UNFILTERED_TRACKNAME],
                                   [nan, nan, nan, nan])
            self.assertArraysEqual(genome["chr1"][128:132,
                                   UNFILTERED_TRACKNAME],
                                   [nan, nan, 0.5, 0.5])
            self.assertArraysEqual(genome["chr1"][168:172,
                                   UNFILTERED_TRACKNAME],
                                   [0.9, 0.9, nan, nan])
            self.assertArraysEqual(genome["chr1"][206:210,
                                   UNFILTERED_TRACKNAME],
                                   [nan, nan, nan, nan])

    def test_delete_tracks(self):
        # Test ability to delete a track
        trackname = "primate"
        old_entry = (290, -2.327)

        # Test value before deleting track
        warnings.simplefilter("ignore")
        with Genome(self.gdfilepath, "r+") as genome:
            chromosome = genome["chr1"]
            self.assertArraysEqual(chromosome[old_entry[0], trackname],
                                   old_entry[1])
            chromosome._erase_data(trackname)

        warnings.resetwarnings()
        # Re-close data
        close_data(self.gdfilepath, verbose=self.verbose)

        # Test value after deleting track
        with Genome(self.gdfilepath) as genome:
            chromosome = genome["chr1"]
            self.assertArraysEqual(chromosome[old_entry[0], trackname],
                                   old_entry[1])

    def test_replace_track(self):
        # Test ability to delete and replace a track

        old_trackname = "primate"
        old_entry = (290, -2.327)
        new_trackname = "placental"
        new_entry = (290, -2.297)

        # Test value before deleting track
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", GenomedataDirtyWarning)
            with Genome(self.gdfilepath) as genome:
                chromosome = genome["chr1"]
                self.assertArraysEqual(chromosome[old_entry[0], old_trackname],
                                       old_entry[1])

            # Remove track
            erase_data(self.gdfilepath, old_trackname, verbose=self.verbose)

        # Now replace it with the data from a different track
        track_index = self.tracknames.index(new_trackname)
        datafile = self.trackfiles[track_index]
        load_data(self.gdfilepath, new_trackname, datafile,
                  verbose=self.verbose)

        # Re-close data
        close_data(self.gdfilepath, verbose=self.verbose)

        # Now test that the data matches the new track data
        with Genome(self.gdfilepath) as genome:
            chromosome = genome["chr1"]
            self.assertArraysEqual(chromosome[new_entry[0], new_trackname],
                                   new_entry[1])


class GenomedataGivenDataTester(GenomedataTesterBase):
    """
    test a given Genomedata file
    """
    def setUp(self):
        GenomedataTesterBase.setUp(self)

        self.tracknames = ["placental", "primate", "vertebrate"]


# XXX: why isn't this a sublcass of GenomeDataTesterBase?
class GenomedataNoDataTester(unittest.TestCase):
    def setUp(self):
        # Defaults
        self.verbose = False
        self.mode = "dir"
        # Track to be added by test_add_track
        self.new_track = ("primate", "chr1.phyloP44way.primate.short.wigFix")

        # Potentially override defaults
        self.init()  # Call to sub-classed method

        # Create Genomedata collection from test files
        seqs = ["chr1.short.fa", "chrY.short.fa.gz"]
        # Placental includes data for chr1 and chrY
        if self.mode == "dir":
            gdfilename = make_temp_dir()

        elif self.mode == "file":
            tempfile, gdfilename = mkstemp(prefix="genomedata")
            os.close(tempfile)
            os.remove(gdfilename)  # Allow load_genomedata to create it
        else:
            self.fail("Unrecognized mode: %s" % self.mode)

        self.gdfilepath = Path(gdfilename).expand()

        # Get resource paths instead of filenames
        seqfiles = [test_filename(file) for file in seqs]
        self.seqfiles = seqfiles

        load_genomedata(self.gdfilepath, None, seqfiles,
                        verbose=self.verbose, mode=self.mode)

        if self.mode == "dir":
            self.chroms = [val.split(".")[0] for val in seqs]
            for chrom in self.chroms:
                filename = os.extsep.join([chrom, "genomedata"])
                filepath = self.gdfilepath.joinpath(filename)
                self.assertTrue(filepath.is_file(),
                                "Chromosome file was not found: %s" % filepath)
        elif self.mode == "file":
            self.assertTrue(self.gdfilepath.is_file(),
                            "Genomedata archive was not created: %r" %
                            self.gdfilepath)
        else:
            self.fail("Unrecognized mode: %s" % self.mode)

    def tearDown(self):
        if self.mode == "dir":
            self.gdfilepath.rmtree()
        elif self.mode == "file":
            self.gdfilepath.remove()
        else:
            self.fail("Unrecognized mode: %s" % self.mode)

    def assertArraysEqual(self, observed, expected):
        expected = array(expected, dtype=observed.dtype)
        not_equal = (observed != expected)
        # Do some special stuff to allow testing equality between NaN's
        both_nan = logical_and(isnan(observed), isnan(observed))
        if logical_and(not_equal, logical_not(both_nan)).any():
            self.fail("%r != %r" % (observed, expected))

    def test_add_track(self):
        new_track_name, new_track_file = self.new_track

        # Open new track
        genome = Genome(self.gdfilepath, mode="r+")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", GenomedataDirtyWarning)
            with genome:
                self.assertEqual(genome.num_tracks_continuous, 0)
                genome.add_track_continuous(new_track_name)

        # Load data for new track
        load_data(self.gdfilepath, new_track_name,
                  test_filename(new_track_file), verbose=self.verbose)

        # Close data with new track
        close_data(self.gdfilepath, verbose=self.verbose)

        # Make sure addition was successful
        genome = Genome(self.gdfilepath)
        with genome:
            # Track ordering should now end with dnase
            self.assertEqual(genome.tracknames_continuous, [new_track_name])

            # Given track ordering, check single track data retrieval
            self.assertArraysEqual(genome["chr1"][305:310, new_track_name],
                                   [-2.65300012, 0.37200001, 0.37200001,
                                    0.37200001, 0.37099999])

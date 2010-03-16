#!/usr/bin/env python
from __future__ import division, with_statement

"""
DESCRIPTION
"""

__version__ = "$Revision: $"

# Copyright 2009 Orion Buske <stasis@uw.edu>

import inspect
import os
import sys
import unittest
import warnings

from numpy import array, isnan, logical_and, logical_not, NAN
from path import path
from tempfile import mkdtemp, mkstemp

from genomedata import Genome
from genomedata.load_genomedata import load_genomedata
from genomedata._load_data import load_data
from genomedata._close_data import close_data
from genomedata._erase_data import erase_data

test_filename = lambda filename: os.path.join("data", filename)

def seq2str(seq):
    return seq.tostring().lower()

class GenomedataTester(unittest.TestCase):
    def init(self):
        pass

    def setUp(self):
        self.verbose = False
        self.mode = "dir"

        self.init()  # Call to sub-classed method

        # Create Genomedata collection from test files
        seqs = ["chr1.short.fa", "chrY.short.fa"]
        tracks = {"vertebrate": "chr1.phyloP44way.vertebrate.short.wigFix",
                  "placental": "chr1.phyloP44way.placental.short.wigFix",
                  "primate": "chr1.phyloP44way.primate.short.wigFix",
                  "dnase": "chr1.wgEncodeDukeDNaseSeqBaseOverlap" \
                      "SignalK562V2.wig"}
        if self.mode == "dir":
            gdfilename = mkdtemp(prefix="genomedata")

        elif self.mode == "file":
            tempfile, gdfilename = mkstemp(prefix="genomedata")
            os.close(tempfile)
            os.remove(gdfilename)  # Allow load_genomedata to create it
        else:
            self.fail("Unrecognized mode: %s" % self.mode)

        self.gdfilepath = path(gdfilename).expand()

        self.tracknames = sorted(tracks.keys())

        # Get resource paths instead of filenames
        seqfiles = [test_filename(file) for file in seqs]
        trackfiles = [test_filename(tracks[trackname])
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
                self.assertTrue(filepath.isfile(),
                                "Chromosome file was not found: %s" % filepath)
        elif self.mode == "file":
            self.assertTrue(self.gdfilepath.isfile(),
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

    def test_interface(self):
        with Genome(self.gdfilepath) as genome:
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
            warnings.simplefilter("ignore")
            self.assertEqual(seq2str(chromosome.seq[30000]), "n")
            warnings.resetwarnings()

            # Track ordering should be: dnase, placental, primate, vertebrate
            self.assertEqual(chromosome.tracknames_continuous, self.tracknames)

            # Given track ordering, check single track data retrieval
            self.assertArraysEqual(chromosome[155:168, "dnase"],
                                   [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0])

            # Given track ordering, check multi-track data retrieval
            self.assertArraysEqual(chromosome[290, 1:4],
                                   [-2.297, -2.327, -2.320])

            # Test filling of unassigned continuous segments
            chromosome = genome["chrY"]
            # Get first supercontig
            for supercontig in chromosome:
                break
            self.assertArraysEqual(supercontig.continuous[0, 3], NAN)

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
        tracknames = genome.tracknames_continuous
        data = chr1[100:1000]  # Used to segfault
        chr2 = genome["chrY"]
        chr2.close()  # Make sure manual close doesn't break it
        self.assertTrue(chr1.isopen)
        self.assertFalse(chr2.isopen)
        genome.close()
        self.assertFalse(chr1.isopen)
        self.assertRaises(Exception, iter(chr1).next)

    def test_open_chromosomes(self):
        genome = Genome(self.gdfilepath)
        with genome:
            chr1 = genome["chr1"]
            chr2 = genome["chr1"]  # Memoized
            self.assertEqual(chr1, chr2)
            chr3 = genome["chrY"]
            self.assertEqual(len(genome.open_chromosomes), 2)

        self.assertEqual(genome.open_chromosomes, {})

    def test_delete_tracks(self):
        # Test ability to delete a track
        trackname = "primate"
        old_entry = (290, -2.327)
        new_entry = (290, NAN)

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
        with Genome(self.gdfilepath) as genome:
            chromosome = genome["chr1"]
            self.assertArraysEqual(chromosome[old_entry[0], old_trackname],
                                   old_entry[1])

        # Remove track
        erase_data(self.gdfilepath, old_trackname,
                      verbose=self.verbose)

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

class TestGenomedataDir(GenomedataTester):
    def init(self):
        self.mode = "dir"

class TestGenomedataFile(GenomedataTester):
    def init(self):
        self.mode = "file"

def suite():
    def is_test_class(member):
        name, value = member
        return inspect.isclass(value) and name.startswith("Test")

    classes = []
    members = inspect.getmembers(sys.modules[__name__])
    for name, value in members:
        if inspect.isclass(value) and name.startswith("Test"):
            classes.append(value)

    tests = map(unittest.TestLoader().loadTestsFromTestCase, classes)
    return unittest.TestSuite(tests)

if __name__ == "__main__":
    unittest.main()

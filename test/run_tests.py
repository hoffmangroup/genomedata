#!/usr/bin/env python
from __future__ import division, with_statement

"""
DESCRIPTION
"""

__version__ = "$Revision: $"

# Copyright 2009 Orion Buske <stasis@uw.edu>
# Copyright 2010 Michael M. Hoffman <mmh1@uw.edu>

import inspect
import os
import sys
from tempfile import mkdtemp, mkstemp
import unittest

from numpy import array, isnan, logical_and, logical_not
from path import path

from genomedata import Genome
from genomedata.load_genomedata import load_genomedata
from genomedata._load_data import load_data
from genomedata._close_data import close_data

import test_genomedata
from test_genomedata import test_filename

class TestGenomedataDir(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "dir"

class TestGenomedataFile(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"

class TestParseWigVar(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"
        self.new_track = ("dnase", "wgEncodeDukeDNaseSeqBase"
                          "OverlapSignalK562V2.wigVar")

class TestParseBed(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"
        self.new_track = ("dnase", "wgEncodeDukeDNaseSeqBase"
                          "OverlapSignalK562V2.bed")

class TestParseBedGraph(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"
        self.new_track = ("dnase", "wgEncodeDukeDNaseSeqBase"
                          "OverlapSignalK562V2.bedGraph")

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
            gdfilename = mkdtemp(prefix="genomedata")

        elif self.mode == "file":
            tempfile, gdfilename = mkstemp(prefix="genomedata")
            os.close(tempfile)
            os.remove(gdfilename)  # Allow load_genomedata to create it
        else:
            self.fail("Unrecognized mode: %s" % self.mode)

        self.gdfilepath = path(gdfilename).expand()

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

    def test_add_track(self):
        new_track_name, new_track_file = self.new_track

        # Open new track
        genome = Genome(self.gdfilepath, mode="r+")
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

class TestNoDataGenomedataDir(GenomedataNoDataTester):
    def init(self):
        self.mode = "dir"

class TestNoDataGenomedataFile(GenomedataNoDataTester):
    def init(self):
        self.mode = "file"

def is_test_class(member):
    name, value = member

    return inspect.isclass(value) and name.startswith("Test")

def suite():
    """
    return only classes that start with "Test"
    """
    classes = []
    members = inspect.getmembers(sys.modules[__name__])
    for name, value in members:
        if inspect.isclass(value) and name.startswith("Test"):
            classes.append(value)

    tests = map(unittest.TestLoader().loadTestsFromTestCase, classes)
    return unittest.TestSuite(tests)

def main():
    cur_dir = os.path.dirname(__file__)
    if cur_dir:
        os.chdir(cur_dir)

    unittest.main()

    # unittest.TextTestRunner(verbosity=2).run(suite())

if __name__ == "__main__":
    main()

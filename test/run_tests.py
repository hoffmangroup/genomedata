#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import os
from os import chdir, remove
import subprocess
import unittest

from path import path

from genomedata.load_genomedata import load_genomedata
from genomedata import Genome

import test_genomedata

"""
DESCRIPTION
"""

__version__ = "$Revision: $"

# Copyright 2010-2011, 2013 Michael M. Hoffman <mmh1@uw.edu>
# Copyright 2009 Orion Buske <stasis@uw.edu>


def isBigWigToBedGraphInstalled():
    """Checks if the UCSC tool bigWigToBedGraph is installed"""
    # NB:
    # bigWigToBedGraph outputs to stderr always
    # bigWigToBedGraph will always return a non-zero exit unless it actually
    # converts a bigWig file. There is no "-version" option for a zero exit.

    # Try to run bigWigToBedGraph
    try:
        subprocess.check_output(["bigWigToBedGraph"], stderr=subprocess.STDOUT)
    # If the program wasn't found
    except OSError:
        # return false
        return False
    # Otherwise if the program was found and returned non-zero exit status
    except subprocess.CalledProcessError:
        # return true
        return True


class TestGenomedataDir(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "dir"


class TestGenomedataFile(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"


class TestFilterWigFixed(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"
        self.filter = "filter_fixed.wig"


class TestParseWigVar(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"
        self.new_track = ("dnase", "wgEncodeDukeDNaseSeqBase"
                          "OverlapSignalK562V2.wigVar")
        self.filter = "filter_variable.wig.gz"


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


@unittest.skipIf(isBigWigToBedGraphInstalled() is False,
                 "bigWigToBedGraph not installed")
class TestParseBigWig(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"
        self.new_track = ("bigW", "ENCFF324BPE.bigWig")


class TestParseWigVarDOS(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"
        self.new_track = ("dnase", "wgEncodeDukeDNaseSeqBase"
                          "OverlapSignalK562V2.dos.wigVar")


class TestParseBedDOS(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"
        self.new_track = ("dnase", "wgEncodeDukeDNaseSeqBase"
                          "OverlapSignalK562V2.dos.bed")


class TestParseBedGraphDOS(test_genomedata.GenomedataTester):
    def init(self):
        self.mode = "file"
        self.new_track = ("dnase", "wgEncodeDukeDNaseSeqBase"
                          "OverlapSignalK562V2.dos.bedGraph")


class TestGivenDataV0(test_genomedata.GenomedataGivenDataTester):
    def init(self):
        self.mode = "file"
        self.write = False
        self.set_gdfilepath("data/v0.genomedata")


class TestGivenDataV1(test_genomedata.GenomedataGivenDataTester):
    def init(self):
        self.mode = "file"
        self.write = False
        self.set_gdfilepath("data/v1.genomedata")


class TestNoDataGenomedataDir(test_genomedata.GenomedataNoDataTester):
    def init(self):
        self.mode = "dir"


class TestNoDataGenomedataFile(test_genomedata.GenomedataNoDataTester):
    def init(self):
        self.mode = "file"


class TestChunks(unittest.TestCase):

    def data_path(self, filename):
        return os.path.join("data", filename)

    def setUp(self):
        self.genomedata_name = self.data_path("chunk_test.genomedata")
        self.track_tuples = [("track{}".format(i),
                             self.data_path(
                                 "chunk_test_track{}.bed".format(i)))
                             for i in range(1, 3)]
        self.sequence_type = "agp"
        self.sequences = [self.data_path("test_chunks.agp")]
        self.mode = "file"
        self.verbose = False

        load_genomedata(self.genomedata_name, self.track_tuples,
                        self.sequences, self.mode, self.sequence_type,
                        verbose=self.verbose)

    def test_chunk_positions(self):
        # Get all chunk starts and ends for each supercontig
        supercontig_chunks = []
        with Genome(self.genomedata_name) as genome:
            # Assume only chr1 is defined
            chromosome = genome["chr1"]
            # For every supercontig (in chr1)
            for supercontig in chromosome:
                supercontig_attrs = supercontig.attrs
                # Append a tuple of a list of chunk starts and ends
                # NB: There should only be one chunk listed per supercontig for
                # this test
                supercontig_chunks.append(
                    (supercontig_attrs.chunk_starts.tolist(),
                     supercontig_attrs.chunk_ends.tolist())
                )

        expected_chunks = [
            # Test first chunk extends from 25 to end of first supercontig
            ([50], [100]),
            # Test second chunk extends from 125 to 150
            ([125], [150]),
            # Test third supercontig contains no chunks
            ([], []),
            # Test fourth supercontig contains a chunk spanning entire
            # supercontig
            ([0], [100]),
            # Test last chunk start at fifth supercontig and ends early
            ([0], [25]),
        ]
        self.assertEqual(expected_chunks, supercontig_chunks)

    def tearDown(self):
        remove(self.genomedata_name)


def main():
    dirpath = path(__file__).dirname()
    if dirpath:
        chdir(dirpath)

    unittest.main()

if __name__ == "__main__":
    main()

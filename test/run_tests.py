#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import os
from os import chdir, remove
import subprocess
import unittest

import numpy as np
from path import Path

from genomedata._load_seq import load_seq
from genomedata._close_data import close_data
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


def test_data_path(filename):
    return os.path.join("data", filename)


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

    def setUp(self):
        self.genomedata_name = test_data_path("chunk_test.genomedata")
        self.track_tuples = [("track{}".format(i),
                             test_data_path(
                                 "chunk_test_track{}.bed".format(i)))
                             for i in range(1, 3)]
        self.sequence_type = "agp"
        self.sequences = [test_data_path("test_chunks.agp")]
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


class TestAGPMerge(unittest.TestCase):

    def setUp(self):
        self.genomedata_name = test_data_path("mergeagp.genomedata")
        self.verbose = False
        self.agp_filenames = [test_data_path("adjacent.agp.gz")]
        self.mode = "file"

        load_seq(self.genomedata_name, self.agp_filenames, self.verbose,
                 self.mode, seqfile_type="agp")
        close_data(self.genomedata_name, self.verbose)

    def test_merged_agp_regions(self):
        with Genome(self.genomedata_name) as genome:
            # Assume only chr10 is defined
            chromosome = genome["chr10"]
            # The following supercontigs should span the following coordinates:
            expected_supercontig_coordinates = [
                (39239118, 39254773),
                (39254793, 39338430),
                (39338450, 39341685),
            ]
            supercontig_coordinates = []
            # For every supercontig (in chr10)
            for supercontig in chromosome:
                # Create a list of start/end coordinates
                supercontig_attrs = supercontig.attrs
                supercontig_coordinates.append(
                    (supercontig_attrs.start,
                     supercontig_attrs.end)
                )

        self.assertEqual(expected_supercontig_coordinates,
                         supercontig_coordinates)

    def tearDown(self):
        remove(self.genomedata_name)


class TestChromosomeNameTranslation(unittest.TestCase):

    def setUp(self):
        self.genomedata_name = test_data_path("mapped_names.genomedata")
        self.verbose = False
        self.agp_filenames = [test_data_path("chrY.agp.gz"),
                              test_data_path("chr1.agp.gz")]
        self.mode = "file"

        with open(test_data_path("assembly_report.txt"), "r") as \
                assembly_report_file:
            load_seq(self.genomedata_name, self.agp_filenames, self.verbose,
                     self.mode, seqfile_type="agp",
                     assembly_report_file=assembly_report_file,
                     chromosome_name_style="UCSC-style-name")

        close_data(self.genomedata_name, self.verbose)

    def test_translated_chromsome_names(self):
        chromosome_names = []

        with Genome(self.genomedata_name) as genome:
            # Check to see if chrY from Genbank and chr1 from RefSeq were
            # translated
            chromosome_names = [chromosome.name for chromosome in genome]

        self.assertIn("chr1", chromosome_names)
        self.assertIn("chrY", chromosome_names)

    def tearDown(self):
        remove(self.genomedata_name)


class GenomedataBigWigTester(unittest.TestCase):
    def setUp(self):
        self.bigwig_name = test_data_path("ENCFF324BPE.bigWig")

    def test_interface(self):
        with Genome(self.bigwig_name) as genome:

            self.assertEqual(repr(genome),
                             "Genome('{}')".format(self.bigwig_name))

            # Test sorted order of chromosomes
            chr_names = [chromosome.name for chromosome in genome]
            self.assertListEqual(chr_names, sorted(chr_names))

            # Test chromosome attrbutes
            self.assertTrue("chr1" in genome)
            self.assertFalse("chrY" in genome)

            # Test for no format version
            self.assertTrue(genome.format_version is None)

            # Test genome is open
            self.assertTrue(genome.isopen)

            # Test failure for attempting to erase data
            with self.assertRaises(NotImplementedError):
                genome.erase_data(self.bigwig_name)

            # Test for trackname "list"
            self.assertListEqual(genome.tracknames_continuous,
                                 [self.bigwig_name])

            self.assertEqual(genome.num_datapoints, [3088137])
            self.assertEqual(genome.mins, [0])
            self.assertEqual(genome.maxs, [93473])
            self.assertEqual(genome.sums, [32720078])
            self.assertEqual(genome.sums_squares, [1280737190372])
            self.assertAlmostEqual(genome.vars[0], 414615.8372, places=4)

            # Test chromosome retrieval
            chr1 = genome["chr1"]  # memoization

            self.assertEqual(repr(chr1),
                             "<Chromosome 'chr1', file='{}'>".format(
                             self.bigwig_name))
            self.assertEqual(str(chr1), "chr1")

            # Chromosome interface testing
            # chr1	569800	569801	0.01108
            # chr1	569801	569829	0.02215

            self.assertAlmostEqualList(
                chr1[569798:569802, "foo.track"],  # should ignore tracknames
                [np.nan, np.nan, 0.01108, 0.02215],
            )

            # Single element should result in a single value
            # Expected shape output
            # c[i, j] == data
            # c[i] == [data]
            # c[i:j, k] = [data]
            # c[i:j, [k]] = [[data]]
            single_val = chr1[569800, 0]
            self.assertAlmostEqual(single_val, 0.01108, places=5)
            self.assertFalse(isinstance(single_val, np.ndarray))

            single_val_array = chr1[569800]
            self.assertTrue(isinstance(single_val_array, np.ndarray))
            self.assertTrue(single_val_array.shape == (1,))

            single_val_array = chr1[569800:569801, 0]
            self.assertTrue(isinstance(single_val_array, np.ndarray))
            self.assertTrue(single_val_array.shape == (1,))

            # Track index shape should shape data output
            self.assertAlmostEqualList(
                np.array([[np.nan], [0.01108], [0.02215]]),
                chr1[569799:569802, [0]]
            )

            # Test the trackname is implictly the filename
            self.assertEqual(self.bigwig_name, chr1.tracknames_continuous[0])

            # Test chromosome attributes
            self.assertEqual(chr1.start, 569800)
            self.assertEqual(chr1.end, 249200695)

            # Test sequence data should error
            with self.assertRaises(NotImplementedError):
                _ = chr1.seq[0:20]

            # Test entire Supercontig interface
            # Get first (and only) supercontig
            for supercontig in chr1:
                break

            # Slice indexing
            self.assertAlmostEqualList(supercontig.continuous[1:2],
                                       [0.02215, 0.02215])
            # Single indexing
            self.assertAlmostEqualList(supercontig.continuous[0], [0.01108])

    def assertAlmostEqualList(self, first, second, places=1):
        """ Compares floating point numbers in a list if close enough.
        Does not apply to nans """

        self.assertTrue(
            np.allclose(first, second, equal_nan=True)
        )


def main():
    dirpath = Path(__file__).dirname()
    if dirpath:
        chdir(dirpath)

    unittest.main()


if __name__ == "__main__":
    main()

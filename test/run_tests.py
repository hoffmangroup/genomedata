#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

from os import chdir
import subprocess
import unittest

from path import path

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
        self.set_gdfilepath("data/v0.genomedata")


class TestGivenDataV1(test_genomedata.GenomedataGivenDataTester):
    def init(self):
        self.mode = "file"
        self.set_gdfilepath("data/v1.genomedata")


class TestNoDataGenomedataDir(test_genomedata.GenomedataNoDataTester):
    def init(self):
        self.mode = "dir"


class TestNoDataGenomedataFile(test_genomedata.GenomedataNoDataTester):
    def init(self):
        self.mode = "file"


def main():
    dirpath = path(__file__).dirname()
    if dirpath:
        chdir(dirpath)

    unittest.main()

if __name__ == "__main__":
    main()

#!/usr/bin/env python
from __future__ import division, with_statement

"""
DESCRIPTION
"""

__version__ = "$Revision: $"

# Copyright 2010-2011 Michael M. Hoffman <mmh1@uw.edu>
# Copyright 2009 Orion Buske <stasis@uw.edu>

from os import chdir
import unittest

from path import path

import test_genomedata

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

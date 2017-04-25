#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

# Copyright 2008-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

from contextlib import closing
from gzip import open as _gzip_open
from os import extsep
from pybedtools import BedTool
from pybedtools.contrib.bigwig import bigwig_to_bedgraph
import struct
import sys

from numpy import array, empty
from tables import Filters

BIG_WIG_FIELD_COUNT = 4
BIG_WIG_SIGNATURE = 0x888FFC26
BIG_WIG_SIGNATURE_BYTE_SIZE = 4

FILTERS_GZIP = Filters(complevel=1)

EXT_GZ = "gz"
SUFFIX_GZ = extsep + EXT_GZ

def die(msg="Unexpected error."):
    print(msg, file=sys.stderr)
    sys.exit(1)

class LightIterator(object):
    def __init__(self, handle):
        self._handle = handle
        self._defline = None

    def __iter__(self):
        return self

    def next(self):
        lines = []
        defline_old = self._defline

        for line in self._handle:
            if not line:
                if not defline_old and not lines:
                    raise StopIteration
                if defline_old:
                    self._defline = None
                    break
            elif line.startswith(">"):
                self._defline = line[1:].rstrip()
                if defline_old or lines:
                    break
                else:
                    defline_old = self._defline
            else:
                lines.append(line.rstrip())

        if not lines:
            raise StopIteration

        if defline_old is None:
            raise ValueError("no definition line found at next position in %r" % self._handle)

        return defline_old, ''.join(lines)

# XXX: suggest as default
def fill_array(scalar, shape, dtype=None, *args, **kwargs):
    if dtype is None:
        dtype = array(scalar).dtype

    res = empty(shape, dtype, *args, **kwargs)
    res.fill(scalar)

    return res

# XXX: suggest as default
def gzip_open(*args, **kwargs):
    return closing(_gzip_open(*args, **kwargs))

def maybe_gzip_open(filename, *args, **kwargs):
    if filename.endswith(SUFFIX_GZ):
        return gzip_open(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)


def init_num_obs(num_obs, continuous):
    curr_num_obs = continuous.shape[1]
    assert num_obs is None or num_obs == curr_num_obs

    return curr_num_obs

def new_extrema(func, data, extrema):
    curr_extrema = func(data, 0)

    return func([extrema, curr_extrema], 0)


def get_float_from_filter_score(score):
    try:
        result = float(score)
    except ValueError:
        print("Found ", score, "in filter value, could not convert to a "
              "number")
        raise
    return result
        

def score_filter(bed_interval, threshold):
    """ Returns true or false if the given  BED interval object
    (pybedtools.Interval) if the value assigned to that interval is greater
    or equal to the threshold given. Attempts to work with both a properly
    formatted bigWig and BED file """
    # NB: a bigwig file will have the "score" in the 4th field while a
    # regular bed file will have it's score in it's normal 5th field

    bed_interval_field_count = len(bed_interval.fields)
    # If there are only 4 fields (bigWig or a bedGraph)
    # Attempt to only use the 4th field for the score to compare against
    if bed_interval_field_count == BIG_WIG_FIELD_COUNT:
        score = get_float_from_filter_score(bed_interval[BIG_WIG_FIELD_COUNT-1])
    # If there are more fields than a bedGraph, assume a bed format
    elif bed_interval_field_count > BIG_WIG_FIELD_COUNT:
        score = get_float_from_filter_score(bed_interval.score)
    # Otherwise if there are too few fields
    else:
        # Raise an error with an appropriate message
        raise ValueError("No value for a threshold found in filter file")

    return score >= threshold


def is_big_wig(filename):
    """ Checks that the given filename refers to a valid bigWig file """
    with open(filename, "rb") as big_wig_file:
        signature_string = big_wig_file.read(BIG_WIG_SIGNATURE_BYTE_SIZE)

    # unpack returns a tuple regardless of length
    # the kent reference checks both little endian and big endian packing
    # of the 4 byte signature
    little_endian_signature = struct.unpack("<L", signature_string)[0]
    big_endian_signature = struct.unpack(">L", signature_string)[0]

    if (little_endian_signature == BIG_WIG_SIGNATURE or
       big_endian_signature == BIG_WIG_SIGNATURE):
        return True

    return False


def get_bed_from_track_file(filename):
    """ Returns a pybedtools.BedTool object from a given filename. Attempts to
    automatically detect the bigWig format and convert"""
    # If the file is a bigWig
    if is_big_wig(filename):
        return bigwig_to_bedgraph(filename)
    # Otherwise read in normal bed file
    else:
        return BedTool(filename)


def main(args=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

"""
load_seq: DESCRIPTION
"""

# Copyright 2008-2014 Michael M. Hoffman <michael.hoffman@utoronto.ca>

import csv
from re import compile, VERBOSE
import sys
import warnings

from argparse import ArgumentParser
from collections import defaultdict

from numpy import frombuffer, uint32
from path import Path
from tabdelim import DictReader
from tables import UInt8Atom

from . import FILE_MODE_CHROMS, FORMAT_VERSION, Genome, __version__
from ._chromosome import SEQ_DTYPE
from ._util import (chromosome_name_map_parser,
                    DEFAULT_CHROMOSOME_NAME_STYLE, FILTERS_GZIP,
                    GENOMEDATA_ENCODING, GenomedataDirtyWarning,
                    ignore_comments, LightIterator, maybe_gzip_open)

MIN_GAP_LEN = 100000
assert not MIN_GAP_LEN % 2  # must be even for division

# sre_constants.MAXREPEAT is 65535, so I have to break repeats into
# two segments
REGEX_SEGMENT_LEN = MIN_GAP_LEN // 2  # max == MAXREPEAT

DNA_LETTERS_UNAMBIG = "ACGTacgt"

SUPERCONTIG_NAME_FMT = "supercontig_%s"

SEQ_ATOM = UInt8Atom()

# https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
AGP_FIELDNAMES = ["object", "object_beg", "object_end", "part_number",
                  "component_type", "col6", "col7", "col8", "col9"]

# See NCBI Assembly reports
ASSEMBLY_REPORT_COMMENT_DELIMITER = "#"
ASSEMBLY_REPORT_FIELD_DELIMITER = "\t"

GAP_COMPONENT_TYPES = frozenset("NU")

GenomicPosition = uint32
GENOMIC_POSITION_0 = GenomicPosition(0)


def read_agp(iterable):
    """
    converts coordinates from 1-based to 0-based
    """
    reader = DictReader(ignore_comments(iterable), AGP_FIELDNAMES)
    for row in reader:
        row["object_beg"] = int(row["object_beg"]) - 1
        row["object_end"] = int(row["object_end"])
        row["part_number"] = int(row["part_number"])

        if row["component_type"] in GAP_COMPONENT_TYPES:
            row["gap_length"] = int(row["col6"])
            row["gap_type"] = row["col7"]
            if row["col8"] == "yes":
                row["linkage"] = True
            else:
                assert row["col8"] == "no"
                row["linkage"] = False
            row["linkage_evidence"] = row["col9"]
        else:
            row["component_id"] = row["col6"]
            row["component_beg"] = int(row["col7"]) - 1
            row["component_end"] = int(row["col8"])
            row["orientation"] = row["col8"]

        yield row


def get_merged_agp_coordinates(agp_iterable):
    """
    merges on an iterator on directly adjacent AGP entries and returns a new
    iterator with start and end coordinates inside a tuple
    assumes AGP entries have been converted to 0-based coordinates and assume
    AGP regions are sorted based on start coordinates
    """
    # Get the first AGP row
    # NB: will raise an exception if the iterator is empty (allow)
    current_row = next(agp_iterable)

    # For all remaining AGP rows
    for agp_row in agp_iterable:
        # If there is an current agp row to potentially merge
        if current_row:
            # If the current row overlaps the region under consideration
            if current_row["object_end"] >= agp_row["object_beg"]:
                # Merge into the current agp row
                # Assumes that the new region's end is larger that our current
                current_row["object_end"] = agp_row["object_end"]
            # Otherwise the regions do not overlap
            else:
                # Return the current region
                yield current_row["object_beg"], current_row["object_end"]
                # Set the new agp row to potentially merge
                current_row = agp_row

    # If there is a remaining row
    if current_row:
        # Return remaining row
        yield current_row["object_beg"], current_row["object_end"]


def create_supercontig(chromosome, index, seq=None, start=None, end=None):
    name = SUPERCONTIG_NAME_FMT % index
    h5file = chromosome.h5file
    where = chromosome.h5group
    supercontig = h5file.create_group(where, name)

    if seq is not None:
        seq_array = frombuffer(seq.encode(GENOMEDATA_ENCODING), SEQ_DTYPE)
        h5file.create_carray(supercontig, "seq", SEQ_ATOM, seq_array.shape)

        # XXXopt: does this result in compression?
        supercontig.seq[...] = seq_array

    if start is None:
        start = chromosome.start
    if end is None:
        end = chromosome.end

    attrs = supercontig._v_attrs
    attrs.start = GenomicPosition(start)
    attrs.end = GenomicPosition(end)


# XXXopt: the all-regex approach is much slower than the hybrid
# approach (3 min to load chr21, 6 min to load chr1), but it is easier
# to understand the code and get it correct
#
# the previous code (r22) might have worked fine. Consider backing down to
# it at some point.
#
# XXXopt: a numpy implementation might be better
re_gap_segment = compile(r"""
(?:([^%s]{%d}[^%s]{%d,})                                  # group(1): ambig
   |                                                      #  OR
   ((?:(?:[%s]+|^)(?:[^%s]{1,%d}[^%s]{,%d}(?![^%s]))*)+)) # group(2): unambig
""" % (DNA_LETTERS_UNAMBIG, REGEX_SEGMENT_LEN,
       DNA_LETTERS_UNAMBIG, REGEX_SEGMENT_LEN,
       DNA_LETTERS_UNAMBIG, DNA_LETTERS_UNAMBIG, REGEX_SEGMENT_LEN,
       DNA_LETTERS_UNAMBIG, REGEX_SEGMENT_LEN-1,
       DNA_LETTERS_UNAMBIG), VERBOSE)


def init_chromosome_start(chromosome):
    """
    doesn't set end, so you can't use create_supercontig(end=None) if you use
    this
    """
    chromosome.attrs.start = GENOMIC_POSITION_0
    chromosome._file_attrs.genomedata_format_version = FORMAT_VERSION


def init_chromosome(chromosome, size):
    init_chromosome_start(chromosome)
    chromosome.attrs.end = GenomicPosition(size)


def read_seq(chromosome, seq):
    supercontig_index = 0

    init_chromosome(chromosome, len(seq))

    for m_segment in re_gap_segment.finditer(seq):
        seq_unambig = m_segment.group(2)
        if seq_unambig:
            span = m_segment.span()
            create_supercontig(chromosome, supercontig_index,
                               seq_unambig, *span)
            supercontig_index += 1
        else:
            assert m_segment.group(1)


def read_assembly(chromosome, infile):
    supercontig_index = 0

    init_chromosome_start(chromosome)

    gap_filtered_agp_rows = (part for part in read_agp(infile)
                             if part["component_type"] not in
                             GAP_COMPONENT_TYPES)

    for start, end in get_merged_agp_coordinates(gap_filtered_agp_rows):
        create_supercontig(chromosome, supercontig_index, None, start, end)
        supercontig_index += 1

    chromosome.attrs.end = GenomicPosition(end)


def size_chromosome(chromosome, size):
    init_chromosome(chromosome, size)
    create_supercontig(chromosome, 0)


def load_sizes(filename):
    res = {}

    with maybe_gzip_open(filename) as infile:
        for line in infile:
            row = line.rstrip().split()
            res[row[0]] = int(row[1])

    return res


def get_num_seq(filenames):
    """
    Run through all files once and count the number of sequences
    """
    res = 0

    for filename in filenames:
        with maybe_gzip_open(filename) as infile:
            for line in infile:
                if line.startswith(">"):
                    res += 1

    return res


def read_chromosome_name_map(assembly_report_file):
    res = {}

    # NB: The format for these files are based on NCBI assembly reports:
    # https://www.ncbi.nlm.nih.gov/assembly/help/#report
    # A copy of the current hg38 assembly report is in the test/data directory

    # NB: These csv files are assumed to be small enough to be held in memory
    assembly_report_lines = assembly_report_file.readlines()

    # Get the header line from the assembly report
    header_index = get_assembly_report_header_index(assembly_report_lines)
    header_line = assembly_report_lines[header_index]

    # Retrive the string after the comment delimiter and split on whitespace
    # for the header fields
    header_fields = header_line.partition(
        ASSEMBLY_REPORT_COMMENT_DELIMITER)[2].strip().split()

    # Return a dictionary with header fields for keys containing a list for
    # their respective column
    for header_field in header_fields:
        res[header_field] = []

    for row in csv.DictReader(assembly_report_lines[header_index+1:],
                              header_fields,
                              delimiter=ASSEMBLY_REPORT_FIELD_DELIMITER):
        for header_field in header_fields:
            res[header_field].append(row[header_field])

    return res


def get_assembly_report_header_index(assembly_report_lines):
    """
    Retrieve the index of the last commented header line from a list of
    assembly report lines which is assumed to be the header for the following
    fields
    If no comment is found, assumes first line is the header line
    """

    for line_index, line in enumerate(assembly_report_lines):
        # If the line does not start with a comment
        if not line.startswith(ASSEMBLY_REPORT_COMMENT_DELIMITER):
            # Stop and use the previous line as header line
            # If the first line was not commented, use that line and assume it
            # has the header information
            return max(0, line_index - 1)

    # If no commented lines were found return the first line index
    return 0


def get_chromosome_name(name, chromosome_name_map, name_style):
    """
    Retrieve alias for a chromsome name based on given map and style

    name: name of the chromsome to translate
    chromsome_name_map: dictionary with keys of styles and ordered lists for
    chromosome names as values
    name_style: the style to translate the given 'name' to
    """

    # If no chromsome name mapping exists
    if not chromosome_name_map:
        # Return given name
        return name
    # Otherwise the chromosome name mapping exists
    else:
        name_position = None

        # Find current name in the map
        for style in chromosome_name_map:
            # If the name was found in the map
            try:
                # Store position in the list
                name_position = chromosome_name_map[style].index(name)
            except ValueError:
                pass

        # If the query was name found
        if name_position is not None:
            # Look up name from style and position in its list
            # Return name lookup
            return chromosome_name_map[name_style][name_position]
        # Otherwise
        else:
            # Raise an exception for an unknown query name
            raise ValueError("Cannot find {} in assembly report".format(name))


def load_seq(gdfilename, filenames, verbose=False, mode=None,
             seqfile_type="fasta", assembly_report_file=None,
             chromosome_name_style=DEFAULT_CHROMOSOME_NAME_STYLE):
    gdpath = Path(gdfilename)

    # load sizes if necessary to figure out number of chromosomes
    if seqfile_type == "sizes":
        assert len(filenames) == 1
        sizes = load_sizes(filenames[0])
    else:
        sizes = None

    # decide whether should use dir or file mode
    if mode is None:
        if seqfile_type == "sizes":
            num_seq = len(sizes)
        elif seqfile_type == "agp":
            num_seq = len(filenames)
        else:
            num_seq = get_num_seq(filenames)

        if num_seq < FILE_MODE_CHROMS:
            mode = "dir"
        else:
            mode = "file"

        if verbose:
            msg = ("Implementation unspecified. Found %d"
                   "chromosomes/scaffolds, so using: %s" % (num_seq, mode))
            print(msg, file=sys.stderr)
    if mode == "dir":
        if gdpath.exists():
            assert gdpath.is_dir()
        else:
            gdpath.makedirs()
    elif mode == "file":
        assert not gdpath.exists()
        gdpath.touch()  # Create file (then Genomedata archive will be a file)
    else:
        raise ValueError("Unsupported mode: %s" % mode)

    # If there's an assembly report file to map chromosome names
    if assembly_report_file:
        # Load the chrosome name mapping
        chromosome_name_map = read_chromosome_name_map(assembly_report_file)
    else:
        chromosome_name_map = None

    # XXX: ignoring all warnings is probably bad. Why is this here?
    # Can we be more specific?
    #
    # is this because we are going to close chromosomes when they are
    # set dirty?
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GenomedataDirtyWarning)
        with Genome(gdpath, mode="w", filters=FILTERS_GZIP) as genome:
            if seqfile_type == "sizes":
                for name, size in sizes.items():
                    name = get_chromosome_name(name, chromosome_name_map,
                                               chromosome_name_style)
                    chromosome = genome._create_chromosome(name, mode)
                    size_chromosome(chromosome, size)
            else:
                assert seqfile_type in frozenset(["agp", "fasta"])
                for filename in filenames:
                    if verbose:
                        print(filename, file=sys.stderr)

                    with maybe_gzip_open(filename) as infile:
                        if seqfile_type == "agp":
                            # Read the entire assembly into a buffer
                            # Filter out comments
                            agp_lines = ignore_comments(infile.readlines())

                            # Split AGP buffer by chromosome entries
                            agp_chromosome_buffer = defaultdict(list)
                            # For every AGP line
                            agp_object_index = AGP_FIELDNAMES.index("object")
                            for agp_line in agp_lines:
                                chr_name = \
                                    agp_line.split("\t")[agp_object_index]

                                # Add the line by chromsome name in the dict
                                agp_chromosome_buffer[chr_name].append(
                                    agp_line)

                            # For each chromosome and its agp lines
                            for chromosome_name in agp_chromosome_buffer:
                                # Create the chromosome in genomedata
                                name = get_chromosome_name(
                                       chromosome_name, chromosome_name_map,
                                       chromosome_name_style)

                                chromosome = genome._create_chromosome(
                                    name, mode)
                                # Read the assembly in to the chromosome entry
                                # in genomedata
                                read_assembly(chromosome,
                                              agp_chromosome_buffer[
                                                  chromosome_name])

                        else:
                            for defline, seq in LightIterator(infile):
                                name = get_chromosome_name(
                                    defline,
                                    chromosome_name_map,
                                    chromosome_name_style)

                                chromosome = genome._create_chromosome(
                                    name, mode)
                                read_seq(chromosome, seq)


def parse_options(args):

    description = ("Start a Genomedata archive at GENOMEDATAFILE with the"
                   " provided sequences. SEQFILEs should be in fasta format,"
                   " and a separate Chromosome will be created for each"
                   " definition line.")

    parser = ArgumentParser(description=description,
                            parents=[chromosome_name_map_parser],
                            prog='genomedata-load-seq')

    parser.add_argument('--version', action='version', version=__version__)

    parser.add_argument('gdarchive', help='genomedata archive')

    parser.add_argument('seqfiles', nargs='+',
                        help='sequences in FASTA format')

    parser.add_argument("-a", "--assembly", action="store_const",
                        const="agp", dest="seqfile_type", default="fasta",
                        help="SEQFILE contains assembly (AGP) files instead of"
                        " sequence")
    parser.add_argument("-s", "--sizes", action="store_const", const="sizes",
                        dest="seqfile_type", default="fasta",
                        help="SEQFILE contains list of sizes instead of"
                        " sequence")
    parser.add_argument("-f", "--file-mode", dest="mode",
                        default=None, action="store_const", const="file",
                        help="If specified, the Genomedata archive will be"
                        " implemented as a single file, with a separate h5"
                        " group for each Chromosome. This is recommended if"
                        " there are a large number of Chromosomes. The default"
                        " behavior is to use a single file if there are at"
                        " least %s Chromosomes being added." %
                        FILE_MODE_CHROMS)
    parser.add_argument("-d", "--directory-mode", dest="mode",
                        action="store_const", const="dir",
                        help="If specified, the Genomedata archive will be"
                        " implemented as a directory, with a separate file for"
                        " each Chromosome. This is recommended if there are a"
                        " small number of Chromosomes. The default behavior is"
                        " to use a directory if there are fewer than %s"
                        " Chromosomes being added." % FILE_MODE_CHROMS)
    parser.add_argument("--verbose", default=False, action="store_true",
                        help="Print status updates and diagnostic messages")

    args = parser.parse_args(args)

    return args


def main(argv=sys.argv[1:]):
    args = parse_options(argv)

    return load_seq(args.gdarchive, args.seqfiles, verbose=args.verbose,
                    mode=args.mode, seqfile_type=args.seqfile_type,
                    assembly_report_file=args.assembly_report,
                    chromosome_name_style=args.name_style)


if __name__ == "__main__":
    sys.exit(main())

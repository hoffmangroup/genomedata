from __future__ import absolute_import
"""
_filter_data_parsers.py: helper generators that produce chromosome regions
based on input files and a given filter function
"""

from string import digits
import re

# NB: None of the file formats officially support comments but there is some
# expectation that it may happen anyway
from ._util import ignore_comments

# See http://genome.ucsc.edu/goldenPath/help/wiggle.html for reference
WIG_VARIABLE_STEP_DEFINITION = "variableStep"
WIG_FIXED_STEP_DEFINITION = "fixedStep"
WIG_DEFAULT_SPAN_VALUE = 1
WIG_VARIABLE_START_INDEX = 0
WIG_VARIABLE_VALUE_INDEX = 1

NUM_BED3_FIELDS = 3
BED_CHROM_NAME_FIELD_INDEX = 0
BED_START_FIELD_INDEX = 1
BED_END_FIELD_INDEX = 2
BED_SCORE_FIELD_INDEX = 4


def passes_filter(filter_function, value):
    """Returns true if the filter doesn't exist or the filter does exist and
    the value evaluates to true on the filter
    """
    return (not filter_function or
            filter_function(value))


def get_wiggle_span(span):
    """If the span exists, return its value otherwise return the default span
    for the wiggle format (default 1)"""
    if span:
        return int(span)
    else:
        return WIG_DEFAULT_SPAN_VALUE


def merged_filter_region_generator(filter_region_generator, filter_file,
                                   filter_function):
    """ Generates merged regions from a given generator that produces
    chromosome regions """

    # NB: These defaults are irrelevant
    current_chromosome = None
    current_start = None
    current_end = None

    # For every filter region
    for chromosome_name, filter_start, filter_end in \
            filter_region_generator(filter_file, filter_function):

        # If we have an initial region
        if current_chromosome:
            # If the next chromsome matches
            # And the start of it overlaps with the current end
            if (chromosome_name == current_chromosome and
               current_end >= filter_start):
                # Merge them into a single region
                current_end = filter_end
            # Otherwise the regions do not overlap
            else:
                # Return the current region
                yield current_chromosome, current_start, current_end
                # Create a new current region with the next region
                current_chromosome = chromosome_name
                current_start = filter_start
                current_end = filter_end
        # Otherwise initialize our initial region
        else:
            current_chromosome = chromosome_name
            current_start = filter_start
            current_end = filter_end

    # Return remaining filter region if there were any
    if current_chromosome:
        yield current_chromosome, current_start, current_end


# All region generators return a tuple (chromosome, start, end)

def get_bed_filter_region(filter_file_handle, filter_function):
    for line in ignore_comments(filter_file_handle):
        valid_line = True
        fields = line.split("\t")

        # If there is a filter function and the line has more than 3 fields
        # e.g. chr1    0       100     A	0.1
        if (len(fields) > NUM_BED3_FIELDS and
           filter_function):
            # Read a score from the BED line
            try:
                score = float(fields[BED_SCORE_FIELD_INDEX])
            # If the score cannot be understood
            except ValueError:
                # Raise an error
                raise ValueError("Could not understand filter score from BED "
                                 "line: {}".format(line))

            valid_line = filter_function(score)

        # If the score passes the filter or there is no filter or score
        if valid_line:
            # Return the result
            yield (fields[BED_CHROM_NAME_FIELD_INDEX],
                   int(fields[BED_START_FIELD_INDEX]),
                   int(fields[BED_END_FIELD_INDEX].rstrip()))


def get_wig_filter_region(filter_file_handle, filter_function):
    # See http://genome.ucsc.edu/goldenPath/help/wiggle.html as reference
    current_span = WIG_DEFAULT_SPAN_VALUE
    # XXX: These initial settings should be irrelevant
    current_wig_definition = WIG_VARIABLE_STEP_DEFINITION
    current_start = 0
    current_step = 1
    current_chromosome = "chr1"

    # variableStep chrom=chrN [span=windowSize]
    variable_definition_regex = re.compile(WIG_VARIABLE_STEP_DEFINITION +
                                           r"\s+chrom=(?P<chromosome>\w+)"
                                           r"(\s+span=(?P<span>\d+))?")

    # fixedStep chrom=chrN start=position step=stepInterval [span=windowSize]
    fixed_definition_regex = re.compile(WIG_FIXED_STEP_DEFINITION +
                                        r"\s+chrom=(?P<chromosome>\w+)"
                                        r"\s+start=(?P<start>\d+)"
                                        r"\s+step=(?P<step>\d+)"
                                        r"(\s+span=(?P<span>\d+))?")

    for line in ignore_comments(filter_file_handle):
        # If the current line starts with a number
        if line[0] in digits:
            # Process a wiggle data line
            # If the current definition is variable
            if current_wig_definition == WIG_VARIABLE_STEP_DEFINITION:
                # Get start coordinate and value
                line_items = line.split()
                current_start = int(line_items[WIG_VARIABLE_START_INDEX])
                value = float(line_items[WIG_VARIABLE_VALUE_INDEX])

                # If a filter exists and the value passes the filter
                if passes_filter(filter_function, value):
                    # NB: Span is the number of elements to include
                    # The end coordinate is exclusive when indexing into
                    # genomedata so it is not necessary to subtract 1
                    # Get the end coordinate based on span
                    end = current_start + current_span

                    # Return chromosome and coordinates
                    yield current_chromosome, int(current_start), int(end)

            # Else (the current definition is a fixed step)
            else:
                # Get the value
                value = float(line.rstrip())
                # If a filter exists and the value passes the filter
                if passes_filter(filter_function, value):
                    # NB: See comment on span above
                    end = current_start + current_span
                    # Return chromosome and coordinates
                    yield current_chromosome, int(current_start), int(end)
                # Update the start coordinate based on fixed step
                current_start += current_step
        # Otherwise the current line is a wiggle definition line
        else:
            # If the current definition is variable step
            re_match = variable_definition_regex.match(line)
            if re_match:
                # Update the current wig definition
                current_wig_definition = WIG_VARIABLE_STEP_DEFINITION
                # Update the current chromsome and span (default 1)
                current_chromosome = re_match.group("chromosome")
                new_span = re_match.group("span")
                # If a span was defined
                # Update the current span
                # Otherwise set the default
                current_span = get_wiggle_span(new_span)

            # If the current definition line is fixed step
            re_match = fixed_definition_regex.match(line)
            if re_match:
                # Update the current wig definition
                current_wig_definition = WIG_FIXED_STEP_DEFINITION
                # Update the current chromsome, start, step and span
                current_chromosome = re_match.group("chromosome")
                current_start = int(re_match.group("start"))
                current_step = int(re_match.group("step"))
                new_span = re_match.group("span")
                # If a span was defined
                # Update the current span
                # Otherwise set the default
                current_span = get_wiggle_span(new_span)

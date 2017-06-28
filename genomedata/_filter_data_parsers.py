import re

# NB: None of the file formats officially support comments but there is some
# expectation that it may happen anyway
from ._util import ignore_comments

WIG_VARIABLE_STEP_DEFINITION = "variableStep"
WIG_FIXED_STEP_DEFINITION = "fixedStep"
WIG_DEFAULT_SPAN_VALUE = 1


def passes_filter(filter_function, value):
    """Returns true if the filter doesn't exist or the filter does exist and
    the value evaluates to true on the filter
    """
    return (not filter_function or (
            filter_function and
            filter_function(value)))


def get_wiggle_span(span):
    """If the span exists, return it's value otherwise return the default span
    for the wiggle format (1)"""
    if span:
        return int(span)
    else:
        return WIG_DEFAULT_SPAN_VALUE


# All region generators return a tuple of chromosome, start, end


def get_bed_filter_region(filter_file_handle, filter_function):
    for line in ignore_comments(filter_file_handle):
        valid_line = True
        fields = line.split("\t")

        # If there is a filter function and the line has more than 3 fields
        # e.g. chr1    0       100     A	0.1
        if (len(fields) > 3 and
           filter_function):
            # Read a score from the BED line
            try:
                # If the score cannot be understood
                score = float(fields[4])
            except:
                # Raise an error
                raise ValueError("Could not understand filter score from BED "
                                 "line: {}".format(line))

            valid_line = filter_function(score)

        # If the score passes the filter or there is no filter or score
        if valid_line:
            # Return the result
            yield fields[0], int(fields[1]), int(fields[2].rstrip())


def get_wig_filter_region(filter_file_handle, filter_function):
    # See http://genome.ucsc.edu/goldenPath/help/wiggle.html as reference
    current_wig_definition = WIG_VARIABLE_STEP_DEFINITION
    current_start = 0
    current_span = 1
    current_step = 1
    current_chromosome = "chr1"  # This default should not matter

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

    digits = map(str, range(10))  # A list of digits from 0 to 9 as str
    for line in ignore_comments(filter_file_handle):
        # If the current line starts with a number
        if line[0] in digits:
            # Process a wiggle data line
            # If the current definition is variable
            if current_wig_definition == WIG_VARIABLE_STEP_DEFINITION:
                # Get start coordinate and value
                line_items = line.split()
                current_start = int(line_items[0])
                value = float(line_items[1])

                # If a filter exists and the value passes the filter
                if passes_filter(filter_function, value):
                    # NB: Span is the number of elements to include
                    # The end coordinate is exclusive when indexing into
                    # genomedata so it is necessary to add by 1
                    # Get the end coordinate based on span
                    end = current_start + current_span + 1

                    # Return chromosome and coordinates
                    yield current_chromosome, int(current_start), int(end)

            # Else (the current definition is a fixed step)
            else:
                # Get the value
                value = float(line.rstrip())
                # If a filter exists and the value passes the filter
                if passes_filter(filter_function, value):
                    # NB: See comment on span above
                    end = current_start + current_span + 1
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
                # Otherwise set the default (1)
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
                # Otherwise set the default (1)
                current_span = get_wiggle_span(new_span)

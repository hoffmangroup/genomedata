#!/usr/bin/env perl

# requires the Inline::Python module from CPAN
# and genomedata installed for Python
use Inline Python => <<'END';

from genomedata import Genome

def open_genomedata(filename):
    # simulate a context manager
    return Genome(filename).__enter__()

def close_genomedata(genome):
    return genome.__exit__(None, None, None)

def get_genomedata(genome, chrom, index, trackname):
    return float(genome[chrom][index, trackname])

END

my $genome = open_genomedata($ARGV[0]);
shift();

while (<STDIN>) {
    chomp();
    my ($chrom, $index) = split();
    foreach my $trackname (@ARGV) {
        print(get_genomedata($genome, $chrom, int($index), $trackname), " ");
    }
    print("\n");
}

close_genomedata($genome);

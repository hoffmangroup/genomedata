#!/usr/bin/env perl

use Inline Python => <<'END';

from genomedata import Genome

def open_genomedata(filename):
    # fake context manager
    return Genome(filename).__enter__()

def close_genomedata(genome)
    return genome.__exit__()

def get_genomedata(genome, chrom, index, trackname):
    return genome[chrom][index, trackname]

END

my $genome = open_genomedata($ARGV[0]);
shift();

while (<>) {
    chomp();
    my ($chrom, $index) = split();
    foreach (my $trackname in @ARGV) {
        print(get_genomedata($genome, $chrom, $index, $trackname), " ");
    }
    print("\n");
}

close_genomedata($genome);

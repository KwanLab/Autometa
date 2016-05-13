#!/usr/bin/perl

# A program to process VizBin temporary files and make a complete table where the point
# coordinates are joined to the contig name
# USAGE: vizbin_process.pl filteredSequences.fa points.txt
# (outputs to stdout)

use strict;
use warnings;

my $filtered_seq_path = $ARGV[0];
my $points_path = $ARGV[1];

# First let's get the order of sequences in the fasta file
my @sequence_names;
open (FASTA, "<$filtered_seq_path") or die "\nCouldn't open $filtered_seq_path\n";
while (<FASTA>)
{
	my $line = $_;
	next unless ( substr( $line, 0, 1 ) eq '>' );

	chomp $line;
	my @seqNameArray = split (' ', substr( $line, 1 ));
	push ( @sequence_names, $seqNameArray[0] );
}
close FASTA;

# Now go through points, and output table
open (POINTS, "<$points_path") or die "\nCouldn't open $points_path\n";
print "contig\tvizbin_x\tvizbin_y\n";
while (<POINTS>)
{
	my $line = $_;
	chomp $line;

	my @coordinateArray = split (',', $line);

	my $seqname = shift( @sequence_names );

	print "$seqname\t$coordinateArray[0]\t$coordinateArray[1]\n";
}
close POINTS;
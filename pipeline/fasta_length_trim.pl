#!/usr/bin/perl

# Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan
#
# This file is part of Autometa.
#
# Autometa is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Autometa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Autometa. If not, see <http://www.gnu.org/licenses/>.

# Program to trim a multifasta according to length of sequence (i.e. cutting out sequences below a certain length)
# USAGE: fasta_length_trim.pl input.fasta <cutoff in bp> output.fasta

use strict;
use warnings;

unless (defined $ARGV[2])
{
	print "\nNot enough arguments\n";
	print "USAGE: fasta_length_trim.pl input.fasta <cutoff in bp> output.fasta\n";
	exit 1;
}

my $input_fasta = $ARGV[0];
my $cutoff = $ARGV[1];
my $output_fasta = $ARGV[2];

open (INPUTFASTA, "<$input_fasta") or die "\nCouldn't open $input_fasta\n";
open (OUTPUTFASTA, ">$output_fasta") or die "\nCouldn't open $output_fasta\n";

my $last_seqname;
my $last_seq;

while (<INPUTFASTA>)
{
	my $line = $_;
	chomp $line;
	
	if ($line =~ m{>}) # The current line is a sequence header
	{
		# First process last sequence
		if (defined $last_seqname)
		{
			my $seq_length = length ($last_seq);
			if ($seq_length >= $cutoff)
			{
				print OUTPUTFASTA "$last_seqname\n$last_seq\n";
			}
		}
	
		$last_seqname = $line;
		$last_seq = ""; #re-initialize sequence
	}
	else
	{
		$last_seq = $last_seq . $line;
	}
}

if (defined $last_seqname)
{
	my $seq_length = length ( $last_seq );
	if ( $seq_length >= $cutoff )
	{
		print OUTPUTFASTA "$last_seqname\n$last_seq\n";
	}
}

close INPUTFASTA;
close OUTPUTFASTA;

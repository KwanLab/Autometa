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

# Program to print a table of sequences in a multifasta file and their gc content
# USAGE: fasta_gc_table.pl <input fasta>
# prints to STDOUT

use strict;
use warnings;

my $input_fasta = $ARGV[0];

my %sequences = loadFasta($input_fasta);

print "Sequence\tlength\n";
foreach my $seq (keys %sequences)
{
	#my $gc_percent = getGC( $sequences{ $seq } );
	my $length = length( $sequences{ $seq } );
	print "$seq\t$length\n";
}

sub getGC
{
	my $seq = $_[0];
	my $length = length ($seq);
	my @seqarray = split (//, $seq);
	my $GC_tot = 0;
	foreach my $nucleotide (@seqarray)
	{
		if (($nucleotide eq "G") || ($nucleotide eq "g") || ($nucleotide eq "C") || ($nucleotide eq "c"))
		{
			$GC_tot++;
		}
	}
	my $GC_percent = ($GC_tot / $length) * 100;
	return $GC_percent;
}

sub loadFasta
{
	my $input_path = $_[0];
	my %fastaHash;
	my $last_seqname;
	my $last_seq;
	open (INPUTFASTA, "<$input_path") or die "\nCouldn't open $input_path\n";
	while (<INPUTFASTA>)
	{
		my $line = $_;
		chomp $line;
		if (substr($line, 0, 1) eq ">")
		{
			if (defined $last_seqname)
			{
				$fastaHash{ $last_seqname } = $last_seq;
			}
			my @linearray = split (" ", $line);
			$last_seqname = substr($linearray[0], 1);
			$last_seq = "";
		} else {
			$last_seq = $last_seq . $line;
		}
	}
	close INPUTFASTA;
	# For last sequence
	$fastaHash{ $last_seqname } = $last_seq;
	return (%fastaHash);
}
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

# Program to calculate and tabulate the average coverage for each contig from a 
# bedtools coverage output (obtained by running genomeCoverageBed -ibam <sorted bam>
# -g genome.txt > coverage.txt, where genome.txt is tab delimited <contig_name>\t<length>)
# The genomeCoverageBed output is a tab delimited file with the following columns:
# 1. Chromosome
# 2. Depth of coverage
# 3. Number of bases on chromosome with that coverage
# 4. Size of chromosome
# 5. Fraction of bases on that chromosome with that coverage
# See also: http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html

# USAGE: contig_coverage_from_bedtools.pl <bedtools output>

use strict;
use warnings;

my $bedtools_output = $ARGV[0];

my %contig_base_coverages;
open (BED, "<$bedtools_output") or die "\nCouldn't open $bedtools_output\n";
while (<BED>)
{
	my $line = $_;
	chomp $line;
	my @lineArray = split ("\t", $line);

	my $current_contig = $lineArray[0];
	my $coverage_total = $lineArray[1] * $lineArray[2];
	my $bases_considered = $lineArray[2];
	# For each line, keep a running total of total coverage
	#  (each time coverage x number of bases)
	#  and total number of bases examined
	if ( defined $contig_base_coverages{ $lineArray[0] } )
	{
		$contig_base_coverages{ $current_contig }->{ 'coverage_total' } += $coverage_total;
		$contig_base_coverages{ $current_contig }->{ 'bases_total' } += $bases_considered;
		} else {
			my %temp_hash = ( 
					'coverage_total' => $coverage_total,
					'bases_total' => $bases_considered );
			$contig_base_coverages{ $current_contig } = \%temp_hash;
		}	
}
close BED;

# Now work out averages and print out coverage table
print "contig\tcoverage\n";
foreach my $contig ( keys %contig_base_coverages )
{
	my %contigHash = %{ $contig_base_coverages{ $contig } };
	my $coverage = $contigHash{ 'coverage_total' } / $contigHash{ 'bases_total' };
	print "$contig\t$coverage\n";
}
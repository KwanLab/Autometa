#!/usr/bin/perl

# Program that takes a MEGAN readname_taxid table (tab delimited) for protein blastp
# results, and maps the taxonomic classifications to contigs (by majority vote, ranks
# considered in descending order from Species).  Outputs a table with the following 
# format:
# contig\tcov\tgc\tlength\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n
# USAGE: megan_protein_to_contig.pl <taxid table> <assembly nucleotide fasta> 
# <taxdump folder> <output file>
# Note: download NCBI taxonomy from:
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

use strict;
use warnings;
use Data::Dumper;
use FindBin;
use POSIX;
use lib "$FindBin::Bin/lib/sequence-toolkit";
use Sequence_toolkit;

my $taxid_table = $ARGV[0];
my $nucleotide_fasta = $ARGV[1];
my $taxdump_dir = $ARGV[2];
my $output_file = $ARGV[3];

# Process taxid tree/names
my $names_dmp = $taxdump_dir . '/names.dmp';
my $nodes_dmp = $taxdump_dir . '/nodes.dmp';

my %taxid;
print "Processing taxid names\n";
open (NAMES, "<$names_dmp") or die "\nCouldn't open $names_dmp\n";
while (<NAMES>)
{
	my $line = $_;
	chomp $line;
	my @lineArray = split (/\|/, $line);
	# Remove trailing and leading spaces
	for (my $i = 0; $i < (@lineArray); $i++)
	{
		$lineArray[$i] =~ s/^\s+//;
		$lineArray[$i] =~ s/\s+$//;
	}

	if ( $lineArray[3] eq 'scientific name' )
	{
		# Replace middle space with underscore
		$lineArray[1] =~ s/ /_/g;
		$taxid{ $lineArray[0] } = {
			'name'	=>	$lineArray[1]
		}
	}
}
close NAMES;

print "Processing taxid nodes\n";
open (NODES, "<$nodes_dmp") or die "\nCouldn't open $nodes_dmp\n";
while (<NODES>)
{
	my $line = $_;
	chomp $line;
	my @lineArray = split (/\|/, $line);
	# Remove trailing and leading spaces
	for (my $i = 0; $i < (@lineArray); $i++)
	{
		$lineArray[$i] =~ s/^\s+//;
		$lineArray[$i] =~ s/\s+$//;
	}

	$taxid{ $lineArray[0] }->{ 'parent' } = $lineArray[1];
	$taxid{ $lineArray[0] }->{ 'rank' } = $lineArray[2];
}
close NODES;

#print Dumper (\%taxid);

print "Parsing megan table\n";
my %protein_classifications;
my %number_of_proteins;
my %canonical_ranks = (
		'superkingdom'	=> 	1,
		'phylum'		=>	1,
		'class'			=>	1,
		'order'			=>	1,
		'family'		=>	1,
		'genus'			=>	1,
		'species'		=>	1,
	);
my @rank_priority = ('species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom');
open (MEGAN, "<$taxid_table") or die "\nCouldn't open $taxid_table\n";
while (<MEGAN>)
{
	my $line = $_;
	chomp $line;
	my @lineTabArray = split ("\t", $line);
	my @seqUnderscoreArray = split ('_', $lineTabArray[0]);
	pop(@seqUnderscoreArray);
	#pop(@seqUnderscoreArray); # <--- Diff between megan_protein_to_contig.pl and megan_protein_to_contig_3.pl
	my $contig = join ('_', @seqUnderscoreArray);
	my $line_taxid = $lineTabArray[1];

	# Protect for edge case - when the taxid retrieved is not yet in taxdump - treat as if taxid 1
	unless (defined $taxid{ $line_taxid })
	{
		$line_taxid = 1;
	}

#### - Below is new code
	
	# New voting scheme - for each taxid in the megan table (for one ORF),
	# each taxid in the whole taxon path is counted (as long as it is a 
	# canonical rank).  At the end the contig classification will be determined
	# as the majority vote of the highest rank in which a majority exists.

	unless ( defined $protein_classifications{ $contig } )
	{
		# Initialize variables
		$protein_classifications{ $contig } = {};
		foreach my $rank ( @rank_priority )
		{
			$protein_classifications{ $contig }->{ $rank } = {};
		}
		$number_of_proteins{ $contig } = 0;
	}

	$number_of_proteins{ $contig }++;

	my $current_taxid = $line_taxid;
	
	until ( $current_taxid == 1 )
	{
		if ( exists $canonical_ranks{ $taxid{ $current_taxid }->{ 'rank' } } )
		{
			my $rank = $taxid{ $current_taxid }->{ 'rank' };

			if ( defined $protein_classifications{ $contig }->{ $rank }->{ $current_taxid } )
			{
				$protein_classifications{ $contig }->{ $rank }->{ $current_taxid }++;
			} else {
				$protein_classifications{ $contig }->{ $rank }->{ $current_taxid } = 1;
			}
		}

		$current_taxid = $taxid{ $current_taxid }->{ 'parent' };
	}
}
close MEGAN;

#print Dumper (\%number_of_proteins);
#print Dumper (\%protein_classifications);
#exit;

# Determine top ranked taxid
print "Ranking taxids\n";
my %top_taxids;
foreach my $contig ( keys %protein_classifications )
{
	my $best_taxid;
	foreach my $rank ( @rank_priority )
	{
		if ( defined $protein_classifications{ $contig }->{ $rank } )
		{
			# Work out if there is a majority
			my $total_proteins = $number_of_proteins{ $contig };
			my $majority = ceil($total_proteins / 2);

			my $majority_vote;

			foreach my $taxid ( keys %{ $protein_classifications{ $contig }->{ $rank } } )
			{
				my $votes = $protein_classifications{ $contig }->{ $rank }->{ $taxid };
				if ( $votes >= $majority )
				{
					$majority_vote = $taxid;
				}
			}
			
			if ( defined $majority_vote )
			{
				$best_taxid = $majority_vote;
				last;
			}
		}
	}

	if ( defined $best_taxid )
	{
		$top_taxids{ $contig } = $best_taxid;
	} else {
		$top_taxids{ $contig } = 1;
	}
}

#print Dumper ( \%top_taxids );
#exit;

# Resolve taxon paths
print "Resolving taxon paths\n";
my %taxon_paths;
foreach my $contig ( keys %top_taxids )
{
	my %taxon_path;
	my $current_taxid = $top_taxids{ $contig };

	#May need to put a check at this point to see if taxid is defined in taxdump

	until ($current_taxid == 1)
	{
		my $rank = $taxid{ $current_taxid }->{ 'rank' };
		$taxon_path{ $rank } = $current_taxid;
		$current_taxid = $taxid{ $current_taxid }->{ 'parent' };
	}

	#Resolve taxon path to names
	my %taxon_path_names;
	foreach my $rank ( keys %taxon_path )
	{
		my $tax_name = $taxid{ $taxon_path{ $rank } }->{ 'name' };
		$taxon_path_names{ $rank } = $tax_name;
	}

	$taxon_paths{ $contig } = \%taxon_path_names;
}

# Load fasta file
print "Loading fasta file\n";
my $fastaHashRef = Sequence_toolkit::fastaToHash($nucleotide_fasta);

# Determine ranks and make output file

print "Writing table\n";
open (OUTPUT, ">$output_file");
print OUTPUT "contig\tcov\tgc\tlength\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\ttaxid\n";
foreach my $contig (keys %{ $fastaHashRef })
{
	my $seq = $fastaHashRef->{ $contig };
	my $gc = Sequence_toolkit::gc_percent( $seq );
	my $length = length( $seq );
	my @contigNameArray = split ('_', $contig);
	my $cov = $contigNameArray[5];

	my %corrected_taxon_path;
	if ( defined $taxon_paths{ $contig } )
	{
		foreach my $rank ( @rank_priority )
		{
			if ( defined $taxon_paths{ $contig }->{ $rank } )
			{
				$corrected_taxon_path{ $rank } = $taxon_paths{ $contig }->{ $rank };
			} else {
				$corrected_taxon_path{ $rank } = 'unclassified';
			}
		}
	} else {
		foreach my $rank ( @rank_priority )
		{
			$corrected_taxon_path{ $rank } = 'unclassified';
		}
	}

	my $output_string = $contig . "\t" . $cov . "\t" . $gc . "\t" . $length . "\t" .
							$corrected_taxon_path{ 'superkingdom' } . "\t" .
							$corrected_taxon_path{ 'phylum' } . "\t" .
							$corrected_taxon_path{ 'class' } . "\t" .
							$corrected_taxon_path{ 'order' } . "\t" .
							$corrected_taxon_path{ 'family' } . "\t" .
							$corrected_taxon_path{ 'genus' } . "\t" .
							$corrected_taxon_path{ 'species' } . "\t";

	if ( defined $top_taxids{ $contig } )
	{
		$output_string = $output_string . $top_taxids{ $contig };
	} else {
		$output_string = $output_string . 'unclassified';
	}

	$output_string = $output_string . "\n";

	print OUTPUT $output_string;
}
close OUTPUT;
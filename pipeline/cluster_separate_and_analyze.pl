#!/usr/bin/perl

# Program that takes a table with clustered contigs and outputs separate fastas
# After making the separate fastas the program will then work out genome completeness
# and create a table of results.

# Dependencies - hmmscan, prodigal (in path)

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $input_fasta;
my $input_table;
my $output_directory;
my $marker_hmm_db;
my $cutoff_file;
my $cluster_column = 'db.cluster';

GetOptions (
	"fasta|f=s"		=> \$input_fasta,
	"table|t=s"		=> \$input_table,
	"outputdir|o=s"	=> \$output_directory,
	"hmmdb|h=s"		=> \$marker_hmm_db,
	"cutoffs|c=s"	=> \$cutoff_file,
	"column|n=s"	=> \$cluster_column,
	);

# First, define clusters
my %cluster_fastas; # Will hold the fastas for each cluster
my %contig_classifications; # Will hold the cluster of each contig in the data set

# Note - by default we expect here a column called db.cluster and a column called contig
open (TABLE, "<$input_table") or die "\nCouldn't open $input_table\n";
my $first_line = <TABLE>;
chomp $first_line;
my @first_line_array = split ("\t", $first_line);
my $cluster_index = 0;
my $contig_index = 0;

for ( my $i = 0; $i < scalar( @first_line_array ); $i++ )
{
	$cluster_index = $i if ( $first_line_array[$i] eq $cluster_column );
	$contig_index = $i if ( $first_line_array[$i] eq 'contig' );
}

while (<TABLE>)
{
	my $line = $_;
	chomp $line;
	my @lineArray = split ("\t", $line);
	my $contig = $lineArray[$contig_index];
	my $cluster = $lineArray[$cluster_index];

	$contig_classifications{ $contig } = $cluster;
	unless ( defined $cluster_fastas{ $cluster } )
	{
		# Make a filename for the cluster fasta
		my $filepath = $output_directory . '/cluster_' . $cluster . '.fasta';
		$cluster_fastas{ $cluster } = $filepath;
	} 
}
close TABLE;

# Now we need to separate out the fasta file

# First open filehandles
my %filehandles;
my %cluster_sequences;
foreach my $cluster ( keys %cluster_fastas )
{
	open ($filehandles{ $cluster }, ">$cluster_fastas{$cluster}") or die "\nCouldn't open $cluster_fastas{$cluster}\n";
}

open (FASTA, "<$input_fasta") or die "\nCouldn't open $input_fasta\n";
my $current_cluster;
my $current_contig;
my $current_sequence;
while (<FASTA>)
{
	my $line = $_;
	if ( substr( $line, 0, 1 ) eq '>' )
	{
		# Work out current cluster
		chomp $line;
		my @lineArray = split (' ', substr( $line, 1 ));
		$current_contig = $lineArray[0];
		unless ( defined $contig_classifications{ $current_contig } )
		{
			# This is an unclaimed contig
			$current_cluster = undef;
			next;
		}

		$current_cluster = $contig_classifications{ $current_contig };

		print {$filehandles{ $current_cluster }} "$line\n";

		$current_sequence = '';
		if ( defined $cluster_sequences{ $current_cluster } )
		{
			$cluster_sequences{ $current_cluster }->{ $current_contig } = '';
		} else {
			$cluster_sequences{ $current_cluster } = {
				$current_contig => '',
			};
		}
	} else {
		next unless ( defined $current_cluster );

		print {$filehandles{ $current_cluster }} $line;

		my $sequence_so_far = $cluster_sequences{ $current_cluster }->{ $current_contig };
		chomp $line;
		$cluster_sequences{ $current_cluster }->{ $current_contig } = $sequence_so_far . $line;
	}
}
close FASTA;

# Close filehandles
foreach my $cluster ( keys %filehandles )
{
	close $filehandles{ $cluster };
}

# For each fasta we need to run prodigal
my %cluster_protein_fastas;
foreach my $cluster ( keys %cluster_fastas )
{
	my $protein_fasta_path = $output_directory . '/cluster_' . $cluster . '_orf.faa';
	my $cluster_fasta = $cluster_fastas{ $cluster };
	my $output_file = $output_directory . '/cluster_' . $cluster . '.gbk';
	my $nucleotide_fasta_path = $output_directory . '/cluster_' . $cluster . '_orf.fasta';
	my $summary_file = $output_directory . '/cluster_' . $cluster . '_prodigal_summary';

	my $prodigal_command = "prodigal -a $protein_fasta_path -i $cluster_fasta -o $output_file -d $nucleotide_fasta_path -s $summary_file -f gbk -p meta";
	my $command_failed = system( $prodigal_command );

	if ( $command_failed )
	{
		print "Error! Prodigal run on cluster $cluster failed\n";
		exit 1;
	}

	$cluster_protein_fastas{ $cluster } = $protein_fasta_path;
}



# Now carry out the hmmscans
my %cluster_hmmscan_tables;
foreach my $cluster ( keys %cluster_protein_fastas )
{
	my $input_protein_fasta = $cluster_protein_fastas{ $cluster };
	my $hmm_output_path = $output_directory . '/cluster_' . $cluster . '_proteins_hmmscan.tbl';

	my $hmmscan_command = "hmmscan --tblout $hmm_output_path $marker_hmm_db $input_protein_fasta";
	my $command_failed = system( $hmmscan_command );

	if ( $command_failed )
	{
		print "Error! Hmmscan run on cluster $cluster failed\n";
		exit 1;
	}

	$cluster_hmmscan_tables{ $cluster } = $hmm_output_path;
}

# Now analyze the hmm output
# Parse cutoffs 
my %cutoffs; # Hash, keyed by accession, stores cutoff scores
open (CUTOFF, "<$cutoff_file") or die "\nCouldn't open $cutoff_file\n";
while (<CUTOFF>)
{
	my $line = $_;
	chomp $line;
	my @linearray = split (' ', $line);
	$cutoffs{ $linearray[0] } = $linearray[1];
}
close CUTOFF;

# Now parse hmm tables
my %cluster_completeness_stats;
foreach my $cluster ( keys %cluster_hmmscan_tables )
{
	my $hmm_table = $cluster_hmmscan_tables{ $cluster };
	open (HMM, "<$hmm_table") or die "\nCouldn't open $hmm_table\n";
	my %hits;

	while (<HMM>)
	{
		my $line = $_;
		next if ( substr( $line, 0, 1 ) eq '#' ); # Header line
		chomp $line;
		my @linearray = split (/ +/, $line);

		my @accessionArray = split(/\./, $linearray[1]);
		my $accession = $accessionArray[0];

		my $score = $linearray[5];
		my $name = $linearray[2];

		if (( defined ($cutoffs{ $accession }) ) && ( $score > $cutoffs{ $accession } ))
		{
			if ( defined $hits{ $accession } )
			{
				push ( @{ $hits{ $accession } }, $name );
				} else {
					$hits{ $accession } = [ $name ];
				}
		}		
	}
	close HMM;

	# Calculate some stats
	my $total_accessions = scalar( keys %cutoffs );
	my $number_found = scalar( keys %hits );
	my $genome_completeness = ( $number_found / $total_accessions ) * 100;
	my $duplicated_accessions = 0;
	my $missing_accessions = 0;

	foreach my $accession ( keys %cutoffs )
	{
		if ( defined ( $hits{ $accession } ) )
		{
			my $number_of_hits = scalar( @{ $hits{ $accession } } );
			$duplicated_accessions++ if ( $number_of_hits > 1 );
		} else {
			$missing_accessions++;
		}
	}

	# Calculate some sequence stats
	my $sequence_stat_hashRef = assess_asm( $cluster_sequences{ $cluster } );

	my %tmp_hash = (
		'completeness'		=> $genome_completeness,
		'number_markers'	=> $number_found,
		'duplicated'		=> $duplicated_accessions,
		'missing'			=> $missing_accessions,
		'size'				=> $sequence_stat_hashRef->{ 'size' },
		'n50'				=> $sequence_stat_hashRef->{ 'n50' },
		'number_contigs'	=> $sequence_stat_hashRef->{ 'number_contigs' },
		'largest_contig'	=> $sequence_stat_hashRef->{ 'largest_contig' },
		'gc_mean'			=> $sequence_stat_hashRef->{ 'gc_mean' },
		'gc_sd'				=> $sequence_stat_hashRef->{ 'gc_sd' },
		'cov_mean'			=> $sequence_stat_hashRef->{ 'cov_mean' },
		'cov_sd'			=> $sequence_stat_hashRef->{ 'cov_sd' },
		);

	$cluster_completeness_stats{ $cluster } = \%tmp_hash;
}

# Now write out the table
my $output_table_path = $output_directory . '/cluster_summary_table';
open (OUTPUT, ">$output_table_path") or die "\nCouldn't open $output_table_path\n";
print OUTPUT "cluster\tsize\tn50\tnumber_contigs\tlargest_contig\tgc_mean\tgc_sd\tcov_mean\tcov_sd\tcompleteness\tnumber_markers\tduplicated\tmissing\n";

foreach my $cluster ( keys %cluster_completeness_stats )
{
	my $completeness = $cluster_completeness_stats{ $cluster }->{ 'completeness' };
	my $size = $cluster_completeness_stats{ $cluster }->{ 'size' };
	my $n50 = $cluster_completeness_stats{ $cluster }->{ 'n50' };
	my $number_contigs = $cluster_completeness_stats{ $cluster }->{ 'number_contigs' };
	my $largest_contig = $cluster_completeness_stats{ $cluster }->{ 'largest_contig' };
	my $gc_mean = $cluster_completeness_stats{ $cluster }->{ 'gc_mean' };
	my $gc_sd = $cluster_completeness_stats{ $cluster }->{ 'gc_sd' };
	my $cov_mean = $cluster_completeness_stats{ $cluster }->{ 'cov_mean' };
	my $cov_sd = $cluster_completeness_stats{ $cluster }->{ 'cov_sd' };
	my $number = $cluster_completeness_stats{ $cluster }->{ 'number_markers' };
	my $duplicated = $cluster_completeness_stats{ $cluster }->{ 'duplicated' };
	my $missing = $cluster_completeness_stats{ $cluster }->{ 'missing' };

	print OUTPUT "$cluster\t$size\t$n50\t$number_contigs\t$largest_contig\t$gc_mean\t$gc_sd\t$cov_mean\t$cov_sd\t$completeness\t$number\t$duplicated\t$missing\n";
}
close OUTPUT;

sub assess_asm
{
	my $inputHashRef = $_[0]; # $clustersequences{ $cluster }

	my @sequences_in_descending_length;
	my @contignames;

	my $total_length = 0;
	for my $contig ( sort { length( $inputHashRef->{ $b } ) <=> length( $inputHashRef->{ $a } ) } keys %{ $inputHashRef } )
	{
		push ( @sequences_in_descending_length, $inputHashRef->{ $contig } );
		push ( @contignames, $contig );
		$total_length += length( $inputHashRef->{ $contig } );
	}

	my $number_of_contigs = scalar( @sequences_in_descending_length );
	my $largest_contig = length( $sequences_in_descending_length[0] );

	# Calculate N50
	my $length_counted = 0;
	my $n50;
	foreach my $seq ( @sequences_in_descending_length )
	{
		$length_counted += length( $seq );
		if ( $length_counted >= ( $total_length / 2 ) )
		{
			$n50 = length( $seq );
			last;
		}
	}

	my @gc_percents;
	my @coverages;
	foreach my $contig ( keys( %{ $inputHashRef} ) )
	{
		# Calculate gc percent
		my $seq = uc( $inputHashRef->{ $contig } );
		my $seq_length = length( $seq );
		my @seqArray = split (//, $seq);
		my $gc_count = 0;

		foreach my $char ( @seqArray )
		{
			$gc_count++ if ( $char eq 'G' || $char eq 'C' || $char eq 'S' );
			$gc_count += 0.5 if ( $char eq 'N' || $char eq 'R' || $char eq 'Y' || $char eq 'K' || $char eq 'M' );
			$gc_count += 0.66 if ( $char eq 'B' || $char eq 'V' );
			$gc_count += 0.33 if ( $char eq 'D' || $char eq 'H' );
		}

		my $gc_percent = ( $gc_count / $seq_length ) * 100;
		push ( @gc_percents, $gc_percent );

		# Get coverage from spades contig name
		my @contigNameArray = split ('_', $contig);
		my $coverage = $contigNameArray[5];
		push ( @coverages, $coverage );
	}

	my ( $gc_mean, $gc_sd ) = mean_sd( \@gc_percents );
	my ( $cov_mean, $cov_sd ) = mean_sd( \@coverages );

	my %tmpHash = (
		'size'				=>	$total_length,
		'n50'				=>	$n50,
		'number_contigs'	=>	$number_of_contigs,
		'largest_contig'	=>	$largest_contig,
		'gc_mean'			=>	$gc_mean,
		'gc_sd'				=>	$gc_sd,
		'cov_mean'			=>	$cov_mean,
		'cov_sd'			=>	$cov_sd,
		);

	return \%tmpHash;
}		

sub mean_sd
{
	my $valueArrayRef = $_[0];

	my $total = 0;
	my $n = scalar( @{ $valueArrayRef } );

	foreach my $value ( @{ $valueArrayRef } )
	{
		$total += $value;
	}

	my $mean = $total / $n;

	if ( $n < 2 )
	{
		return ( $mean, 'NA' );
	}

	my $sum_of_squares = 0;

	foreach my $value ( @{ $valueArrayRef } )
	{
		$sum_of_squares += ( $value - $mean ) ** 2;
	}

	my $standard_deviation = sqrt( $sum_of_squares / ( $n - 1 ) );

	return ( $mean, $standard_deviation );
}

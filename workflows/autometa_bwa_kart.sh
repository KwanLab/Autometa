#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH --error=autometa.%J.err
#SBATCH --output=autometa.%J.out

Help()
{
   # Display Help
   echo "Required flags"
   echo
   echo "Syntax: scriptTemplate [-d|a|s|f|r|n|m|l|c]"
   echo "-h     Display this help menu"
   echo "Required flags:"
   echo "-d     Path to where your input files are and where the output will go"
   echo "-a     Full path to where your scaffold/contig file is (inlclude file name)"
   echo "-s     Simple name to prefix output files"
   echo "-f     DNA forward reads"
   echo "-r     DNA reverse reads"
   echo "-n     Path to NCBI databases"
   echo "-m     Path to marker files"
   echo "-l     Contig length cutoff in bp. Change this according to the N50 you got from running metaQuast"
   echo "-c     How many cpus you would like to use"
   echo
   echo "Example usage:"
   echo "sbatch autometa_bwa_kart.sh -d /home/user/Trial/Sample1 -a /home/user/Trial/Sample1/scaffolds.fasta -s Sample1 -f /home/user/Trial/Sample1/sample1_fwd.fastq -r /home/user/Trial/Sample1/sample1_rev.fastq -n /home/user/Databases -m /home/user/Autometa/autometa/databases/markers -l 1000 -c 16"
}

#Define flags
while getopts "d:a:s:f:r:n:m:l:c" flag
do
    case "${flag}" in
        d) work_dir=${OPTARG};;
        a) assembly=${OPTARG};;
        s) simple_name=${OPTARG};;
        f) fwd_reads=${OPTARG};;
        r) rev_reads=${OPTARG};;
        n) ncbi=${OPTARG};;
        m) marker_db=${OPTARG};;
        l) length_cutoff=${OPTARG};;
        c) cpus=${OPTARG};;
                h) # display Help
                        Help
                        exit;;
    esac
done

# Step 0: Do some Path handling with the provided `assembly` filepath
outdir="${work_dir}/${simple_name}_Autometa_Output"
if [ ! -d $outdir ]
then mkdir -p $outdir
fi

# Autometa Parameters - If you don't know what you're doing, leave these alone!!
kmer_size=5
norm_method="am_clr" # choices: am_clr, clr, ilr
pca_dimensions=50 # NOTE: must be greater than $embed_dimensions
embed_dimensions=2 # NOTE: must be less than $pca_dimensions
embed_method="bhsne" # choices: bhsne, sksne, umap, densne, trimap
# Binning clustering method
cluster_method="hdbscan" # choices: hdbscan, dbscan
# Binning metrics cutoffs
completeness=20.0
purity=95.0
cov_stddev_limit=25.0
gc_stddev_limit=5.0
seed=42

######### BEGIN #########

# Step 00: Report autometa version
autometa --version

# Step 1: filter assembly by length and retrieve contig lengths as well as GC content

# output:
filtered_assembly="${outdir}/${simple_name}.filtered.fna"
gc_content="${outdir}/${simple_name}.gc_content.tsv"

# script:
autometa-length-filter \
    --assembly $assembly \
    --cutoff $length_cutoff \
    --output-fasta $filtered_assembly \
    --output-gc-content $gc_content

# Step 2: Determine coverages from assembly read alignments

alignment_sam_file="${outdir}/${simple_name}.alignment.sam"
alignment_bam_file="${outdir}/${simple_name}.alignment.bam"

# First index metagenome assembly
bwa index \\
    -b 550000000000 \\ # block size for the bwtsw algorithm (effective with -a bwtsw) [default=10000000]
    $filtered_assembly     # Path to input metagenome

# Now perform alignments (we are using kart, but you can use another alignment tool if you'd like)
kart \\
    -i $filtered_assembly                   \\ # Path to input metagenome
    -t $cpus                               \\ # Number of cpus to use
    -f $fwd_reads  \\ # Path to forward paired-end reads
    -f2 $rev_reads \\ # Path to reverse paired-end reads
    -o $alignment_sam_file   # Path to alignments output

# Now sort alignments and convert to bam format
samtools sort \\
    -@ $cpus              \\ # Number of cpus to use
    -m 10G             \\ # Amount of memory to use
    $alignment_sam_file     \\ # Input alignments file path
    -o $alignment_bam_file     # Output alignments file path

# output:
coverages="${outdir}/${simple_name}.coverages.tsv"

# script:
autometa-coverage \
    --assembly $filtered_assembly \
    --bam $alignment_bam_file \
    --out $coverages \
    --cpus $cpus

# Step 3: Annotate and filter markers
# input:
# $orfs --> User input
# $cpus --> User input
# $seed --> User input

orfs="${outdir}/${simple_name}.orfs.faa"
orfs_nuc="${outdir}/${simple_name}.orfs.fna"

#Get ORFs
autometa-orfs \
    --assembly $filtered_assembly \
    --output-nucls $orfs_nuc \
    --output-prots $orfs \
    --cpus $cpus

autometa-config \
	--section databases \
	--option markers \
	--value $marker_db

#Check if marker files are pressed, if not, press them
if [ -f "$marker_db/bacteria.single_copy.hmm.h3i" ]; then
    echo "Hmm files are pressed, moving on..."
else
    hmmpress -f $marker_db/bacteria.single_copy.hmm
fi


if [ -f "$marker_db/archaea.single_copy.hmm.h3i" ]; then
    echo "Hmm files are pressed, moving on..."
else
    hmmpress -f $marker_db/archaea.single_copy.hmm
fi

kingdoms=(bacteria archaea)

# NOTE: We iterate through both sets of markers for binning both bacterial and archeal kingdoms
for kingdom in ${kingdoms[@]};do
    # kingdom-specific output:
    hmmscan="${outdir}/${simple_name}.${kingdom}.hmmscan.tsv"
    markers="${outdir}/${simple_name}.${kingdom}.markers.tsv"

    # script:
    autometa-markers \
        --orfs $orfs \
        --hmmscan $hmmscan \
        --out $markers \
        --kingdom $kingdom \
        --parallel \
        --cpus $cpus \
        --seed $seed
done

autometa-config \
    --section databases --option ncbi \
    --value $ncbi

#Check if ncbi files are present, if not, download and prep them
if [ -f "$ncbi/nr.dmnd" ]; then
    echo "NCBI files are present where specified, moving on..."
else
    autometa-update-databases --update-ncbi
fi

# Step 4.1: Determine ORF lowest common ancestor (LCA) amongst top hits

# input:
# $blast --> User Input
# $ncbi --> User Input

#Run blastp
blast="${outdir}/${simple_name}.blastp.tsv" #Generate output file name

diamond blastp \
    --query $orfs \
    --db "${ncbi}/nr.dmnd" \
    --evalue 1e-5 \
    --max-target-seqs 200 \
    --threads $cpus \
    --outfmt 6 \
    --out $blast

# output:
lca="${outdir}/${simple_name}.orfs.lca.tsv"
sseqid2taxid="${outdir}/${simple_name}.orfs.sseqid2taxid.tsv"
errorTaxids="${outdir}/${simple_name}.orfs.errortaxids.tsv"

# script:
autometa-taxonomy-lca \
    --blast $blast \
    --dbdir $ncbi \
    --lca-output $lca \
    --sseqid2taxid-output $sseqid2taxid \
    --lca-error-taxids $errorTaxids

# Step 4.2: Perform Modified Majority vote of ORF LCAs for all contigs that returned hits in blast search

# input:
# $lca --> Generated by step 4.1
# $ncbi --> User Input

# output:
votes="${outdir}/${simple_name}.taxids.tsv"

# script:
autometa-taxonomy-majority-vote --lca $lca --output $votes --dbdir $ncbi

# Step 4.3: Split assigned taxonomies into kingdoms

# input:
# $votes --> Generated by step 4.2
# $outdir --> Generated by step 0
# $ncbi --> User Input
# $assembly --> User Input

# output:
# Will write recovered superkingdoms to ${outdir}
# e.g. ${outdir}/${simple_name}.bacteria.fna
# e.g. ${outdir}/${simple_name}.archaea.fna
# e.g. ${outdir}/${simple_name}.taxonomy.tsv

# script:
autometa-taxonomy \
    --votes $votes \
    --output $outdir \
    --prefix $simple_name \
    --split-rank-and-write superkingdom \
    --assembly $filtered_assembly \
    --ncbi $ncbi

# Step 5: Perform k-mer counting on respective kingdoms

# input:
# $kmer_size --> User input
# $norm_method --> User input
# $embed_method --> User input
# $embed_dimensions --> User input
# $cpus --> User input
# $seed --> User input

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do
    # kingdom-specific input:
    # NOTE: $fasta --> Generated by step 4.3
    fasta="${outdir}/${simple_name}.${kingdom}.fna"

    # kingdom-specific output:
    counts="${outdir}/${simple_name}.${kingdom}.${kmer_size}mers.tsv"
    normalized="${outdir}/${simple_name}.${kingdom}.${kmer_size}mers.${norm_method}.tsv"
    embedded="${outdir}/${simple_name}.${kingdom}.${kmer_size}mers.${norm_method}.${embed_method}.tsv"

    if [ ! -f $fasta ]
    then
        echo "${fasta} does not exist, skipping..."
        continue
    fi
    # script:
    autometa-kmers \
        --fasta $fasta \
        --kmers $counts \
        --size $kmer_size \
        --norm-output $normalized \
        --norm-method $norm_method \
        --pca-dimensions $pca_dimensions \
        --embedding-output $embedded \
        --embedding-method $embed_method \
        --embedding-dimensions $embed_dimensions \
        --cpus $cpus \
        --seed $seed
done

# Step 6: Perform binning on each set of bacterial and archaeal contigs

# input:
# $cpus --> User input
# $seed --> User input
taxonomy="${outdir}/${simple_name}.taxonomy.tsv" # Generated by step 4.3
# $gc_content --> Generated by step 1

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do
    # kingdom-specific input:
    kmers="${outdir}/${simple_name}.${kingdom}.${kmer_size}mers.${norm_method}.${embed_method}.tsv" # Generated by step 5
    markers="${outdir}/${simple_name}.${kingdom}.markers.tsv" # Generated by step 3

    # kingdom-specific output:
    output_binning="${outdir}/${simple_name}.${kingdom}.${cluster_method}.tsv"
    output_main="${outdir}/${simple_name}.${kingdom}.${cluster_method}.main.tsv"

    if [ -f $output_main ] && [ -s $output_main ];then
        echo "$(basename $output_main) already exists. continuing..."
        continue
    fi
    # script:
    autometa-binning \
        --kmers $kmers \
        --coverages $coverages \
        --gc-content $gc_content \
        --markers $markers \
        --output-binning $output_binning \
        --output-main $output_main \
        --clustering-method $cluster_method \
        --completeness $completeness \
        --purity $purity \
        --cov-stddev-limit $cov_stddev_limit \
        --gc-stddev-limit $gc_stddev_limit \
        --taxonomy $taxonomy \
        --starting-rank superkingdom \
        --rank-filter superkingdom \
        --rank-name-filter $kingdom
done

# Step 7: Create binning summary files

# input:
# $ncbi -> User input
# $assembly -> User input

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do

    # kingdom-specific input:
    binning_main="${outdir}/${simple_name}.${kingdom}.${cluster_method}.main.tsv" # Generated by step 6
    markers="${outdir}/${simple_name}.${kingdom}.markers.tsv" # Generated by step 3

    # kingdom-specific output:
    output_stats="${outdir}/${simple_name}_${kingdom}_metabin_stats.tsv"
    output_taxonomy="${outdir}/${simple_name}_${kingdom}_metabin_taxonomy.tsv"
    output_metabins="${outdir}/${simple_name}_${kingdom}_metabins"

    if [ ! -f $binning_main ]
    then
        echo "${binning_main} does not exist, skipping..."
        continue
    fi

    # script:
    autometa-binning-summary \
        --binning-main $binning_main \
        --markers $markers \
        --metagenome $filtered_assembly \
        --ncbi $ncbi \
        --output-stats $output_stats \
        --output-taxonomy $output_taxonomy \
        --output-metabins $output_metabins
done

#########  END  #########

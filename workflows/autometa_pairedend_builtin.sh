#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH --error=autometa.%J.err
#SBATCH --output=autometa.%J.out
#SBATCH --mail-user=samche42@gmail.com
#SBATCH --mail-type=ALL

#Change the following to suit YOUR files:

work_dir="/home/sam/Autometa_test" #work_dir is where your input files are and where the output will go

assembly="/home/sam/Autometa_test/N30_scaffolds.fasta" #assembly is the metagenome fasta file with a full path to where this file is located

fwd_reads="/home/sam/Autometa_test/N30_1P.fastq" #Forward reads used for assembly

rev_reads="/home/sam/Autometa_test/N30_2P.fastq" #Reverse reads used for assembly

ncbi="/media/bigdrive1/Databases/autometa_databases" #ncbi is the full path to where the NCBI databases are found. This does not include the file name.

marker_db="/home/sam/Tools/Autometa/autometa/databases/markers" #marker spcifies were the marker files are. They should have been downloaded with autometa.

simpleName="N30" #A prefix for all output files specific to your sample

length_cutoff=1000 # in bp. Change this according to the N50 you got from running metaQuast.

cpus=16 #How many cpus you would like to use

#IMPORTANT: YOU DO NOT NEED TO CHANGE ANYTHING ELSE FROM THIS POINT ONWARD!

# Step 0: Do some Path handling with the provided `assembly` filepath
outdir="${work_dir}/${simpleName}_Autometa_Output"
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
filtered_assembly="${outdir}/${simpleName}.filtered.fna"
gc_content="${outdir}/${simpleName}.gc_content.tsv"

# script:
autometa-length-filter \
    --assembly $assembly \
    --cutoff $length_cutoff \
    --output-fasta $filtered_assembly \
    --output-gc-content $gc_content

# Step 2: Determine coverages from assembly read alignments

# output:
coverages="${outdir}/${simpleName}.coverages.tsv"

# script:
autometa-coverage \
    --assembly $filtered_assembly \
    --fwd-reads $fwd_reads \
    --rev-reads $rev_reads \
    --out $coverages\
    --cpus $cpus

# Step 3: Annotate and filter markers
# input:
# $orfs --> User input
# $cpus --> User input
# $seed --> User input

orfs="${outdir}/${simpleName}.orfs.faa"
orfs_nuc="${outdir}/${simpleName}.orfs.fna"

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
    hmmscan="${outdir}/${simpleName}.${kingdom}.hmmscan.tsv"
    markers="${outdir}/${simpleName}.${kingdom}.markers.tsv"

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
blast="${outdir}/${simpleName}.blastp.tsv" #Generate output file name

diamond blastp \
    --query $orfs \
    --db "$ncbi/nr.dmnd" \
    --evalue 1e-5 \
    --max-target-seqs 200 \
    --threads $cpus \
    --outfmt 6 \
    --out $blast

# output:
lca="${outdir}/${simpleName}.orfs.lca.tsv"
sseqid2taxid="${outdir}/${simpleName}.orfs.sseqid2taxid.tsv"
errorTaxids="${outdir}/${simpleName}.orfs.errortaxids.tsv"

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
votes="${outdir}/${simpleName}.taxids.tsv"

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
# e.g. ${outdir}/${simpleName}.bacteria.fna
# e.g. ${outdir}/${simpleName}.archaea.fna
# e.g. ${outdir}/${simpleName}.taxonomy.tsv

# script:
autometa-taxonomy \
    --votes $votes \
    --output $outdir \
    --prefix $simpleName \
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
    fasta="${outdir}/${simpleName}.${kingdom}.fna"

    # kingdom-specific output:
    counts="${outdir}/${simpleName}.${kingdom}.${kmer_size}mers.tsv"
    normalized="${outdir}/${simpleName}.${kingdom}.${kmer_size}mers.${norm_method}.tsv"
    embedded="${outdir}/${simpleName}.${kingdom}.${kmer_size}mers.${norm_method}.${embed_method}.tsv"

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
taxonomy="${outdir}/${simpleName}.taxonomy.tsv" # Generated by step 4.3
# $gc_content --> Generated by step 1

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do
    # kingdom-specific input:
    kmers="${outdir}/${simpleName}.${kingdom}.${kmer_size}mers.${norm_method}.${embed_method}.tsv" # Generated by step 5
    markers="${outdir}/${simpleName}.${kingdom}.markers.tsv" # Generated by step 3

    # kingdom-specific output:
    output_binning="${outdir}/${simpleName}.${kingdom}.${cluster_method}.tsv"
    output_main="${outdir}/${simpleName}.${kingdom}.${cluster_method}.main.tsv"

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
    binning_main="${outdir}/${simpleName}.${kingdom}.${cluster_method}.main.tsv" # Generated by step 6
    markers="${outdir}/${simpleName}.${kingdom}.markers.tsv" # Generated by step 3

    # kingdom-specific output:
    output_stats="${outdir}/${simpleName}_${kingdom}_metabin_stats.tsv"
    output_taxonomy="${outdir}/${simpleName}_${kingdom}_metabin_taxonomy.tsv"
    output_metabins="${outdir}/${simpleName}_${kingdom}_metabins"

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

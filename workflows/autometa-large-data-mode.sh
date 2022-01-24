#!/usr/bin/env bash
#SBATCH -p partition
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

## First create environment to run Autometa, (and optionally GTDB-Tk and CheckM)
# git clone git@github.com:KwanLab/Autometa
# cd Autometa
# make create_environment
# conda activate autometa

## Install autometa from source
## NOTE: For list of available commands try `make` in the Autometa directory
# cd Autometa
# make install
# hmmpress -f autometa/databases/markers/bacteria.single_copy.hmm
# hmmpress -f autometa/databases/markers/archaea.single_copy.hmm

## Install GTDB-Tk for post-processing autometa bins
## For more info see: https://ecogenomics.github.io/GTDBTk/

## Install CheckM for post-processing autometa bins
## For more info see: https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm

# Filepaths
assembly="Path to metagenome assembly fasta file"
bam="Path to metagenome read alignments.bam"
orfs="Path to orfs used as input to diamond blastp"
blast="Path to diamond blastp results (outfmt 6)"
ncbi="Path to NCBI databases directory" # should contain prot.accession2taxid.gz, names.dmp, nodes.dmp and merged.dmp
markers_dbdir="Path to autometa markers databases" # should contain marker hmms and cutoffs files (can be downloaded from Autometa repo under `Autometa/autometa/databases/markers`)

# Autometa Parameters
length_cutoff=3000 # in bp
# kmer counting, normalization transform and embedding method
kmer_size=5
norm_method="am_clr" # choices: am_clr, clr, ilr
pca_dimensions=50 # NOTE: must be greater than $embed_dimensions
embed_dimensions=2 # NOTE: must be less than $pca_dimensions
embed_method="umap" # choices: bhsne, sksne, umap, densne, trimap
# Binning clustering method
clustering_method="hdbscan" # choices: hdbscan, dbscan
# Binning metrics cutoffs
completeness=20.0
purity=95.0
cov_stddev_limit=25.0
gc_stddev_limit=5.0
max_partition_size=10000
cpus=16
seed=42

# Step 0: Do some path handling with the provided `assembly` filepath
simpleName="TemplateAssemblyName"
outdir="AutometaOutdir"
if [ ! -d $outdir ]
then mkdir -p $outdir
fi

######### BEGIN #########

# Step 1: filter assembly by length and retrieve contig lengths as well as GC content

# input:
# $assembly --> User input
# $length_cutoff --> User input

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

# input:
# NOTE: $bam is defined at top and the rest of the inputs are generated by autometa

# output:
bed="${outdir}/${simpleName}.coverages.bed.tsv"
coverages="${outdir}/${simpleName}.coverages.tsv"

# script:
autometa-bedtools-genomecov --ibam $bam --bed $bed --output $coverages

# Step 3: Annotate and filter markers
# input:
# $orfs --> User input
# $cpus --> User input
# $seed --> User input
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
        --dbdir $markers_dbdir \
        --out $markers \
        --kingdom $kingdom \
        --parallel \
        --cpus 4 \
        --seed $seed
done

# Step 4.1: Determine ORF lowest common ancestor (LCA) amongst top hits

# input:
# $blast --> User Input
# $ncbi --> User Input

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
# e.g. ${outdir}/${prefix}.bacteria.fna
# e.g. ${outdir}/${prefix}.archaea.fna
# e.g. ${outdir}/${prefix}.taxonomy.tsv

# script:
autometa-taxonomy \
    --votes $votes \
    --output $outdir \
    --prefix $simpleName \
    --split-rank-and-write superkingdom \
    --assembly $assembly \
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
    if [ ! -f $fasta ]
    then
        echo "${fasta} does not exist, skipping..."
        continue
    fi

    # kingdom-specific output:
    counts="${outdir}/${simpleName}.${kingdom}.${kmer_size}mers.tsv"

    # script:
    autometa-kmers \
        --fasta $fasta \
        --kmers $counts \
        --size $kmer_size \
        --cpus $cpus
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
    kmers="${outdir}/${simpleName}.${kingdom}.${kmer_size}mers.tsv" # Generated by step 5 (counts)
    markers="${outdir}/${simpleName}.${kingdom}.markers.tsv" # Generated by step 3

    # kingdom-specific output:
    cache="${outdir}/${simpleName}_${kingdom}_cache"
    output_binning="${outdir}/${simpleName}.${kingdom}.${clustering_method}.tsv"
    output_main="${outdir}/${simpleName}.${kingdom}.${clustering_method}.main.tsv"

    if [ ! -f $kmers ]
    then
        echo "${kmers} does not exist, skipping..."
        continue
    fi

    if [ -f $output_main ] && [ -s $output_main ];then
        echo "$(basename $output_main) already exists. continuing..."
        continue
    fi

    # script:
    autometa-ldm-binning \
        --kmers $kmers \
        --coverages $coverages \
        --gc-content $gc_content \
        --markers $markers \
        --taxonomy $taxonomy \
        --output-binning $output_binning \
        --output-main $output_main \
        --clustering-method $clustering_method \
        --completeness $completeness \
        --purity $purity \
        --cov-stddev-limit $cov_stddev_limit \
        --gc-stddev-limit $gc_stddev_limit \
        --norm-method $norm_method \
        --pca-dims $pca_dimensions \
        --embed-method $embed_method \
        --embed-dims $embed_dimensions \
        --max-partition-size $max_partition_size \
        --starting-rank superkingdom \
        --cache $cache \
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
    binning_main="${outdir}/${simpleName}.${kingdom}.${clustering_method}.main.tsv" # Generated by step 6
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
        --metagenome $assembly \
        --ncbi $ncbi \
        --output-stats $output_stats \
        --output-taxonomy $output_taxonomy \
        --output-metabins $output_metabins
done

# Step 8: (OPTIONAL) Annotate Autometa bins' taxonomies and completeness/contamination metrics using GTDB-Tk & CheckM

# input:
# $output_metabins --> Generated by step 7

# kingdoms=(bacteria archaea)

# for kingdom in ${kingdoms[@]};do

#     # kingdom-specific input
#     output_metabins="${outdir}/${simpleName}_${kingdom}_metabins"

#     # kingdom-specific outputs
#     gtdbtk_outdir="${outdir}/${simpleName}_${kingdom}_gtdbtk_classify_wf"
#     checkm_outdir="${outdir}/${simpleName}_${kingdom}_checkm_lineage_wf"
#     checkm_stats="${outdir}/${simpleName}_${kingdom}_checkm_stats.tsv"

#     if [ ! -d $output_metabins ]
#     then
#         echo "${output_metabins} does not exist, skipping..."
#         continue
#     fi

#     # NOTE: autometa-binning-summary writes out all bins with the <cluster>.fna extension
#     # Unclustered sequences are written to one file with the .fasta extension. e.g. outdir/unclustered.fasta

#     # First run `gtdbtk check_install` to ensure the following command will run
#     gtdbtk classify_wf \
#         --genome_dir $output_metabins \
#         --extension fna \
#         --out_dir $gtdbtk_outdir \
#         --cpus $cpus \
#         --scratch_dir $LOCAL \
#         --pplacer_cpus 1

#     checkm lineage_wf \
#         --extension fna \
#         --threads $cpus \
#         --tmpdir $LOCAL \
#         --pplacer_threads 1 \
#         --tab_table \
#         --file $checkm_stats \
#         $output_metabins \
#         $checkm_outdir

# done

#########  END  #########
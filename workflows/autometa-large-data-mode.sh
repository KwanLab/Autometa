#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --error=autometa.%J.err
#SBATCH --output=autometa.%J.out

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
ncbi="Path to NCBI databases directory" # For more info see: https://autometa.readthedocs.io/en/latest/databases.html#ncbi
gtdb="Path to GTDB databases directory" # (Optional) For more info see: https://autometa.readthedocs.io/en/latest/databases.html#genome-taxonomy-database-gtdb
markers_dbdir="Path to autometa markers databases" # should contain marker hmms and cutoffs files (can be downloaded from Autometa repo under `Autometa/autometa/databases/markers`)

# Autometa Parameters
length_cutoff=3000 # in bp
# Taxon Assignment Parameters
taxa_routine="ncbi" # Choices are "ncbi" or "ncbi_gtdb"
# NOTE: When using the "ncbi_gtdb" option, blastP will be performed against the GTDB database
# using the kingdom-specific ORFs retrieved from the NCBI taxon-assignment sub-workflow.
# K-mer Counting, Normalization and Embedding Parameters
kmer_size=5
norm_method="am_clr" # choices: am_clr, clr, ilr
pca_dimensions=50 # NOTE: must be greater than $embed_dimensions
embed_dimensions=2 # NOTE: must be less than $pca_dimensions
embed_method="umap" # choices: bhsne, sksne, umap, densmap, trimap
# Binning Parameters (clustering methods and MAG quality thresholds)
cluster_method="hdbscan" # choices: hdbscan, dbscan
# Binning metrics cutoffs
completeness=20.0 # Accept MAGs greater than this value
purity=95.0 # Accept MAGs greater than this value
cov_stddev_limit=25.0 # Accept MAGs less than this value
gc_stddev_limit=5.0 # Accept MAGs less than this value
max_partition_size=10000
# Runtime Parameters
cpus=16
seed=42

if [[ $taxa_routine != "ncbi" ]] && [[ $taxa_routine != "ncbi_gtdb" ]]
then
    echo "ERROR: Invalid Taxonomic routine value. Please choose between ncbi or ncbi_gtdb. Current selection: ${taxa_routine}"
    exit 1
fi

# Step 0: Do some path handling with the provided `assembly` filepath
simple_name="TemplateAssemblyName"
outdir="AutometaOutdir"
if [ ! -d $outdir ]
then mkdir -p $outdir
fi

######### BEGIN #########

# Step 00: Report autometa version
set -x
autometa --version
{ set +x; } 2>/dev/null

# Step 1: filter assembly by length and retrieve contig lengths as well as GC content

# input:
# $assembly --> User input
# $length_cutoff --> User input

# output:
filtered_assembly="${outdir}/${simple_name}.filtered.fna"
gc_content="${outdir}/${simple_name}.gc_content.tsv"

# script:
set -x
autometa-length-filter \
    --assembly $assembly \
    --cutoff $length_cutoff \
    --output-fasta $filtered_assembly \
    --output-gc-content $gc_content
{ set +x; } 2>/dev/null

# Step 2: Determine coverages from assembly read alignments

# input:
# NOTE: $bam is defined at top and the rest of the inputs are generated by autometa

# output:
bed="${outdir}/${simple_name}.coverages.bed.tsv"
coverages="${outdir}/${simple_name}.coverages.tsv"

# script:
set -x
autometa-bedtools-genomecov --ibam $bam --bed $bed --output $coverages
{ set +x; } 2>/dev/null

# Step 3: Annotate and filter markers
# input:
# $orfs --> User input
# $cpus --> User input
# $seed --> User input
kingdoms=(bacteria archaea)

# NOTE: We iterate through both sets of markers for binning both bacterial and archeal kingdoms
for kingdom in ${kingdoms[@]};do
    # kingdom-specific output:
    hmmscan="${outdir}/${simple_name}.${kingdom}.hmmscan.tsv"
    markers="${outdir}/${simple_name}.${kingdom}.markers.tsv"

    # script:
    set -x
    autometa-markers \
        --orfs $orfs \
        --hmmscan $hmmscan \
        --dbdir $markers_dbdir \
        --out $markers \
        --kingdom $kingdom \
        --parallel \
        --cpus 4 \
        --seed $seed
    { set +x; } 2>/dev/null
done

# Step 4.1: Determine ORF lowest common ancestor (LCA) amongst top hits

# input:
# $blast --> User Input
# $ncbi --> User Input
# $dbtype --> Updated according to $taxa_routine
dbtype="ncbi"
prefix="${simple_name}.${dbtype}"

# output:
lca="${outdir}/${prefix}.orfs.lca.tsv"
sseqid_to_taxid="${outdir}/${prefix}.orfs.sseqid2taxid.tsv"
error_taxids="${outdir}/${prefix}.orfs.errortaxids.tsv"

# script:
set -x
autometa-taxonomy-lca \
    --blast $blast \
    --dbdir $ncbi \
    --dbtype $dbtype \
    --lca-output $lca \
    --sseqid2taxid-output $sseqid_to_taxid \
    --lca-error-taxids $error_taxids
{ set +x; } 2>/dev/null

# Step 4.2: Perform Modified Majority vote of ORF LCAs for all contigs that returned hits in blast search

# input:
# $lca --> Generated by step 4.1
# $ncbi --> User Input

# output:
votes="${outdir}/${prefix}.taxids.tsv"

# script:
set -x
autometa-taxonomy-majority-vote --lca $lca --output $votes --dbdir $ncbi --dbtype $dbtype
{ set +x; } 2>/dev/null
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
set -x
autometa-taxonomy \
    --votes $votes \
    --output $outdir \
    --prefix $prefix \
    --split-rank-and-write superkingdom \
    --assembly $assembly \
    --dbdir $ncbi \
    --dbtype $dbtype
{ set +x; } 2>/dev/null

# Step 5: Taxon-assignment using the GTDB database 
# NOTE: only performed if `taxa_routine` is 'ncbi_gtdb'

# Step 5.1: Extract bacterial ORFs and run GTDB
# input:
# $kingdom_fasta --> Generated by step 4.3
# $orfs --> User Input

# output:
# contig_ids --> text file containing metagenome contig IDs classified within NCBI bacteria and archaea
# orf_ids --> text file containing contig ORF IDs classified within NCBI bacteria and archaea
# kingdom_orfs --> fasta file containing metagenome ORFs classified within NCBI bacteria or archaea
# gtdb_orfs --> metagenome orfs classified within NCBI bacteria *and* archaea

if [[ "$taxa_routine" == "ncbi_gtdb" ]]
then
    echo "Running GTDB taxon assignment step."

    # output
    gtdb_orfs="${outdir}/${prefix}.orfs.faa"

    for kingdom in ${kingdoms[@]};do

        kingdom_fasta="${outdir}/${prefix}.${kingdom}.fna"
        
        contig_ids="${outdir}/${prefix}.${kingdom}.contigIDs.txt"
        orf_ids="${outdir}/${prefix}.${kingdom}.orfIDs.txt"
        kingdom_orfs="${outdir}/${prefix}.${kingdom}.orfs.faa"

        if [ ! -f $kingdom_fasta ]
        then
            echo "${kingdom_fasta} does not exist, skipping..."
            continue
        fi

        # Retrieve contig IDs from kingdom fasta file
        set -x
        grep ">" $kingdom_fasta | \
            sed 's/^>//' | \
            sed 's/$/_/' | \
            cut -f1 -d" " > $contig_ids
        # Retrieve ORF IDs from contig IDs
        grep -f $contig_ids $orfs | sed 's/^>//' | cut -f1 -d" " > $orf_ids
        # Retrieve ORF seqs from ORF IDs
        seqkit grep \
            --pattern-file $orf_ids \
            --out-file $kingdom_orfs \
            $orfs
        # Concatenate kingdom ORFs to single file for GTDB blastp
        cat $kingdom_orfs >> $gtdb_orfs
        { set +x; } 2>/dev/null
    done
    dbtype="gtdb"
    prefix="${simple_name}.${dbtype}"

    # Step 5.2: Run blastp
    # input:
    # $gtdb_orfs --> Generated from step 5.1
    gtdb_dmnd_db="${gtdb}/gtdb.dmnd" # generated using autometa-setup-gtdb (Must be performed prior to using this script)
    # output
    blast="${outdir}/${prefix}.blastp.tsv"

    # script
    set -x
    diamond blastp \
        --query $gtdb_orfs \
        --db $gtdb_dmnd_db \
        --evalue 1e-5 \
        --max-target-seqs 200 \
        --threads $cpus \
        --outfmt 6 \
        --out $blast
    { set +x; } 2>/dev/null

    # Step 5.3: Determine LCA
    # input:
    # $blast --> Generated from step 5.2
    # $gtdb --> User Input
    # $dbtype --> Updated according to $taxa_routine
    dbtype="gtdb"

    # output:
    lca="${outdir}/${prefix}.orfs.lca.tsv"
    sseqid_to_taxid="${outdir}/${prefix}.orfs.sseqid2taxid.tsv"
    error_taxids="${outdir}/${prefix}.orfs.errortaxids.tsv"

    # script:
    set -x
    autometa-taxonomy-lca \
        --blast $blast \
        --dbdir $gtdb \
        --dbtype $dbtype \
        --lca-output $lca \
        --sseqid2taxid-output $sseqid_to_taxid \
        --lca-error-taxids $error_taxids 
    { set +x; } 2>/dev/null

    # Step 5.4: Perform Modified Majority vote of ORF LCAs for all contigs that returned hits in blast search
    # input:
    # $lca --> Generated from step 5.3
    # $gtdb --> User Input

    # output:
    votes="${outdir}/${prefix}.taxids.tsv"

    # script:
    set -x
    autometa-taxonomy-majority-vote \
        --lca $lca \
        --output $votes \
        --dbdir $gtdb \
        --dbtype $dbtype
    { set +x; } 2>/dev/null

    # Step 5.5: Split assigned taxonomies into kingdoms
    # input:
    # $votes --> Generated from step 5.4
    # $outdir --> Generated from step 0
    # $prefix --> Generated as input from steps 5.2 to 6
    # $filtered_assembly --> Generated from step 1
    # $gtdb --> User Input

    # output:
    # Will write recovered superkingdoms to $outdir
    # e.g. ${outdir}/${prefix}.bacteria.fna
    # e.g. ${outdir}/${prefix}.archaea.fna
    # e.g. ${outdir}/${prefix}.taxonomy.tsv

    # script:
    set -x
    autometa-taxonomy \
        --votes $votes \
        --output $outdir \
        --prefix $prefix \
        --split-rank-and-write superkingdom \
        --assembly $filtered_assembly \
        --dbdir $gtdb \
        --dbtype $dbtype
    { set +x; } 2>/dev/null
fi

# Step 6: Perform k-mer counting on respective kingdoms

# input:
# $kmer_size --> User input
# $cpus --> User input

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do
    # kingdom-specific input:
    # NOTE: $fasta --> Generated by step 4.3
    fasta="${outdir}/${prefix}.${kingdom}.fna" # NOTE: $prefix is updated according to taxa_routine above
    # kingdom-specific output:
    counts="${outdir}/${prefix}.${kingdom}.${kmer_size}mers.tsv"
    
    if [ ! -f $fasta ]
    then
        echo "${fasta} does not exist, skipping..."
        continue
    fi
    # script:
    set -x
    autometa-kmers \
        --fasta $fasta \
        --kmers $counts \
        --size $kmer_size \
        --cpus $cpus
    { set +x; } 2>/dev/null
done

# Step 7: Perform binning on each set of bacterial and archaeal contigs

# input:
# $coverages --> Generated by step 2
# $gc_content --> Generated by step 1
taxonomy="${outdir}/${prefix}.taxonomy.tsv" # NOTE: $prefix is updated according to taxa_routine above
# $taxonomy is generated by either steps 4.3 or 5.5 depending on whether taxa_routine is 'ncbi' or 'ncbi_gtdb', respectively
# $cluster_method --> User input
# $completeness --> User input
# $purity --> User input
# $cov_stddev_limit --> User input
# $gc_stddev_limit --> User input
# $norm_method --> User input
# $pca_dimensions --> User input
# $embed_method --> User input
# $max_partition_size --> User input

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do
    # kingdom-specific input:
    kmers="${outdir}/${prefix}.${kingdom}.${kmer_size}mers.tsv" # Generated by step 6 (counts)
    markers="${outdir}/${simple_name}.${kingdom}.markers.tsv" # Generated by step 3 (before taxon-assignment sub-workflows)

    # kingdom-specific output:
    cache="${outdir}/${prefix}_${dbtype}_${kingdom}_cache"
    output_binning="${outdir}/${prefix}.${kingdom}.${cluster_method}.tsv"
    output_main="${outdir}/${prefix}.${kingdom}.${cluster_method}.main.tsv"

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
    set -x
    autometa-binning-ldm \
        --kmers $kmers \
        --coverages $coverages \
        --gc-content $gc_content \
        --markers $markers \
        --taxonomy $taxonomy \
        --output-binning $output_binning \
        --output-main $output_main \
        --clustering-method $cluster_method \
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
    { set +x; } 2>/dev/null
done

# Step 8: Create binning summary files

# input:
# $ncbi -> User input
# $gtdb -> User input
# $assembly -> User input
# $dbtype -> # NOTE: $prefix is updated according to taxa_routine above

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do

    # kingdom-specific input:
    binning_main="${outdir}/${prefix}.${kingdom}.${cluster_method}.main.tsv" # Generated by step 7
    markers="${outdir}/${simple_name}.${kingdom}.markers.tsv" # Generated by step 3

    # kingdom-specific output:
    output_stats="${outdir}/${prefix}_${kingdom}_metabin_stats.tsv"
    output_taxonomy="${outdir}/${prefix}_${kingdom}_metabin_taxonomy.tsv"
    output_metabins="${outdir}/${prefix}_${kingdom}_metabins"

    if [ ! -f $binning_main ]
    then
        echo "${binning_main} does not exist, skipping..."
        continue
    fi

    if [[ "$taxa_routine" == "ncbi_gtdb" ]]
    then
        dbdir=$gtdb
    else
        dbdir=$ncbi
    fi
    set -x
    autometa-binning-summary \
        --binning-main $binning_main \
        --markers $markers \
        --metagenome $assembly \
        --dbdir $dbdir \
        --dbtype $dbtype \
        --output-stats $output_stats \
        --output-taxonomy $output_taxonomy \
        --output-metabins $output_metabins
    { set +x; } 2>/dev/null
done

# Step 9: (OPTIONAL) Annotate Autometa bins' taxonomies and completeness/contamination metrics using GTDB-Tk & CheckM

# input:
# $output_metabins --> Generated by step 8

# kingdoms=(bacteria archaea)

# for kingdom in ${kingdoms[@]};do

#     # kingdom-specific input
#     # output_metabins="${outdir}/${simple_name}_${dbtype}_${kingdom}_metabins"

#     # kingdom-specific outputs
#     gtdbtk_outdir="${outdir}/${simple_name}_${dbtype}_${kingdom}_gtdbtk_classify_wf"
#     checkm_outdir="${outdir}/${simple_name}_${dbtype}_${kingdom}_checkm_lineage_wf"
#     checkm_stats="${outdir}/${simple_name}_${dbtype}_${kingdom}_checkm_stats.tsv"

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

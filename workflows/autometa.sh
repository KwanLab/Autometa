#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --error=autometa.%J.err
#SBATCH --output=autometa.%J.out

# NOTE: To create the conda environment for autometa you can supply the Makefile command:
# make create_environment

# Now activate the created conda env
# conda activate autometa

# NOTE: To install autometa in the created conda environment you can supply the Makefile command:
# make install

# Filepaths
assembly="Path to metagenome assembly fasta file"
bam="Path to metagenome read alignments.bam"
orfs="Path to orfs used as input to diamond blast"
blast="Path to diamond output file (outfmt 6)." # BlastP should be done against the NCBI `nr` database.
ncbi="Path to NCBI databases directory"
gtdb="Path to GTDB databases directory"

# Autometa Parameters
length_cutoff=3000 # in bp
# kmer counting, normalization transform and embedding method
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
cpus=16
seed=42
taxa_routine="Choose between ncbi or ncbi_gtdb" # Choices are "ncbi" or "ncbi_gtdb". When using the ncbi_gtdb option, blastP will be repeated against the GTDB database using the kingdom specific ORFs retrieved taxaon assignment after blastP against the `nr` database.

# Step 0: Do some Path handling with the provided `assembly` filepath
simpleName="TemplateAssemblyName"
outdir="AutometaOutdir"
if [ ! -d $outdir ]
then mkdir -p $outdir
fi

######### BEGIN #########

# Step 00: Report autometa version
autometa --version

# Step 1: filter assembly by length and retrieve contig lengths as well as GC content

# input:
# $assembly --> User input
# $length_cutoff --> User input

# output:
filtered_assembly="${outdir}/${simpleName}.filtered.fna"
gc_content="${outdir}/${simpleName}.gc_content.tsv"

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
bed="${outdir}/${simpleName}.coverages.bed.tsv"
coverages="${outdir}/${simpleName}.coverages.tsv"

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
    hmmscan="${outdir}/${simpleName}.${kingdom}.hmmscan.tsv"
    markers="${outdir}/${simpleName}.${kingdom}.markers.tsv"

    # script:
    set -x
    autometa-markers \
        --orfs $orfs \
        --hmmscan $hmmscan \
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

# output:
lca="${outdir}/${simpleName}.orfs.lca.tsv"
sseqid2taxid="${outdir}/${simpleName}.orfs.sseqid2taxid.tsv"
errorTaxids="${outdir}/${simpleName}.orfs.errortaxids.tsv"

# script:
set -x
autometa-taxonomy-lca \
    --blast $blast \
    --dbdir $ncbi \
    --lca-output $lca \
    --sseqid2taxid-output $sseqid2taxid \
    --lca-error-taxids $errorTaxids \
    --dbtype ncbi
{ set +x; } 2>/dev/null

# Step 4.2: Perform Modified Majority vote of ORF LCAs for all contigs that returned hits in blast search

# input:
# $lca --> Generated by step 4.1
# $ncbi --> User Input

# output:
votes="${outdir}/${simpleName}.taxids.tsv"

# script:
set -x
autometa-taxonomy-majority-vote --lca $lca --output $votes --dbdir $ncbi --dbtype ncbi
{ set +x; } 2>/dev/null

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
set -x
autometa-taxonomy \
    --votes $votes \
    --output $outdir \
    --prefix $simpleName \
    --split-rank-and-write superkingdom \
    --assembly $assembly \
    --dbdir $ncbi \
    --dbtype ncbi
{ set +x; } 2>/dev/null

# Run the taxonomy steps using GTDB database
# Step 5.1: Extract bacterial ORFs and run GTDB


gtdbOrfs="${outdir}/${simpleName}.orfs.gtdb.faa"

if [[ "$taxa_routine" == "ncbi_gtdb" ]]
then
    echo "Running GTDB taxon assignment step."

    for kingdom in ${kingdoms[@]};do

        contigIds="${outdir}/${simpleName}.${kingdom}.contigIDs.txt"
        orfIds="${outdir}/${simpleName}.${kingdom}.orfIDs.txt"
        kingdomFasta="${outdir}/${simpleName}.${kingdom}.fna"
        kingdomOrfs="${outdir}/${simpleName}.${kingdom}.orfs.faa"

        if [ ! -f $kingdomFasta ]
        then
            echo "${kingdomFasta} does not exist, skipping..."
            continue
        fi

        # Retrieve contig IDs from kingdom fasta file
        set -x
        grep ">" $kingdomFasta | \
        sed 's/^>//' | \
        sed 's/$/_/' | \
        cut -f1 -d " " > $contigIds
        # Retrieve ORF IDs from contig IDs
        grep -f $contigIds $orfs |  sed 's/^>//' | cut -f1 -d " " > $orfIds
        # Retrieve ORF seqs from ORF IDs
        seqkit grep --pattern-file $orfIds $orfs --out-file $kingdomOrfs
        # Concatenate kingdom ORFs to single file for GTDB blastp
        cat $kingdomOrfs >> $gtdbOrfs
        { set +x; } 2>/dev/null
    done

    #Step 5.2: Run blastp
    blast="${outdir}/${simpleName}.blastp.gtdb.tsv" #Generate output file name

    set -x
    diamond blastp \
        --query $gtdbOrfs \
        --db "$gtdb/gtdb.dmnd" \
        --evalue 1e-5 \
        --max-target-seqs 200 \
        --threads $cpus \
        --outfmt 6 \
        --out $blast
    { set +x; } 2>/dev/null

    # output:
    lca="${outdir}/${simpleName}.orfs.lca.gtdb.tsv"
    sseqid2taxid="${outdir}/${simpleName}.orfs.sseqid2taxid.gtdb.tsv"
    error_taxids="${outdir}/${simpleName}.orfs.errortaxids.gtdb.tsv"

    # script:
    set -x
    autometa-taxonomy-lca \
        --blast $blast \
        --dbdir $gtdb \
        --dbtype gtdb \
        --lca-output $lca \
        --sseqid2taxid-output $sseqid2taxid \
        --lca-error-taxids $error_taxids 
    { set +x; } 2>/dev/null

    # Step 5.3: Perform Modified Majority vote of ORF LCAs for all contigs that returned hits in blast search

    votes="${outdir}/${simpleName}.taxids.gtdb.tsv"
    set -x
    autometa-taxonomy-majority-vote \
        --lca $lca \
        --output $votes \
        --dbdir $gtdb \
        --dbtype gtdb
    { set +x; } 2>/dev/null

    # Step 5.4: Split assigned taxonomies into kingdoms
    set -x
    autometa-taxonomy \
        --votes $votes \
        --output "${outdir}/gtdb_taxa/" \
        --prefix $simpleName \
        --split-rank-and-write superkingdom \
        --assembly $filtered_assembly \
        --dbdir $gtdb \
        --dbtype gtdb
    { set +x; } 2>/dev/null
fi

# Step 6: Perform k-mer counting on respective kingdoms

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
    if [[ "$taxa_routine" == "ncbi_gtdb" ]]
    then
        fasta="${outdir}/gtdb_taxa/${simpleName}.${kingdom}.fna"
    else
        fasta="${outdir}/${simpleName}.${kingdom}.fna"
    fi
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
    set -x
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
    { set +x; } 2>/dev/null
done

# Step 7: Perform binning on each set of bacterial and archaeal contigs

# input:
# $cpus --> User input
# $seed --> User input
# $gc_content --> Generated by step 1



if [[ "$taxa_routine" == "ncbi_gtdb" ]]
then
    taxonomy="${outdir}/gtdb_taxa/${simpleName}.taxonomy.tsv"
else 
    taxonomy="${outdir}/${simpleName}.taxonomy.tsv"
fi

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do
    # kingdom-specific input:
    kmers="${outdir}/${simpleName}.${kingdom}.${kmer_size}mers.${norm_method}.${embed_method}.tsv" # 
    markers="${outdir}/${simpleName}.${kingdom}.markers.tsv" # Generated by step 3

    # kingdom-specific output:
    output_binning="${outdir}/${simpleName}.${kingdom}.${cluster_method}.tsv"
    output_main="${outdir}/${simpleName}.${kingdom}.${cluster_method}.main.tsv"

    if [ -f $output_main ] && [ -s $output_main ];then
        echo "$(basename $output_main) already exists. continuing..."
        continue
    fi

    if [ ! -f $kmers ]
        then echo "Required ${kingdom} files not found, skipping."
        continue
    fi

    # script:
    set -x
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
    { set +x; } 2>/dev/null
done

# Step 8: Create binning summary files

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

    if [[ "$taxa_routine" == "ncbi_gtdb" ]]
    then
        set -x
        autometa-binning-summary \
            --binning-main $binning_main \
            --markers $markers \
            --metagenome $assembly \
            --dbdir $gtdb \
            --dbtype gtdb \
            --output-stats $output_stats \
            --output-taxonomy $output_taxonomy \
            --output-metabins $output_metabins
        { set +x; } 2>/dev/null
    else    
        set -x
        autometa-binning-summary \
            --binning-main $binning_main \
            --markers $markers \
            --metagenome $assembly \
            --dbdir $ncbi \
            --dbtype ncbi \
            --output-stats $output_stats \
            --output-taxonomy $output_taxonomy \
            --output-metabins $output_metabins
        { set +x; } 2>/dev/null
    fi
done

#########  END  #########

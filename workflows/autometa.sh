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
ncbi="Path to NCBI databases directory" # For more info see: https://autometa.readthedocs.io/en/latest/databases.html#ncbi
gtdb="Path to GTDB databases directory" # (Optional) For more info see: https://autometa.readthedocs.io/en/latest/databases.html#genome-taxonomy-database-gtdb

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
embed_method="bhsne" # choices: bhsne, sksne, umap, densmap, trimap
# Binning Parameters (clustering methods and MAG quality thresholds)
cluster_method="hdbscan" # choices: hdbscan, dbscan
# Binning metrics cutoffs
completeness=20.0 # Accept MAGs greater than this value
purity=95.0 # Accept MAGs greater than this value
cov_stddev_limit=25.0 # Accept MAGs less than this value
gc_stddev_limit=5.0 # Accept MAGs less than this value
# Runtime Parameters
cpus=16
seed=42

if [[ $taxa_routine != "ncbi" ]] && [[ $taxa_routine != "ncbi_gtdb" ]]
then
    echo "ERROR: Invalid Taxonomic routine value. Please choose between ncbi or ncbi_gtdb. Current selection: ${taxa_routine}"
    exit 1
fi

# Step 0: Do some Path handling with the provided `assembly` filepath
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

# output:
lca="${outdir}/${simple_name}.orfs.lca.tsv"
sseqid_to_taxid="${outdir}/${simple_name}.orfs.sseqid2taxid.tsv"
error_taxids="${outdir}/${simple_name}.orfs.errortaxids.tsv"

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
# $dbtype --> Updated according to $taxa_routine

# output:
votes="${outdir}/${simple_name}.taxids.tsv"

# script:
set -x
autometa-taxonomy-majority-vote --lca $lca --output $votes --dbdir $ncbi --dbtype $dbtype
{ set +x; } 2>/dev/null

# Step 4.3: Split assigned taxonomies into kingdoms

# input:
# $votes --> Generated by step 4.2
# $outdir --> Generated by step 0
prefix="${simple_name}.${dbtype}"
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
    gtdb_orfs="${outdir}/${simple_name}.orfs.gtdb.faa"

    for kingdom in ${kingdoms[@]};do

        kingdom_fasta="${outdir}/${prefix}.${kingdom}.fna"
        
        contig_ids="${outdir}/${simple_name}.${kingdom}.contigIDs.txt"
        orf_ids="${outdir}/${simple_name}.${kingdom}.orfIDs.txt"
        kingdom_orfs="${outdir}/${simple_name}.${kingdom}.orfs.faa"

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

    # Step 5.2: Run blastp
    # input:
    # $gtdb_orfs --> Generated from step 5.1
    gtdb_dmnd_db="${gtdb}/gtdb.dmnd" # generated using autometa-setup-gtdb (Must be performed prior to using this script)
    # output
    blast="${outdir}/${simple_name}.blastp.gtdb.tsv"

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

    #Step 5.3: Determine LCA
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
    votes="${outdir}/${simple_name}.taxids.gtdb.tsv"

    # script:
    set -x
    autometa-taxonomy-majority-vote \
        --lca $lca \
        --output $votes \
        --dbdir $gtdb \
        --dbtype gtdb
    { set +x; } 2>/dev/null

    # Step 5.5: Split assigned taxonomies into kingdoms
    # input:
    # $votes --> Generated from step 5.4
    # $outdir --> Generated from step 0
    prefix="${simple_name}.${dbtype}"
    # $filtered_assembly --> Generated from step 1
    # $gtdb --> User Input

    # output:
    # Will write recovered superkingdoms to $outdir
    # e.g. ${outdir}/${simple_name}.gtdb.bacteria.fna
    # e.g. ${outdir}/${simple_name}.gtdb.archaea.fna
    # e.g. ${outdir}/${simple_name}.gtdb.taxonomy.tsv

    # script:
    set -x
    autometa-taxonomy \
        --votes $votes \
        --output $outdir \
        --prefix $prefix \
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
    fasta="${outdir}/${prefix}.${kingdom}.fna" # NOTE: $prefix is updated according to taxa_routine above
    counts="${outdir}/${prefix}.${kingdom}.${kmer_size}mers.tsv"
    normalized="${outdir}/${prefix}.${kingdom}.${kmer_size}mers.${norm_method}.tsv"
    embedded="${outdir}/${prefix}.${kingdom}.${kmer_size}mers.${norm_method}.${embed_method}.tsv"

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
taxonomy="${outdir}/${prefix}.taxonomy.tsv" # NOTE: $prefix is updated according to taxa_routine above
# $taxonomy is generated by either steps 4.3 or 5.5 depending on whether taxa_routine is 'ncbi' or 'ncbi_gtdb', respectively

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do
    # kingdom-specific input:
    kmers="${outdir}/${prefix}.${kingdom}.${kmer_size}mers.${norm_method}.${embed_method}.tsv" # Generated by step 6
    markers="${outdir}/${simple_name}.${kingdom}.markers.tsv" # Generated by step 3 (before taxon-assignment sub-workflows)

    # kingdom-specific output:
    output_binning="${outdir}/${prefix}.${kingdom}.${cluster_method}.tsv"
    output_main="${outdir}/${prefix}.${kingdom}.${cluster_method}.main.tsv"

    if [ -f $output_main ] && [ -s $output_main ];then
        echo "$(basename $output_main) already exists. continuing..."
        continue
    fi

    if [ ! -f $kmers ]
        then echo "${kingdom} file not found: kmers: ${kmers}), skipping."
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
# $gtdb -> User input
# $assembly -> User input
# $dbtype -> # NOTE: $prefix is updated according to taxa_routine above

kingdoms=(bacteria archaea)

for kingdom in ${kingdoms[@]};do

    # kingdom-specific input:
    binning_main="${outdir}/${prefix}.${kingdom}.${cluster_method}.main.tsv" # Generated by step 7
    markers="${outdir}/${simple_name}.${kingdom}.markers.tsv" # Generated by step 3

    # kingdom-specific output:
    output_stats="${outdir}/${simple_name}_${dbtype}_${kingdom}_metabin_stats.tsv"
    output_taxonomy="${outdir}/${simple_name}_${dbtype}_${kingdom}_metabin_taxonomy.tsv"
    output_metabins="${outdir}/${simple_name}_${dbtype}_${kingdom}_metabins"

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

#########  END  #########

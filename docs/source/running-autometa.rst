================
Running Autometa
================

Before running anything make sure you have activated the conda environment using
``conda activate autometa``.

See the :ref:`Install` page for details on setting up your conda environment.

Data preparation
================

Before you run Autometa, you need to have assembled your shotgun metagenome. The following workflow is recommended:

#. Trim adapter sequences from the reads. We prefer to use Trimmomatic_, but you can go ahead and any tool of your preference.
#. Quality check of reads to make sure that the adapters have been removed, we use FastQC_ for this. Also make sure that most reads are of good quality, if not you can trim the poor quality reads using Trimmomatic_.
#. Assemble the trimmed reads. We recommend using MetaSPAdes which is a part of the SPAdes_ package to assemble the trimmed reads but you can use any other assembler as well.
#. An optional thing to do here would be to check the quality of your assembly as well. This would give you N50 which could be useful in selecting the value of length-filter. We tend to use metaQuast_ for this (use ``--min-contig 1`` option to get an accurate N50).

.. note::

    If you use end up using SPAdes then Autometa can use the coverage information in the contig names. If you have used any other assembler, then you first have to make a coverage table.

    Fortunately, Autometa can construct this table for you with: ``python -m autometa.commmon.coverage``. Use ``--help`` to get the complete usage.

For the impatient (Quickstart)
==============================

You can run the entire workflow using nexflow. This will run all the different entrypoints for you and would restart from the lastest entrypoint in case of any error. 

.. code-block:: bash

    ./nextflow run autometa.nf

Step by step tutorial
=====================

Here is the step by step tutorial of the entire pipeline. This is helpful in case you have your own files or just want to run a specific step.

1. Length filter
----------------

The first step when running autometa is the legth filtering. This would remove any contigs that are below the length cutoff and would generate a new assembly without the filtered contigs. This is useful in removing the noise from the data, as small contigs may have ambigious kmer frequencies. The default cutoff if 3,000bp, ie. any contig that is smaller than 3,000bp would be removed.

.. note::
    It is important that you alter the cutoff based on your N50. If your N50 is really smaller, eg 500bp (pretty common for soil assemblies), then you might want to lower your cutoff to somehwere near N50.
    
Use the following command to run the length-filter step:

.. code-block:: bash

    autometa-length-filter <path/to/assembly.fasta> \
    <path/to/filtered/assembly.fasta> --cutoff <length_cutoff_value>

The above command would generate the following files:

- A new assembly without the contigs below the length filter called ``<path/to/filtered/assembly.fasta>``.

You can view the complete command line opions using ``autometa-length-filter -h``

2. K-mer counting
-----------------

A k-mer is just a sequence of k characters in a string (or nucleotides in a DNA sequence)(`ref <https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/>`_). It is known that contigs that belong to the same genome have similar k-mer composition. 

This step does the following:

#. Create a  k-mer matrix of 512 dimensions using 5-mer frequencies (ie. k-mer of size 5). k-mer size can be altered using the ``--size`` flag
#. Normaization of the k-mer matrix.
#. Reduce the dimensions of 5-mer frequencies to fifty using principle component analysis (PCA) and then to two using Barnes-Hut Stochastic Neighbor Embedding (BH-tSNE) for each contig. BH-tSNE is the default. Other embedding methods that are available are Uniform Manifold Approximation and Projection (UMAP) and SKSNE. This is done to allow the ease of visualization and manual binning of the contigs (see `ViZBin <https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-014-0066-1>`_ paper) 

If you think you need your final output to be in more than 2 dimensions you can alter it using ``--pca-dimensions`` and ``--embed-dimensions`` parameters. You can do all this using the following command:

Use the following command to run the k-mer counting step:

.. code-block:: bash

    autometa-kmers --fasta <path/to/filtered/assembly.fasta> \
    --kmers <path/to/kmer_table.tsv> \
    --normalized <path/to/kmers_norm.tsv> \
    --embedded <path/to/embedded_kmers.tsv> \
    --size <kmer_size> --norm-method <normalization_method> \
    --do-pca --pca-dimensions <pca_dimensions> \
    --embed-dimensions <final_dimensions> \
    --embed-method <kmer_embedding_method> --multiprocess \
    --cpus <num_of_cpus_to_use>

The above command would generate the following files:

#. Raw k-mer matrix table called ``<path/to/kmer_table.tsv>``
#. A normalized k-mer matrix called ``<path/to/kmers_norm.tsv>``
#. An embedded k-mer matrix called ``<path/to/embedded_kmers.tsv>``

You can view the complete command line opions using ``autometa-kmers -h``

3. Coverage calculation
-----------------------

Coverage calculation for each contig is done to provide another parametere to use while clustering contigs. In case you have used SPades to assemble your metagenome, you can use the following command to generate the coverage table:

.. code-block:: bash

    autometa-coverage --from-spades --assembly <path/to/filtered/assembly.fasta> \ 
    --out <path/to/coverage_table.tsv>

In case you have assembled your metagenome using some other assembler you can use the any one of the following commands to generate the coverage table.

.. code-block:: bash

    # If you have the bed file
    autometa-coverage --assembly <path/to/filtered/assembly.fasta> \ 
    --bed <path/to/alignments.bed> --out <path/to/coverage_table.tsv>

    # If you have the bam file
    autometa-coverage --assembly <path/to/filtered/assembly.fasta> \ 
    --bed <path/to/alignments.bam> --out <path/to/coverage_table.tsv>

    # If you have the sam file
    autometa-coverage --assembly <path/to/filtered/assembly.fasta> \ 
    --bed <path/to/alignments.sam> --out <path/to/coverage_table.tsv>

    # If you have the forward and reverse reads
    autometa-coverage --assembly <path/to/filtered/assembly.fasta> \ 
    --fwd-reads <path/to/fwd_reads.fq> --rev-reads <path/to/rev_reads.fq> \
    --out <path/to/coverage_table.tsv>

    # In case you have multiple fwd and rev reads. There should be no space between different read pairs
    autometa-coverage --assembly <path/to/filtered/assembly.fasta> \ 
    --fwd-reads <path/to/fwd_reads_1.fq>,<path/to/fwd_reads_2.fq> \ --rev-reads <path/to/rev_reads_1.fq>,<path/to/rev_reads_2.fq> \
    --out <path/to/coverage_table.tsv>

The above command would generate the following files:

- A coverage table having the coverage and length of each contig called ``<path/to/coverage_table.tsv>``

You can view the complete command line opions using ``autometa-coverage -h``

4. Taxonomy assignment
----------------------

This step does the following:

#. Identify genes in each contig with Prodigal.
#. Search gene protein sequences against nr with DIAMOND.
#. Determine the lowest common ancestor (LCA) of blast hits within 10% of the top bitscore.
#. Determine the taxonomy of each contig by examining the LCA of each component protein

We found that in host-associated metagenomes, this step vastly improves the binning performance of Autometa (and other pipelines) because less eukaryotic or viral contigs will be binned into bacterial bins. 

Use the following command to run the Taxonomy assignment step:

.. code-block:: bash

    autometa-taxonomy --assembly <path/to/filtered/assembly.fasta> \
    --nucl-orfs <path/to/nucleotide/ORFs.ffn> \
    --prot-orfs <path/to/amino_acid/ORFs.faa> \
    --blast <path/to/blastP.tsv> \
    --lca <path/to/lca.tsv> \
    --method majority_vote \
    --split-rank-and-write superkingdom \
    <path/to/taxonomy.tsv>

The above command would generate the following files:

#. Nucleotide sequence of the genes present in the given assembly called ``<path/to/nucleotide/ORFs.ffn>``
#. Amino acid sequence of the genes present in the given assembly called ``<path/to/amino_acid/ORFs.ffn>``
#. blastP output of the ORFs called ``<path/to/blastP.tsv>``
#. File having the lowest common ancestor of the ORFs called ``<path/to/lca.tsv>``
#. If --split-rank-and-write is specified then it will split contigs by provided canonical-rank column then write a file corresponding that rank. Eg. Bacteria.fasta, Archaea.fasta, etc for superkingdom.

You can view the complete command line opions using ``autometa-taxonomy -h``

5. Single copy markers
----------------------

Autometa uses single-copy markers to guide clustering, and does not assume that recoverable genomes will necessarily be ‘complete’.

Use the following command to run the assign single copy marker genes:

.. code-block:: bash

    # Archaeal markers
    --orfs <path/to/amino_acid/ORFs.faa> --kingdom archaea --hmmscan <path/to/archaea.hmmscan.tsv> \
    --out <path/to/archaea.markers.tsv>

    # Bacterial markers
    --orfs <path/to/amino_acid/ORFs.faa> --kingdom bacteria --hmmscan <path/to/bacteria.hmmscan.tsv> \
    --out <path/to/bacteria.markers.tsv>

The above command would generate the following files:

- Table having the archaeal and bacterial marker genes identified on each contig called ``<path/to/archaea.markers.tsv>`` and ``<path/to/bacteria.markers.tsv>`` for archaea and bacteria respectively

You can view the complete command line opions using ``autometa-markers -h``

6. Binning
----------

This is the step where contigs are binned into genomes. There are two binning algorithms to chose from Density-Based Spatial Clustering of Applications with Noise (DBSCAN) and Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN). The default is DBSCAN.

Autometa assesses clusters by examining both their completeness (number of expected single copy markers) and purity (number of single copy markers that are unique in the cluster).

If we supply a taxonomy table, then that is also used to help with clustering. Otherwise, Autometa clusters solely on 5-mer frequency and coverage. 

This step does the following:

#. Find single-copy marker genes in the input contigs with HMMER
#. Cluster contigs based on BH-tSNE coordinates (or any other embedding method that you have used), coverage and (optionally) taxonomy
#. Accept clusters that are estimated to be over 20% complete and 90% pure based on single-copy marker genes. These are default papameteres and can be altered to suit your needs.
#. Unclustered contigs leftover will be re-clustered until no more acceptable clusters are yielded

If you include a taxonomy table in the, Autometa will attempt to further partition the data based on ascending taxonomic specificity (i.e. in the order phylum, class, order, family, genus, species) when clustering unclustered contigs from a previous attempt. We found that this is mainly useful if you have a highly complex metagenome (lots of species), or you have several related species at similar coverage level.

Use the following command to run the binning:

.. code-block:: bash

    # Archaeal binning
    autometa-binning <path/to/kmers_norm.tsv> \
    <path/to/coverage_table.tsv> <path/to/archaea.markers.tsv> \
    <path/to/archaea_binning.tsv> --embedded-kmers <path/to/embedded_kmers.tsv> \
    --taxonomy <path/to/taxonomy.tsv> --clustering-method <dbscan or hdbscan> --domain archaea

    # Bacterial binning
    autometa-binning <path/to/kmers_norm.tsv> \
    <path/to/coverage_table.tsv> <path/to/bacterial.markers.tsv> \
    <path/to/bacteria_binning.tsv> --embedded-kmers <path/to/embedded_kmers.tsv> \
    --taxonomy <path/to/taxonomy.tsv> --clustering-method <dbscan or hdbscan> --domain bacteria

The above command would generate the following files:

- Final binning of each contig into a genome called ``<path/to/archaea_binning.tsv>`` and ``<path/to/bacteria_binning.tsv>`` for archaea and bacteria respectively

You can view the complete command line opions using ``autometa-binning -h``

7. Unclustered recruitment (Optional)
-------------------------------------

Supervised machine learning is used to classify the unclustered contigs to the bins that we have produced. This steop is optional and the results should be verified (see Note below) before going ahead with it.

.. note::
    The machine learning step has been seen to pick up contigs that not necessary belong to the genome. Careful inscpection of coverage and taxonomy should be done before you go ahead and use results from this step.

Use the following command to run the unclustered recruitment step:

.. code-block:: bash

    # Archaea
    autometa-unclustered-recruitment <path/to/kmers_norm.tsv> \
    <path/to/coverage_table.tsv> <path/to/archaea_binning.tsv> \
    <path/to/archaea.markers.tsv> <path/to/arachaea_unclustered_recruitment.tsv> \
    --taxonomy <path/to/taxonomy.tsv> --classifier decision_tree

    # Bacteria
    autometa-unclustered-recruitment <path/to/kmers_norm.tsv> \
    <path/to/coverage_table.tsv> <path/to/bacteria_binning.tsv> \
    <path/to/bacteria.markers.tsv> <path/to/bacteria_unclustered_recruitment.tsv> \
    --taxonomy <path/to/taxonomy.tsv> --classifier decision_tree

The above command would generate the following files:

- Recruitment of unclustered contig into a bins called ``<path/to/archaea.markers.tsv>`` and ``<path/to/bacteria_unclustered_recruitment.tsv>`` for archaea and bacteria respectively

You can view the complete command line opions using ``autometa-unclustered-recruitment -h``

Running modules
===============

Many of the Autometa modules may be run standalone.

Simply pass in the ``-m`` flag when calling a script to signify to python you are
running an Autometa *module*.

I.e. ``python -m autometa.common.kmers -h``

Running functions
=================

Many of the Autometa functions may be run standalone as well. It is same as importing any other python
function.

.. code-block:: python

    from autometa.common.external import samtools

    samtools.sort(sam=<path/to/sam/file>, out=<path/to/output/file>, nproc=4)


.. _SPAdes: http://cab.spbu.ru/software/spades/
.. _Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _metaQuast: http://quast.sourceforge.net/metaquast

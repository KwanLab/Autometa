=====================
Step by step tutorial
=====================

Here is the step by step tutorial of the entire pipeline. This is helpful in case you have your own files or just want to run a specific step.

Before running anything make sure you have activated the conda environment using
``conda activate autometa``.

See the :ref:`Install` page for details on setting up your conda environment.

I will be going through this tutorial using the 78Mbp test dataset which can be found here `<https://drive.google.com/drive/u/2/folders/1McxKviIzkPyr8ovj8BG7n_IYk-QfHAgG>`_. You only need to download ``metagenome.fna.gz`` from the above link and save it at a directory as per your liking. I'm saving it in ``$HOME/autometa_run/test_data/``. Now unzip the above metagenome assembly using ``gunzip $HOME/autometa_run/test_data/metagenome.fna.gz``. For instructions on how to download the dataset using command line see the "Using command line" section on :ref:`Benchmarking` page.

I'm also creating an interim and processed directory to separately store the intermediate results and the final results. This is completely optional and is done to properly segregate the final output files.

.. code-block:: bash

    mkdir $HOME/autometa_run/interim
    mkdir $HOME/autometa_run/processed

1. Length filter
----------------

The first step when running autometa is the legth filtering. This would remove any contigs that are below the length cutoff and would generate a new assembly without the filtered contigs. This is useful in removing the noise from the data, as small contigs may have ambigious kmer frequencies. The default cutoff if 3,000bp, ie. any contig that is smaller than 3,000bp would be removed.

.. note::
    It is important that you alter the cutoff based on your N50. If your N50 is really small, e.g. 500bp (pretty common for soil assemblies), then you might want to lower your cutoff to somehwere near N50. The tradeoff with lowing the length cutoff, however, is a greater number of contigs which may make it more difficult for the dataset to be binned.
    
Use the following command to run the length-filter step:

.. code-block:: bash

    autometa-length-filter --assembly $HOME/autometa_run/test_data/78mbp_metagenome.fna \
    --cutoff 3000 --output-fasta $HOME/autometa_run/interim/78mbp_metagenome.filtered.fna \
    --output-stats $HOME/autometa_run/interim/78mbp_metagenome.stats.tsv \
    --output-gc-content $HOME/autometa_run/interim/78mbp_metagenome.gc_content.tsv

Let us dissect the above command:

+---------------------+----------------------------------------------------------------------+
| Flag                |                            Function                                  |
+=====================+======================================================================+
| --assembly          | Path to metagenome assembly (nucleotide fasta) file                  |
+---------------------+----------------------------------------------------------------------+
| --cutoff            | Length cutoff for the filtered assembly                              |
+---------------------+----------------------------------------------------------------------+
| --output-fasta      | Path to filtered metagenomic assembly that would be used for binning |
+---------------------+----------------------------------------------------------------------+
| --output-stats      | Path to assembly statistics table                                    |
+---------------------+----------------------------------------------------------------------+
| --output-gc-content | Path to assembly contigs' GC content and length stats table          |
+---------------------+----------------------------------------------------------------------+

You can view the complete command line opions using ``autometa-length-filter -h``

The above command generates the following files:

+---------------------------------+---------------------------------------------------------------+
| File                            | Description                                                   |
+=================================+===============================================================+
| 78mbp_metagenome.filtered.fna   | Length filtered metagenomic assembly to be used for binning   |
+---------------------------------+---------------------------------------------------------------+
| 78mbp_metagenome.stats.tsv      | Table describing the metagenome assembly statistics           |
+---------------------------------+---------------------------------------------------------------+
| 78mbp_metagenome.gc_content.tsv | Table of GC content and length of each contig in the assembly |
+---------------------------------+---------------------------------------------------------------+


2. Coverage calculation
-----------------------

Coverage calculation for each contig is done to provide another parameter to use while clustering contigs. If you have used SPades to assemble your metagenome, you can use the following command to generate the coverage table:

.. code-block:: bash

    autometa-coverage --assembly $HOME/autometa_run/interim/78mbp_metagenome.fna --from-spades \
    --out $HOME/autometa_run/interim/78mbp_metagenome.coverages.tsv --cpus 40

If you have assembled your metagenome using some other assembler you can use one of the following commands to generate the coverage table.

.. code-block:: bash

    # If you have already made a bed file
    autometa-coverage --assembly $HOME/autometa_run/interim/78mbp_metagenome.filtred.fna \ 
    --bed 78mbp_metagenome.bed --out $HOME/autometa_run/interim/78mbp_metagenome.coverages.tsv --cpus 40

    # If you have already made an alignment (bam file)
    autometa-coverage --assembly $HOME/autometa_run/interim/78mbp_metagenome.filtred.fna \ 
    --bam 78mbp_metagenome.bam --out $HOME/autometa_run/interim/78mbp_metagenome.coverages.tsv \
    --cpus 40

    # If you have already made an alignment (sam file)
    autometa-coverage --assembly $HOME/autometa_run/interim/78mbp_metagenome.filtred.fna \ 
    --sam 78mbp_metagenome.sam --out $HOME/autometa_run/interim/78mbp_metagenome.coverages.tsv \
    --cpus 40

    # If you just have forward and reverse reads
    autometa-coverage --assembly $HOME/autometa_run/interim/78mbp_metagenome.filtred.fna \ 
    --fwd-reads fwd_reads_1.fastq--rev-reads rev_reads_1.fastq \
   --out $HOME/autometa_run/interim/78mbp_metagenome.coverages.tsv --cpus 40

    # In case you have multiple fwd and rev read pairs supply a comma-delimited list (no spaces, fwd and rev lists should be in the same order)
    autometa-coverage --assembly $HOME/autometa_run/interim/78mbp_metagenome.filtred.fna \ 
    --fwd-reads fwd_reads_1.fastq,fwd_reads_2.fastq \ 
    --rev-reads rev_reads_1.fastq,rev_reads_2.fastq \
    --out $HOME/autometa_run/interim/78mbp_metagenome.coverages.tsv --cpus 40

Let us dissect the above commands:

+---------------+----------------------------------------------------------------------------------------------+
| Flag          | Function                                                                                     |
+===============+==============================================================================================+
| --assembly    | Path to length filtered metagenome assembly                                                  |
+---------------+----------------------------------------------------------------------------------------------+
| --from-spades | If the input assembly is generated using SPades then extract k-mer coverages from contig IDs |
+---------------+----------------------------------------------------------------------------------------------+
| --bed         | Path to pre-prepared bed file                                                                |
+---------------+----------------------------------------------------------------------------------------------+
| --bam         | Path to pre-prepared bam file                                                                |
+---------------+----------------------------------------------------------------------------------------------+
| --lengths     | Path to table having length of each contig                                                   |
+---------------+----------------------------------------------------------------------------------------------+
| --sam         | Path to pre-prepared sam file                                                                |
+---------------+----------------------------------------------------------------------------------------------+
| --fwd-reads   | Path to forward reads                                                                        |
+---------------+----------------------------------------------------------------------------------------------+
| --rev-reads   | Path to reverse reads                                                                        |
+---------------+----------------------------------------------------------------------------------------------+
| --cpus        | Number of CPUs to use (default is to use all available CPUs)                                 |
+---------------+----------------------------------------------------------------------------------------------+
| --out         | Path to coverage table of each contig                                                        |
+---------------+----------------------------------------------------------------------------------------------+

You can view the complete command line opions using ``autometa-coverage -h``

The above command would generate the following files:

+--------------------------------+------------------------------------------------------------------+
| File                           | Description                                                      |
+================================+==================================================================+
| 78mbp_metagenome.coverages.tsv | Table with read or k-mer coverage of each contig in the assembly |
+--------------------------------+------------------------------------------------------------------+

3. Generate Open Reading Frames (ORFs)
--------------------------------------

ORF calling using prodigal is performed here. The ORFs are needed for single copy marker gene detection and for taxonomic assignment.

Use the following command to run the ORF calling step:

.. code-block:: bash

    autometa-orfs --assembly $HOME/autometa_run/interim/78mbp_metagenome.filtred.fna \
    --nucls_out $HOME/autometa_run/interim/78mbp_metagenome.orfs.fna --prots_out \
    $HOME/autometa_run/interim/a78mbp_metagenome.orfs.faa --parallel --cpus 90

Let us dissect the above command:

+-------------+--------------------------------------------------------------+
| Flag        | Function                                                     |
+=============+==============================================================+
| --assembly  | Path to length filtered metagenome assembly                  |
+-------------+--------------------------------------------------------------+
| --nucls_out | Path to nucleic acid sequence of ORFs                        |
+-------------+--------------------------------------------------------------+
| --prots_out | Path to amino acid sequence of ORFs                          |
+-------------+--------------------------------------------------------------+
| --parallel  | Enable GNU parallel (deafult is False)                       |
+-------------+--------------------------------------------------------------+
| --cpus      | Number of CPUs to use (default is to use all available CPUs) |
+-------------+--------------------------------------------------------------+

You can view the complete command line opions using ``autometa-orfs -h``

The above command would generate the following files:

+---------------------------+-------------------------------+
| File                      | Description                   |
+===========================+===============================+
| 78mbp_metagenome.orfs.fna | Nucleic acid sequence of ORFs |
+---------------------------+-------------------------------+
| 78mbp_metagenome.orfs.faa | Amino acid sequence of ORFs   |
+---------------------------+-------------------------------+

4. Single copy markers
----------------------

Autometa uses single-copy markers to guide clustering, and does not assume that recoverable genomes will necessarily be ‘complete’. You first need to download the single-copy markers.

.. code-block:: bash

    #Create a markers directory to hold the marker genes
    mkdir /Autometa/autometa/databases/markers
    # Change the default download path to the directory created above
    autometa-config --section databases --option markers --value /Autometa/autometa/databases/markers
    # Download single-copy marker genes
    autometa-update-databases --update-markers
    # hmmpress the marker genes
    hmmpress -f /Autometa/autometa/databases/markers/bacteria.single_copy.hmm
    hmmpress -f /Autometa/autometa/databases/markers/archaea.single_copy.hmm

Use the following command to annotate contigs containing single copy marker genes:

.. code-block:: bash

    autometa-markers --orfs $HOME/autometa_run/interim/78mbp_metagenome.orfs.faa \
    --kingdom bacteria --hmmscan $HOME/autometa_run/interim/78mbp_metagenome.hmmscan.tsv \
    --parallel --cpus 90 --seed 42 --out $HOME/autometa_run/interim/78mbp_metagenome.markers.tsv

Let us disect the above command:

+------------+-----------------------------------------------------------------------------------------------+
| Flag       | Function                                                                                      |
+============+===============================================================================================+
| --orfs     | Path to fasta file containing amino acid sequences of ORFS                                    |
+------------+-----------------------------------------------------------------------------------------------+
| --kingdom  | Kingdom to search for markers (default: bacteria). Choices bacteria and archaea               |
+------------+-----------------------------------------------------------------------------------------------+
| --hmmscan  | Path to hmmscan output table containing the respective kingdom single-copy marker annotations |
+------------+-----------------------------------------------------------------------------------------------+
| --parallel | Use hmmscan parallel option (default: False)                                                  |
+------------+-----------------------------------------------------------------------------------------------+
| --cpus     | Number of CPUs to use (default is to use all available CPUs)                                  |
+------------+-----------------------------------------------------------------------------------------------+
| --seed     | Seed to set random state for hmmscan. (default: 42)                                           |
+------------+-----------------------------------------------------------------------------------------------+
| --out      | Path to write filtered annotated markers corresponding to kingdom                             |
+------------+-----------------------------------------------------------------------------------------------+

You can view the complete command line opions using ``autometa-markers -h``

The above command would generate the following files:

+------------------------------+---------------------------------------------------------------------------------------+
| File                         | Description                                                                           |
+==============================+=======================================================================================+
| 78mbp_metagenome.hmmscan.tsv | hmmscan output table containing the respective kingdom single-copy marker annotations |
+------------------------------+---------------------------------------------------------------------------------------+
| 78mbp_metagenome.markers.tsv | Annotated marker table corresponding to the particular kingdom                        |
+------------------------------+---------------------------------------------------------------------------------------+

5. Taxonomy assignment: BLASTP
------------------------------

Autometa assigns a taxonomic rank to each contig and then takes only the contig belonging to the specified kingdom (either bacteria or archaea) for binning. We found that in host-associated metagenomes, this step vastly improves the binning performance of Autometa (and other pipelines) because less eukaryotic or viral contigs will be binned into bacterial bins. 

The first step for contig taxonomy assignment is a local alignment search of the ORFs against a reference database. This can be accelerated using `diamond <https://github.com/bbuchfink/diamond>`_.

Create a diamond formatted database of the NCBI non-redundant (nr) protein database.

.. code-block:: bash

    diamond makedb --in /Autometa/autometa/databases/ncbi/nr --db /Autometa/autometa/databases/ncbi/nr -p 40

Breaking down the above command:

+------+--------------------------------------+
| Flag | Function                             |
+======+======================================+
| --in | Path to nr database                  |
+------+--------------------------------------+
| --db | Path to diamond formated nr database |
+------+--------------------------------------+
| -p   | Number of processors to use          |
+------+--------------------------------------+

Run diamond blastp using the following command:

.. code-block:: bash

    diamond blastp --query $HOME/autometa_run/interim/78mbp_metagenome.orfs.faa \
    --db /Autometa/autometa/databases/ncbi/nr.dmnd --evalue 1e-5 \
    --max-target-seqs 200 --threads 40 --outfmt 6 \
    --out $HOME/autometa_run/interim/78mbp_metagenome.blastp.tsv

Breaking down the above command:

+-------------------+-----------------------------------------------------------------------+
| Flag              | Function                                                              |
+===================+=======================================================================+
| --query           | Path to query sequence. Here, amino acid sequence of ORFs             |
+-------------------+-----------------------------------------------------------------------+
| --db              | Path to diamond formatted nr database                                 |
+-------------------+-----------------------------------------------------------------------+
| --evalue          | Maximum expected value to report an alignment                         |
+-------------------+-----------------------------------------------------------------------+
| --max-target-seqs | Maximum number of target sequences per query to report alignments for |
+-------------------+-----------------------------------------------------------------------+
| --threads         | Number of processors to use                                           |
+-------------------+-----------------------------------------------------------------------+
| --outfmt          | Output format of BLASTP results                                       |
+-------------------+-----------------------------------------------------------------------+
| --out             | Path to BLASTP results                                                |
+-------------------+-----------------------------------------------------------------------+

To see the complete list of acceptable output formats see Diamond `GitHub Wiki <https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#output-options>`__. A complete list of all command line options for Diamond can be found on its `GitHub Wiki <https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options>`__.

.. note:: 
    Autometa only parses output format 6 provided above as: ``--outfmt 6``


The above command would generate the following files:

#. If --split-rank-and-write is specified then it will split contigs by provided canonical-rank column then write a file corresponding that rank. Eg. Bacteria.fasta, Archaea.fasta, etc for superkingdom.

6. Taxonomy assignment: LCA
---------------------------

The second step in taxon assignment is finding out the lowest common ancestor (LCA). This step uses the blastp results generated in the previous step to generate a table having the LCA of each ORF.

Use the following command to run the LCA:

.. code-block:: bash

    autometa-taxonomy-lca --blast $HOME/autometa_run/interim/78mbp_metagenome.blastp.tsv \
    --dbdir /Autometa/autometa/databases/ncbi/ \
    --output $HOME/autometa_run/interim/78mbp_metagenome.lca.tsv

Let us dissect the above command:

+----------+-----------------------------------------+
| Flag     | Function                                |
+==========+=========================================+
| --blast  | Path to diamond balstp output           |
+----------+-----------------------------------------+
| --dbdir  | Path to directory having ncbi databases |
+----------+-----------------------------------------+
| --output | Path to write LCA results               |
+----------+-----------------------------------------+

You can view the complete command line opions using ``autometa-taxonomy-lca -h``

The above command would generate a table (``78mbp_metagenome.lca.tsv``) having the name, rank and taxid of the LCA for each ORF.

7. Taxonomy assignment: Majority vote
-------------------------------------

The next step in taxone assignment is doing a majority vote to decide the taxonomy of each contig. A vote system helps in minimizing the effect of horizontal gene transfer (HGT) as even if some ORFs on the contig are divergent there will be other that belong to the organism, thus preventing a complete misclassification of HGT contigs.

You can run the majority vote step using the following command:

.. code-block:: bash

    autometa-taxonomy-majority-vote --lca $HOME/autometa_run/interim/78mbp_metagenome.lca.tsv \
    --output $HOME/autometa_run/interim/78mbp_metagenome.votes.tsv \
    --dbdir /Autometa/autometa/databases/ncbi/

Let us dissect the above command:

+----------+-----------------------------------+
| Flag     | Function                          |
+==========+===================================+
| --lca    | Path to LCA table                 |
+----------+-----------------------------------+
| --output | Path to write majority vote table |
+----------+-----------------------------------+
| --dbdir  | Path to ncbi database directory   |
+----------+-----------------------------------+

You can view the complete command line opions using ``autometa-taxonomy-majority-vote -h``

The above command would generate a table (``78mbp_metagenome.votes.tsv``) having the taxid of each contig identified as per majority vote.

8. Taxonomy assignment: Split kingdoms
--------------------------------------

In this final step of taxon assignment we use the voted taxid of each contig to split the contigs in different kingdoms and write them as per the provided canonical rank.

.. code-block:: bash

    autometa-taxonomy --input $HOME/autometa_run/interim/78mbp_metagenome.votes.tsv \
    --output $HOME/autometa_run/interim/ \
    --assembly $HOME/autometa_run/interim/78mbp_metagenome.filtered.fna \
    --prefix 78mbp_metagenome --split-rank-and-write superkingdom \
    --ncbi /Autometa/autometa/databases/ncbi/

Let us dissect the above command:

+------------------------+--------------------------------------------------------------------------------+
| Flag                   | Function                                                                       |
+========================+================================================================================+
| --input                | Path to voted taxids table                                                     |
+------------------------+--------------------------------------------------------------------------------+
| --output               | Directory to output fasta files of split canonical ranks and taxonomy.tsv      |
+------------------------+--------------------------------------------------------------------------------+
| --assembly             | Path to filtered metagenome assembly                                           |
+------------------------+--------------------------------------------------------------------------------+
| --prefix               | prefix to use for each file written                                            |
+------------------------+--------------------------------------------------------------------------------+
| --split-rank-and-write | Split contigs by provided canonical-rank column then write to output directory |
+------------------------+--------------------------------------------------------------------------------+
| --ncbi                 | Path to ncbi database directory                                                |
+------------------------+--------------------------------------------------------------------------------+

Other options available for ``--split-rank-and-write`` are phylum, class, order, family, genus and species

You can view the complete command line opions using ``autometa-taxonomy -h``

+-----------------------------------+------------------------------------------------------------------------------------------+
| File                              | Description                                                                              |
+===================================+==========================================================================================+
| 78mbp_metagenome.taxonomy.tsv     | Table with taxonomic classification of each contig                                       |
+-----------------------------------+------------------------------------------------------------------------------------------+
| 78mbp_metagenome.bacteria.fna     | Fasta file having the nucleic acid sequence of all bacterial contigs                     |
+-----------------------------------+------------------------------------------------------------------------------------------+
| 78mbp_metagenome.unclassified.fna | Fasta file having the nucleic acid sequence of all contigs unclassified at kingdom level |
+-----------------------------------+------------------------------------------------------------------------------------------+

In my case there are no non-bacterial contigs. For your dataset, ``autometa-taxonomy`` will produce other fasta files, for example Eukaryota.fasta and Viruses.fasta.

9. K-mer counting
-----------------

A k-mer (`ref <https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/>`_) is just a sequence of k characters in a string (or nucleotides in a DNA sequence). It is known that contigs that belong to the same genome have similar k-mer composition. Here, we compute k-mer frequencies of only the bacterial contigs.

This step does the following:

#. Create a  k-mer matrix of k^4/2 dimensions using the specified k-mer frequency (default is k-mer of size 5 bp). k-mer size can be altered using the ``--size`` flag
#. Normaization of the k-mer matrix (default embedding method is am_clr). Normalization method can be altered using ``--norm-method`` flag
#. Reduce the dimensions of k-mer frequencies using principle component analysis (PCA). Default PCA dimensions are 50. This can be altered using the ``--pca-dimensions`` flag 
#. Embedding the PCA dimensions into two dimensions (default embedding method is BH-tSNE) to allow the ease of visualization and manual binning of the contigs (see `ViZBin <https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-014-0066-1>`_ paper). . Embedding method can be altered using ``--embedding-method`` flag

Use the following command to run the k-mer counting step:

.. code-block:: bash

    autometa-kmers --fasta $HOME/autometa_run/interim/78mbp_metagenome.bacteria.fna \
    --kmers $HOME/autometa_run/interim/78mbp_metagenome.bacteria.kmers.tsv --size 5 \
    --norm-output $HOME/autometa_run/interim/78mbp_metagenome.bacteria.kmers.normalized.tsv \
    --norm-method am_clr --pca-dimensions 50 \
    --embedding-output $HOME/autometa_run/processed/78mbp_metagenome.bacteria.kmers.embedded.tsv \
    --embedding-method bhsne --cpus 40 --seed 42

If you noticed I stored the ``78mbp_metagenome.bacteria.kmers.embedded.tsv`` file in the ``processed`` directory as no further analysis is required on the final.

.. note::
    In case you put ``--pca-dimensions`` as zero then autometa will skip PCA.

Let us dissect the above command:

+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| Flag               | Function                                                                                                                 |
+====================+==========================================================================================================================+
| --fasta            | Path to the fasta file having only bacterial contigs                                                                     |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| --kmers            | Path to k-mer frequency table                                                                                            |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| --size             | k-mer size in bp (default 5bp)                                                                                           |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| --norm-output      | Path to normalized k-mer table                                                                                           |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| --norm-method      | Normalization method to transform kmer counts prior to PCA and embedding (default am_clr). Choices : ilr, clr and am_clr |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| --pca-dimensions   | Number of dimensions to reduce to PCA feature space after normalization and prior to embedding (default: 50)             |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| --embedding-output | Path to embedded k-mer table                                                                                             |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| --embedding-method | Embedding method to reduce the k-mer frequencies. Choices: shsne, bhsne, umap.                                           |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| --cpus             | Number of CPUs to use (default is to use all available CPUs)                                                             |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+
| --seed             | Set random seed for dimension reduction determinism (default 42). Useful in replicating the results                      |
+--------------------+--------------------------------------------------------------------------------------------------------------------------+

You can view the complete command line options using ``autometa-kmers -h``

In the above command k-mer normalization is being done using Autometa's implementation of center
log-ratio transform (am_clr). Other available normalization methods are isometric log-ratio transform (ilr, scikit-bio implementation) and center log-ratio transform (clr, scikit-bio implementation)

In the above command k-mer embedding is being done using Barnes-Hut Stochastic Neighbor Embedding (BH-tSNE). Other embedding methods that are available are Uniform Manifold Approximation and Projection (UMAP) and SKSNE (BH-tSNE is the default). bhsne and sksne are two different implementations of BH-tSNE from tsne and scikit-learn respectively, that appear to give very different results. We recommend using the former.

The above command generates the following files:

+---------------------------------------+--------------------------------------------------------+
| File                                  | Description                                            |
+=======================================+========================================================+
| 78mbp_metagenome.kmers.tsv            | Table with raw k-mer frequencies of each contig        |
+---------------------------------------+--------------------------------------------------------+
| 78mbp_metagenome.kmers.normalized.tsv | Table with normalized k-mer frequencies of each contig |
+---------------------------------------+--------------------------------------------------------+
| 78mbp_metagenome.kmers.embedded.tsv   | Table with embedded k-mer frequencies of each contig   |
+---------------------------------------+--------------------------------------------------------+

10. Binning
-----------

This is the step where contigs are binned into genomes. Autometa assesses clusters by examining both their completeness (number of expected single copy markers) and purity (number of single copy markers that are unique in the cluster). If we supply a taxonomy table, then that is also used to help with clustering. Otherwise, Autometa clusters solely on 5-mer frequency and coverage. 

This step does the following:

#. Find single-copy marker genes in the input contigs with HMMER
#. Bin contigs based on embedded k-mer coordinates, coverage and (optionally) taxonomy
#. Accept clusters that are estimated to be over 20% complete and 95% pure based on single-copy marker genes. These are default papameteres and can be altered to suit your needs. Completeness can be altered using the ``--completeness`` flag and purity using the ``--purity`` flag
#. Unclustered contigs leftover will be re-clustered until no more acceptable clusters are yielded

If you include a taxonomy table Autometa will attempt to further partition the data based on ascending taxonomic specificity (i.e. in the order superkingdom, phylum, class, order, family, genus, species) when binning unbinned contigs from a previous attempt. We found that this is mainly useful if you have a highly complex metagenome (lots of species), or you have several related species at similar coverage level.

Use the following command to run the binning:

.. code-block:: bash

    autometa-binning --kmers $HOME/autometa_run/interim/78mbp_metagenome.bacteria.kmers.normalized.tsv \
    --coverages $HOME/autometa_run/interim/78mbp_metagenome.coverages.tsv \
    --gc-content $HOME/autometa_run/interim/78mbp_metagenome.gc_content.tsv \
    --markers $HOME/autometa_run/interim/78mbp_metagenome.markers.tsv \
    --embedded-kmers $HOME/autometa_run/interim/78mbp_metagenome.bacteria.kmers.embedded.tsv \
    --clustering-method dbscan --completeness 20 --purity 90 --cov-stddev-limit 25 \
    --gc-stddev-limit 5 --taxonomy $HOME/autometa_run/interim/78mbp_metagenome.taxonomy.tsv \
    --output-binning $HOME/autometa_run/processed/78mbp_metagenome.binning.tsv \
    --output-main $HOME/autometa_run/processed/78mbp_metagenome.main.tsv \
    --starting-rank superkingdom --domain bacteria 

Since these are the final binning results we store them in the ``processed`` directory.

Let us dissect the above command:

+---------------------+-----------------------------------------------------------------------------------------+
| Flag                | Function                                                                                |
+=====================+=========================================================================================+
| --kmers             | Path to normalized k-mer frequencies table                                              |
+---------------------+-----------------------------------------------------------------------------------------+
| --coverages         | Path to metagenome coverages table                                                      |
+---------------------+-----------------------------------------------------------------------------------------+
| --gc-content        | Path to metagenome GC contents table                                                    |
+---------------------+-----------------------------------------------------------------------------------------+
| --markers           | Path to Autometa annotated markers table                                                |
+---------------------+-----------------------------------------------------------------------------------------+
| --embedded-kmers    | Path to provide embedded k-mer frequencies table                                        |
+---------------------+-----------------------------------------------------------------------------------------+
| --clustering-method | Clustering algorithm to use for recursive binning. Choices dbscan (default) and hdbscan |
+---------------------+-----------------------------------------------------------------------------------------+
| --completeness      | completeness cutoff to retain cluster (default 20)                                      |
+---------------------+-----------------------------------------------------------------------------------------+
| --purity            | purity cutoff to retain cluster (default 95)                                            |
+---------------------+-----------------------------------------------------------------------------------------+
| --cov-stddev-limit  | coverage standard deviation limit to retain cluster (default 25)                        |
+---------------------+-----------------------------------------------------------------------------------------+
| --gc-stddev-limit   | GC content standard deviation limit to retain cluster (default 5)                       |
+---------------------+-----------------------------------------------------------------------------------------+
| --taxonomy          | Path to Autometa assigned taxonomies table                                              |
+---------------------+-----------------------------------------------------------------------------------------+
| --output-binning    | Path to write Autometa binning results                                                  |
+---------------------+-----------------------------------------------------------------------------------------+
| --output-main       | Path to write Autometa main table                                                       |
+---------------------+-----------------------------------------------------------------------------------------+
| --starting-rank     | Canonical rank at which to begin subsetting taxonomy (default: superkingdom)            |
+---------------------+-----------------------------------------------------------------------------------------+
| --domain            | Kingdom to consider. Choices bacteria (default) and archaea                             |
+---------------------+-----------------------------------------------------------------------------------------+

There are two binning algorithms to chose from Density-Based Spatial Clustering of Applications with Noise (DBSCAN) and Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN). The default is DBSCAN.

You can view the complete command line options using ``autometa-binning -h``

The above command generates the following files:

#. ``78mbp_metagenome.binning.tsv`` contains the final binning results along with a few more metrics regarding each genome bin.
#. ``78mbp_metagenome.main.tsv`` which contains the feature table that was utilized during the genome binning process as well as the corresponding output predictions.

The following table describes each column for the resulting binning outputs. We'll start with the columns present in ``78mbp_metagenome.binning.tsv`` then describe the additional columns that are present in ``78mbp_metagenome.main.tsv``.

+-------------------+---------------------------------------------------------------------------------------------------------------------+
| Column            | Description                                                                                                         |
+===================+=====================================================================================================================+
| Contig            | Name of the contig in the input fasta file                                                                          |
+-------------------+---------------------------------------------------------------------------------------------------------------------+
| Cluster           | Cluster assigned by autometa to the contig                                                                          |
+-------------------+---------------------------------------------------------------------------------------------------------------------+
| Completeness      | Estimated completeness of the cluster, based on single-copy marker genes                                            |
+-------------------+---------------------------------------------------------------------------------------------------------------------+
| Purity            | Estimated purity of the cluster, based on the number of single-copy marker genes that are duplicated in the cluster |
+-------------------+---------------------------------------------------------------------------------------------------------------------+
| coverage_stddev   | Coverage standard deviation of the cluster                                                                          |
+-------------------+---------------------------------------------------------------------------------------------------------------------+
| gc_content_stddev | GC content standard deviation of the cluster                                                                        |
+-------------------+---------------------------------------------------------------------------------------------------------------------+

Description of additional columns in ``78mbp_metagenome.main.tsv``:

+--------------+-------------------------------------------------+
| Column       | Description                                     |
+==============+=================================================+
| Coverage     | Estimated coverage of the contig                |
+--------------+-------------------------------------------------+
| gc_content   | Estimated GC content of the contig              |
+--------------+-------------------------------------------------+
| Length       | Estimated length of the contig                  |
+--------------+-------------------------------------------------+
| Species      | Assigned taxonomic species for the contig       |
+--------------+-------------------------------------------------+
| Genus        | Assigned taxonomic genus for the contig         |
+--------------+-------------------------------------------------+
| Family       | Assigned taxonomic family for the contig        |
+--------------+-------------------------------------------------+
| Order        | Assigned taxonomic order for the contig         |
+--------------+-------------------------------------------------+
| Class        | Assigned taxonomic class for the contig         |
+--------------+-------------------------------------------------+
| Phylum       | Assigned taxonomic phylum for the contig        |
+--------------+-------------------------------------------------+
| Superkingdom | Assigned taxonomic superkingdom for the contig  |
+--------------+-------------------------------------------------+
| taxid        | Assigned NCBI taxonomy ID number for the contig |
+--------------+-------------------------------------------------+
| x_1          | The first coordinate after dimension reduction  |
+--------------+-------------------------------------------------+
| x_2          | The second coordinate after dimension reduction |
+--------------+-------------------------------------------------+

You can attempt to improve your genome bins with an unclustered recruitment step which uses features from existing genome bins to recruit unbinned contigs. Alternatively you can use these initial genome bin predictions and continue to the :ref:`Examining Results` section.

11. Unclustered recruitment (Optional)
--------------------------------------

Supervised machine learning is used to classify the unclustered contigs to the bins that we have produced. This steop is optional and the results should be verified (see Note below) before going ahead with it.

.. note::
    The machine learning step has been observed to bin contigs that do not necessary belong to the predicted genome. Careful inspection of coverage and taxonomy should be done before proceed to use these results.

Use the following command to run the unclustered recruitment step:

.. code-block:: bash

    autometa-unclustered-recruitment \
    --kmers $HOME/autometa_run/interim/78mbp_metagenome.bacteria.kmers.normalized.tsv \
    --coverage $HOME/autometa_run/interim/78mbp_metagenome.coverages.tsv \
    --binning $HOME/autometa_run/interim/78mbp_metagenome.binning.tsv \
    --markers $HOME/autometa_run/interim/78mbp_metagenome.markers.tsv \
    --taxonomy $HOME/autometa_run/interim/78mbp_metagenome.taxonomy.tsv \
    --output-binning $HOME/autometa_run/processed/78mbp_metagenome.recruitment.tsv \
    --output-main $HOME/autometa_run/processed/78mbp_metagenome.recruitment.main.tsv \
    --classifier decision_tree --seed 42

Since these are the final binning results we store them in the ``processed`` directory.

Let us dissect the above command:

+------------------+-------------------------------------------------------------------------------------------------+
| Flag             | Function                                                                                        |
+==================+=================================================================================================+
| --kmers          | Path to normalized k-mer frequencies table                                                      |
+------------------+-------------------------------------------------------------------------------------------------+
| --coverages      | Path to metagenome coverages table                                                              |
+------------------+-------------------------------------------------------------------------------------------------+
| --binning        | Path to genome bin assignments                                                                  |
+------------------+-------------------------------------------------------------------------------------------------+
| --markers        | Path to Autometa annotated markers table                                                        |
+------------------+-------------------------------------------------------------------------------------------------+
| --taxonomy       | Path to Autometa assigned taxonomies table                                                      |
+------------------+-------------------------------------------------------------------------------------------------+
| --output-binning | Path to output unclustered recruitment table                                                    |
+------------------+-------------------------------------------------------------------------------------------------+
| --output-main    | Path to write Autometa main table                                                               |
+------------------+-------------------------------------------------------------------------------------------------+
| --classifier     | classifier to use for recruitment of contigs. Choices decision_tree (default) and random_forest |
+------------------+-------------------------------------------------------------------------------------------------+
| --seed           | Seed to use for RandomState when initializing classifiers                                       |
+------------------+-------------------------------------------------------------------------------------------------+

You can view the complete command line options using ``autometa-unclustered-recruitment -h``

The above command would generate ``78mbp_metagenome.recruitment.tsv`` and ``78mbp_metagenome.recruitment.main.tsv``.

``78mbp_metagenome.recruitment.tsv`` contains the final predictions of ``autometa-unclustered-recruitment``. ``78mbp_metagenome.recruitment.main.tsv`` is the feature table with corresponding predictions utilized during/after the unclustered recruitment algorithm. This represents unbinned contigs with their respective annotations and output predictions of their recruitment into a genome bin. The taxonomic features have been encoded using “one-hot encoding” or a presence/absence matrix where each column is a canonical taxonomic rank and its respective value for each row represents its presence or absence. Presence and absence are denoted with 1 and 0, respectively. Hence ‘one-hot’ encoding being an encoding of presence and absence of the respective annotation type. In our case taxonomic designation.

.. todo::
    Add file description of ``78mbp_metagenome.recruitment.tsv`` after evan edits it.

Running modules
---------------

Many of the Autometa modules may be run standalone.

Simply pass in the ``-m`` flag when calling a script to signify to python you are
running an Autometa *module*.

I.e. ``python -m autometa.common.kmers -h``

Running functions
-----------------

Many of the Autometa functions may be run standalone as well. It is same as importing any other python
function.

.. code-block:: python

    from autometa.common.external import samtools

    samtools.sort(sam=<path/to/sam/file>, out=<path/to/output/file>, nproc=4)

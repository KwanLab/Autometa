.. _step-by-step-tutorial:

=================================
📓 Bash Step by Step Tutorial 📓
=================================

Here is the step by step tutorial of on running the entire pipeline manually through Bash. This is helpful in case you have your own files or just want to run a specific step.
If you would like to set up a run of the whole pipeline through Bash, see the :ref:`Bash Workflow<🐚 Bash Workflow 🐚>` section.

Before running anything make sure you have activated the conda environment using
``conda activate autometa``.

See the :ref:`Autometa Package Installation` page for details on setting up your conda environment.

I will be going through this tutorial using the 78Mbp test dataset which can be found here `<https://drive.google.com/drive/u/2/folders/1McxKviIzkPyr8ovj8BG7n_IYk-QfHAgG>`_.
You only need to download ``metagenome.fna.gz`` from the above link and save it at a directory as per your liking. I'm saving it in ``$HOME/tutorial/test_data/``.
For instructions on how to download the dataset using command-line see the "Using command-line" section on :ref:`Benchmarking` page.

1. Length filter
----------------

The first step when running Autometa is the length filtering. This would remove any contigs that are below the length cutoff. This is useful in removing the noise from the data,
as small contigs may have ambiguous kmer frequencies. The default cutoff if 3,000bp, ie. any contig that is smaller than 3,000bp would be removed.

.. note::
    It is important that you alter the cutoff based on your N50. If your N50 is really small, e.g. 500bp (pretty common for soil assemblies),
    then you might want to lower your cutoff to somewhere near N50. The tradeoff with lowering the length cutoff, however, is a greater number of
    contigs which may make it more difficult for the dataset to be binned. As was shown in the `Autometa <https://academic.oup.com/nar/article/47/10/e57/5369936>`_ paper,
    as assembly quality degrades so does the binning performance.

Use the following command to run the length-filter step:

.. code-block:: bash

    autometa-length-filter \
        --assembly $HOME/tutorial/test_data/78mbp_metagenome.fna \
        --cutoff 3000 \
        --output-fasta $HOME/tutorial/78mbp_metagenome.filtered.fna \
        --output-stats $HOME/tutorial/78mbp_metagenome.stats.tsv \
        --output-gc-content $HOME/tutorial/78mbp_metagenome.gc_content.tsv

Let us dissect the above command:

+-------------------------+----------------------------------------------------------------------+-------------+
| Flag                    |                            Input arguments                           | Requirement |
+=========================+======================================================================+=============+
| ``--assembly``          | Path to metagenome assembly (nucleotide fasta) file                  | Required    |
+-------------------------+----------------------------------------------------------------------+-------------+
| ``--cutoff``            | Length cutoff for the filtered assembly. Default is 3,000bp          | Optional    |
+-------------------------+----------------------------------------------------------------------+-------------+
| ``--output-fasta``      | Path to filtered metagenomic assembly that would be used for binning | Required    |
+-------------------------+----------------------------------------------------------------------+-------------+
| ``--output-stats``      | Path to assembly statistics table                                    | Optional    |
+-------------------------+----------------------------------------------------------------------+-------------+
| ``--output-gc-content`` | Path to assembly contigs' GC content and length stats table          | Optional    |
+-------------------------+----------------------------------------------------------------------+-------------+

You can view the complete command-line options using ``autometa-length-filter -h``

The above command generates the following files:

+---------------------------------+------------------------------------------------------------------------+
| File                            | Description                                                            |
+=================================+========================================================================+
| 78mbp_metagenome.filtered.fna   | Length filtered metagenomic assembly to be used for binning            |
+---------------------------------+------------------------------------------------------------------------+
| 78mbp_metagenome.stats.tsv      | Table describing the filtered metagenome assembly statistics           |
+---------------------------------+------------------------------------------------------------------------+
| 78mbp_metagenome.gc_content.tsv | Table of GC content and length of each contig in the filtered assembly |
+---------------------------------+------------------------------------------------------------------------+

.. _coverage-calculation:

2. Coverage calculation
-----------------------

Coverage calculation for each contig is done to provide another parameter to use while clustering contigs.

from SPAdes
^^^^^^^^^^^

If you have used SPAdes to assemble your metagenome, you can use the following command to generate the coverage table:

.. code-block:: bash

    autometa-coverage \
        --assembly $HOME/tutorial/78mbp_metagenome.fna \
        --out $HOME/tutorial/78mbp_metagenome.coverages.tsv \
        --from-spades

from alignments.bed
^^^^^^^^^^^^^^^^^^^

If you have assembled your metagenome using some other assembler you can use one of the following commands to generate the coverage table.

.. code-block:: bash

    # If you have already made a bed file
    autometa-coverage \
        --assembly $HOME/tutorial/78mbp_metagenome.filtered.fna \
        --bed 78mbp_metagenome.bed \
        --out $HOME/tutorial/78mbp_metagenome.coverages.tsv \
        --cpus 40

from alignments.bam
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # If you have already made an alignment (bam file)
    autometa-coverage \
        --assembly $HOME/tutorial/78mbp_metagenome.filtered.fna \
        --bam 78mbp_metagenome.bam \
        --out $HOME/tutorial/78mbp_metagenome.coverages.tsv \
        --cpus 40

from alignments.sam
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # If you have already made an alignment (sam file)
    autometa-coverage \
        --assembly $HOME/tutorial/78mbp_metagenome.filtered.fna \
        --sam 78mbp_metagenome.sam \
        --out $HOME/tutorial/78mbp_metagenome.coverages.tsv \
        --cpus 40

from paired-end reads
^^^^^^^^^^^^^^^^^^^^^

You may calculate coverage using forward and reverse reads with the assembled metagenome.

.. code-block:: bash

    autometa-coverage \
        --assembly $HOME/tutorial/78mbp_metagenome.filtered.fna \
        --fwd-reads fwd_reads_1.fastq \
        --rev-reads rev_reads_1.fastq \
        --out $HOME/tutorial/78mbp_metagenome.coverages.tsv \
        --cpus 40

In case you have multiple forward and reverse read pairs supply a comma-delimited list.

.. code-block:: bash

    autometa-coverage \
        --assembly $HOME/tutorial/78mbp_metagenome.filtered.fna \
        --fwd-reads fwd_reads_1.fastq,fwd_reads_2.fastq \
        --rev-reads rev_reads_1.fastq,rev_reads_2.fastq \
        --out $HOME/tutorial/78mbp_metagenome.coverages.tsv \
        --cpus 40

.. note::

    1. No spaces should be used when providing the forward and reverse reads.
    2. The lists of forward and reverse reads should be in the order corresponding to their respective reads pair.

Let us dissect the above commands:

+-------------------+----------------------------------------------------------------------------------------------+
| Flag              | Function                                                                                     |
+===================+==============================================================================================+
| ``--assembly``    | Path to length filtered metagenome assembly                                                  |
+-------------------+----------------------------------------------------------------------------------------------+
| ``--from-spades`` | If the input assembly is generated using SPades then extract k-mer coverages from contig IDs |
+-------------------+----------------------------------------------------------------------------------------------+
| ``--bed``         | Path to alignments BED file                                                                  |
+-------------------+----------------------------------------------------------------------------------------------+
| ``--bed``         | Path to alignments BAM file                                                                  |
+-------------------+----------------------------------------------------------------------------------------------+
| ``--sam``         | Path to alignments SAM file                                                                  |
+-------------------+----------------------------------------------------------------------------------------------+
| ``--fwd-reads``   | Path to forward reads                                                                        |
+-------------------+----------------------------------------------------------------------------------------------+
| ``--rev-reads``   | Path to reverse reads                                                                        |
+-------------------+----------------------------------------------------------------------------------------------+
| ``--cpus``        | Number of CPUs to use (default is to use all available CPUs)                                 |
+-------------------+----------------------------------------------------------------------------------------------+
| ``--out``         | Path to coverage table of each contig                                                        |
+-------------------+----------------------------------------------------------------------------------------------+

You can view the complete command-line options using ``autometa-coverage -h``

The above command would generate the following files:

+--------------------------------+--------------------------------------------------------------------+
| File                           | Description                                                        |
+================================+====================================================================+
| 78mbp_metagenome.coverages.tsv | Table with read or k-mer coverage of each contig in the metagenome |
+--------------------------------+--------------------------------------------------------------------+

3. Generate Open Reading Frames (ORFs)
--------------------------------------

ORF calling using prodigal is performed here. The ORFs are needed for single copy marker gene detection and for taxonomic assignment.

Use the following command to run the ORF calling step:

.. code-block:: bash

    autometa-orfs \
        --assembly $HOME/tutorial/78mbp_metagenome.filtered.fna \
        --output-nucls $HOME/tutorial/78mbp_metagenome.orfs.fna \
        --output-prots $HOME/tutorial/a78mbp_metagenome.orfs.faa \
        --cpus 40

Let us dissect the above command:

+--------------------+--------------------------------------------------------------+
| Flag               | Function                                                     |
+====================+==============================================================+
| ``--assembly``     | Path to length filtered metagenome assembly                  |
+--------------------+--------------------------------------------------------------+
| ``--output-nucls`` | Path to nucleic acid sequence of ORFs                        |
+--------------------+--------------------------------------------------------------+
| ``--output-prots`` | Path to amino acid sequence of ORFs                          |
+--------------------+--------------------------------------------------------------+
| ``--cpus``         | Number of CPUs to use (default is to use all available CPUs) |
+--------------------+--------------------------------------------------------------+

You can view the complete command-line options using ``autometa-orfs -h``

The above command would generate the following files:

+---------------------------+---------------------------------+
| File                      | Description                     |
+===========================+=================================+
| 78mbp_metagenome.orfs.fna | Nucleic acid fasta file of ORFs |
+---------------------------+---------------------------------+
| 78mbp_metagenome.orfs.faa | Amino acid fasta file of ORFs   |
+---------------------------+---------------------------------+

4. Single copy markers
----------------------

Autometa uses single-copy markers to guide clustering, and does not assume that recoverable genomes will necessarily be "complete". You first need to download the single-copy markers.

.. code-block:: bash

    # Create a markers directory to hold the marker genes
    mkdir -p $HOME/Autometa/autometa/databases/markers

    # Change the default download path to the directory created above
    autometa-config \
        --section databases \
        --option markers \
        --value $HOME/Autometa/autometa/databases/markers

    # Download single-copy marker genes
    autometa-update-databases --update-markers

    # hmmpress the marker genes
    hmmpress -f $HOME/Autometa/autometa/databases/markers/bacteria.single_copy.hmm
    hmmpress -f $HOME/Autometa/autometa/databases/markers/archaea.single_copy.hmm

Use the following command to annotate contigs containing single-copy marker genes:

.. code-block:: bash

    autometa-markers \
        --orfs $HOME/tutorial/78mbp_metagenome.orfs.faa \
        --kingdom bacteria \
        --hmmscan $HOME/tutorial/78mbp_metagenome.hmmscan.tsv \
        --out $HOME/tutorial/78mbp_metagenome.markers.tsv \
        --parallel \
        --cpus 4 \
        --seed 42

Let us dissect the above command:

+----------------+-----------------------------------------------------------------------------------------------+-------------+
| Flag           | Function                                                                                      | Requirement |
+================+===============================================================================================+=============+
| ``--orfs``     | Path to fasta file containing amino acid sequences of ORFS                                    | Required    |
+----------------+-----------------------------------------------------------------------------------------------+-------------+
| ``--kingdom``  | Kingdom to search for markers. Choices bacteria (default) and archaea                         | Optional    |
+----------------+-----------------------------------------------------------------------------------------------+-------------+
| ``--hmmscan``  | Path to hmmscan output table containing the respective kingdom single-copy marker annotations | Required    |
+----------------+-----------------------------------------------------------------------------------------------+-------------+
| ``--out``      | Path to write filtered annotated markers corresponding to kingdom                             | Required    |
+----------------+-----------------------------------------------------------------------------------------------+-------------+
| ``--parallel`` | Use hmmscan parallel option (default: False)                                                  | Optional    |
+----------------+-----------------------------------------------------------------------------------------------+-------------+
| ``--cpus``     | Number of CPUs to use (default is to use all available CPUs)                                  | Optional    |
+----------------+-----------------------------------------------------------------------------------------------+-------------+
| ``--seed``     | Seed to set random state for hmmscan. (default: 42)                                           | Optional    |
+----------------+-----------------------------------------------------------------------------------------------+-------------+

You can view the complete command-line options using ``autometa-markers -h``

The above command would generate the following files:

+------------------------------+---------------------------------------------------------------------------------------+
| File                         | Description                                                                           |
+==============================+=======================================================================================+
| 78mbp_metagenome.hmmscan.tsv | hmmscan output table containing the respective kingdom single-copy marker annotations |
+------------------------------+---------------------------------------------------------------------------------------+
| 78mbp_metagenome.markers.tsv | Annotated marker table corresponding to the particular kingdom                        |
+------------------------------+---------------------------------------------------------------------------------------+

5. Taxonomy assignment
----------------------

5.1 BLASTP
^^^^^^^^^^

Autometa assigns a taxonomic rank to each contig and then takes only the contig belonging to the specified kingdom (either bacteria or archaea) for binning.
We found that in host-associated metagenomes, this step vastly improves the binning performance of Autometa (and other pipelines) because less eukaryotic
or viral contigs will be placed into bacterial bins.

The first step for contig taxonomy assignment is a local alignment search of the ORFs against a reference database. This can be accelerated using `diamond <https://github.com/bbuchfink/diamond>`_.

Create a diamond formatted database of the NCBI non-redundant (nr.gz) protein database.

.. code-block:: bash

    diamond makedb \
        --in $HOME/Autometa/autometa/databases/ncbi/nr.gz \
        --db $HOME/Autometa/autometa/databases/ncbi/nr \
        --threads 40

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

.. note::

    ``diamond makedb`` will append ``.dmnd`` to the provided path of ``--db``.

    i.e. ``--db /path/to/nr`` will become ``/path/to/nr.dmnd``

Run diamond blastp using the following command:

.. code-block:: bash

    diamond blastp \
        --query $HOME/tutorial/78mbp_metagenome.orfs.faa \
        --db $HOME/Autometa/autometa/databases/ncbi/nr.dmnd \
        --evalue 1e-5 \
        --max-target-seqs 200 \
        --threads 40 \
        --outfmt 6 \
        --out $HOME/tutorial/78mbp_metagenome.blastp.tsv

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

To see the complete list of acceptable output formats see Diamond `GitHub Wiki <https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#output-options>`__. A complete list of all command-line options for Diamond can be found on its `GitHub Wiki <https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options>`__.

.. caution::

    Autometa only parses output format 6 provided above as: ``--outfmt 6``

The above command would generate the blastP table (``78mbp_metagenome.blastp.tsv``) in output format 6

5.2 Lowest Common Ancestor (LCA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The second step in taxon assignment is determining each ORF's lowest common ancestor (LCA).
This step uses the blastp results generated in the previous step to generate a table having the LCA of each ORF. As a default only
the blastp hits (subject accessions) which are within 10% of the top bitscore are used. These subject accessions are translated to
their respective taxids (``prot.accession2taxid.gz``) to be looked up in NCBI's taxonomy database (``nodes.dmp``). Each ORFs' list of taxids
are then reduced to its lowest common ancestor via a range minimum query.

.. note::

    For more details on the range minimum query algorithm, see `the closed issue (#170) on Github <https://github.com/KwanLab/Autometa/issues/170>`_
    and a `walkthrough on topcoder <https://www.topcoder.com/thrive/articles/Range%20Minimum%20Query%20and%20Lowest%20Common%20Ancestor>`_



Use the following command to get the LCA of each ORF:

.. code-block:: bash

    autometa-taxonomy-lca \
        --blast $HOME/tutorial/78mbp_metagenome.blastp.tsv \
        --dbdir $HOME/Autometa/autometa/databases/ncbi/ \
        --lca-output $HOME/tutorial/78mbp_metagenome.lca.tsv \
        --sseqid2taxid-output $HOME/tutorial/78mbp_metagenome.lca.sseqid2taxid.tsv \
        --lca-error-taxids $HOME/tutorial/78mbp_metagenome.lca.errorTaxids.tsv

Let us dissect the above command:

+---------------------------+-------------------------------------------------------------------------------------------+----------------+
| Parameter                 | Function                                                                                  | Required (Y/N) |
+===========================+===========================================================================================+================+
| ``--blast``               | Path to diamond blastp output                                                             | Y              |
+---------------------------+-------------------------------------------------------------------------------------------+----------------+
| ``--dbdir``               | Path to NCBI databases directory                                                          | Y              |
+---------------------------+-------------------------------------------------------------------------------------------+----------------+
| ``--lca-output``          | Path to write lca output                                                                  | Y              |
+---------------------------+-------------------------------------------------------------------------------------------+----------------+
| ``--sseqid2taxid-output`` | Path to write qseqids sseqids to taxids translations table                                | N              |
+---------------------------+-------------------------------------------------------------------------------------------+----------------+
| ``--lca-error-taxids``    | Path to write table of blast table qseqids that were assigned root due to a missing taxid | N              |
+---------------------------+-------------------------------------------------------------------------------------------+----------------+

You can view the complete command-line options using ``autometa-taxonomy-lca -h``

The above command would generate a table (``78mbp_metagenome.lca.tsv``) having the name, rank and taxid of the LCA for each ORF.

5.3 Majority vote
^^^^^^^^^^^^^^^^^

The next step in taxon assignment is doing a modified majority vote to decide the taxonomy of each contig. This was developed to help minimize the effect of horizontal gene transfer (HGT). Briefly, the voting system helps assign the correct taxonomy to the contig from its component ORF classification. Even with highly divergent ORFs this allows for accurate kingdom level classification, enabling us to remove any eukaryotic contaminants or host DNA.

You can run the majority vote step using the following command:

.. code-block:: bash

    autometa-taxonomy-majority-vote \
        --lca $HOME/tutorial/78mbp_metagenome.lca.tsv \
        --output $HOME/tutorial/78mbp_metagenome.votes.tsv \
        --dbdir $HOME/Autometa/autometa/databases/ncbi/

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

You can view the complete command-line options using ``autometa-taxonomy-majority-vote -h``

The above command would generate a table (``78mbp_metagenome.votes.tsv``) having the taxid of each contig identified as per majority vote.

5.4 Split kingdoms
^^^^^^^^^^^^^^^^^^

In this final step of taxon assignment we use the voted taxid of each contig to split the contigs into different kingdoms and write them as per the provided canonical rank.

.. code-block:: bash

    autometa-taxonomy \
        --votes $HOME/tutorial/78mbp_metagenome.votes.tsv \
        --output $HOME/tutorial/ \
        --assembly $HOME/tutorial/78mbp_metagenome.filtered.fna \
        --prefix 78mbp_metagenome \
        --split-rank-and-write superkingdom \
        --ncbi $HOME/Autometa/autometa/databases/ncbi/

Let us dissect the above command:

+----------------------------+--------------------------------------------------------------------------------+-------------+
| Flag                       | Function                                                                       | Requirement |
+============================+================================================================================+=============+
| ``--votes``                | Path to voted taxids table                                                     | Required    |
+----------------------------+--------------------------------------------------------------------------------+-------------+
| ``--output``               | Directory to output fasta files of split canonical ranks and taxonomy.tsv      | Required    |
+----------------------------+--------------------------------------------------------------------------------+-------------+
| ``--assembly``             | Path to filtered metagenome assembly                                           | Required    |
+----------------------------+--------------------------------------------------------------------------------+-------------+
| ``--prefix``               | prefix to use for each file written                                            | Optional    |
+----------------------------+--------------------------------------------------------------------------------+-------------+
| ``--split-rank-and-write`` | Split contigs by provided canonical-rank column then write to output directory | Optional    |
+----------------------------+--------------------------------------------------------------------------------+-------------+
| ``--ncbi``                 | Path to ncbi database directory                                                | Optional    |
+----------------------------+--------------------------------------------------------------------------------+-------------+

Other options available for ``--split-rank-and-write`` are ``phylum``, ``class``, ``order``, ``family``, ``genus`` and ``species``

If ``--split-rank-and-write`` is specified then it will split contigs by provided canonical-rank column then write a file corresponding that rank. Eg. Bacteria.fasta, Archaea.fasta, etc for ``superkingdom``.

You can view the complete command-line options using ``autometa-taxonomy -h``

+-----------------------------------+------------------------------------------------------------------------------------------+
| File                              | Description                                                                              |
+===================================+==========================================================================================+
| 78mbp_metagenome.taxonomy.tsv     | Table with taxonomic classification of each contig                                       |
+-----------------------------------+------------------------------------------------------------------------------------------+
| 78mbp_metagenome.bacteria.fna     | Fasta file having the nucleic acid sequence of all bacterial contigs                     |
+-----------------------------------+------------------------------------------------------------------------------------------+
| 78mbp_metagenome.unclassified.fna | Fasta file having the nucleic acid sequence of all contigs unclassified at kingdom level |
+-----------------------------------+------------------------------------------------------------------------------------------+

In my case there are no non-bacterial contigs. For other datasets, ``autometa-taxonomy`` may produce other fasta files, for example Eukaryota.fasta and Viruses.fasta.

6. K-mer counting
-----------------

A k-mer (`ref <https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/>`_) is just a sequence of k characters in a string (or nucleotides in a DNA sequence). It is known that contigs that belong to the same genome have similar k-mer composition (`ref1 <https://sfamjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1462-2920.2004.00624.x?sid=nlm%3Apubmed>`_ and `ref2 <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-8-r85>`_) . Here, we compute k-mer frequencies of only the bacterial contigs.

This step does the following:

#. Create a k-mer count matrix of :math:`k^4/2` dimensions using the specified k-mer length
#. Normalization of the k-mer count matrix to a normalized k-mer frequency matrix
#. Reduce the dimensions of k-mer frequencies using principal component analysis (PCA).
#. Embed the PCA dimensions into two dimensions to allow the ease of visualization and manual binning of the contigs (see `ViZBin <https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-014-0066-1>`_ paper).

Use the following command to run the k-mer counting step:

.. code-block:: bash

    autometa-kmers \
        --fasta $HOME/tutorial/78mbp_metagenome.bacteria.fna \
        --kmers $HOME/tutorial/78mbp_metagenome.bacteria.kmers.tsv \
        --size 5 \
        --norm-method am_clr \
        --norm-output $HOME/tutorial/78mbp_metagenome.bacteria.kmers.normalized.tsv \
        --pca-dimensions 50 \
        --embedding-method bhsne \
        --embedding-output $HOME/tutorial/78mbp_metagenome.bacteria.kmers.embedded.tsv \
        --cpus 40 \
        --seed 42

Let us dissect the above command:

+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| Flag                   | Input arguments                                                                                                          | Requirement |
+========================+==========================================================================================================================+=============+
| ``--fasta``            | Path to length filtered metagenome assembly                                                                              | Required    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| ``--kmers``            | Path to k-mer frequency table                                                                                            | Required    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| ``--size``             | k-mer size in bp (default 5bp)                                                                                           | Optional    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| ``--norm-output``      | Path to normalized k-mer table                                                                                           | Required    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| ``--norm-method``      | Normalization method to transform kmer counts prior to PCA and embedding (default am_clr). Choices : ilr, clr and am_clr | Optional    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| ``--pca-dimensions``   | Number of dimensions to reduce to PCA feature space after normalization and prior to embedding (default: 50)             | Optional    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| ``--embedding-output`` | Path to embedded k-mer table                                                                                             | Required    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| ``--embedding-method`` | Embedding method to reduce the k-mer frequencies. Choices: sksne, bhsne (default), umap, densmap and trimap.             | Optional    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| ``--cpus``             | Number of CPUs to use (default is to use all available CPUs)                                                             | Optional    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+
| ``--seed``             | Set random seed for dimension reduction determinism (default 42). Useful in replicating the results                      | Optional    |
+------------------------+--------------------------------------------------------------------------------------------------------------------------+-------------+

You can view the complete command-line options using ``autometa-kmers -h``

The above command generates the following files:

+---------------------------------------+--------------------------------------------------------+
| File                                  | Description                                            |
+=======================================+========================================================+
| 78mbp_metagenome.kmers.tsv            | Table with raw k-mer counts of each contig             |
+---------------------------------------+--------------------------------------------------------+
| 78mbp_metagenome.kmers.normalized.tsv | Table with normalized k-mer frequencies of each contig |
+---------------------------------------+--------------------------------------------------------+
| 78mbp_metagenome.kmers.embedded.tsv   | Table with embedded k-mer frequencies of each contig   |
+---------------------------------------+--------------------------------------------------------+

.. _advanced-usage-kmers:

Advanced Usage
^^^^^^^^^^^^^^

In the command used above k-mer normalization is being done using Autometa's implementation of
the center log-ratio transform (am_clr). Other available normalization methods are isometric
log-ratio transform (ilr, scikit-bio implementation) and center log-ratio transform (clr, scikit-bio implementation).
Normalization method can be altered using the ``--norm-method`` flag.

In the above command k-mer embedding is being done using Barnes-Hut t-distributed Stochastic Neighbor Embedding (BH-tSNE).
Other embedding methods that are available are Uniform Manifold Approximation and Projection (UMAP), densMAP (a density-preserving tool based
on UMAP) and TriMap, a method that uses triplet constraints to form a low-dimensional embedding of a set of points.
Two implementations of BH-tSNE are available, ``bhsne`` and ``sksne`` corresponding to the tsne and scikit-learn libraries, respectively.
Embedding method can be altered using the ``--embedding-method`` flag.

Autometa uses a k-mer size of 5 and then embeds the resulting k-mer frequency table
into 50 PCA dimensions which are then reduced to two dimentions. k-mer size can be
altered using the ``--size`` flag, number of dimensions to reduce to PCA feature
space after normalization and prior to embedding can be altered using the ``--pca-dimensions``
flag and the number of dimensions of which to reduce k-mer frequencies can be altered using the ``--embedding-dimensions`` flag.

.. note::

    1. Even though ``bhsne`` and ``sksne`` are the same embedding method (but different implementations)
    they appear to give very different results. We recommend using the former.

    2. Providing a ``0`` to ``--pca-dimensions`` will skip the PCA step.

7. Binning
-----------

This is the step where contigs are binned into genomes via clustering.
Autometa assesses genome bins by examining their completeness, purity,
GC content std.dev. and coverage std.dev. A taxonomy table may also be used
to selectively iterate through contigs based on their profiled taxon.

This step does the following:

#. Optionally iterate through contigs based on taxonomy
#. Bin contigs based on embedded k-mer coordinates and coverage
#. Accept genome bins that pass the following metrics:
    #. Above completeness threshold (``default=20.0``)
    #. Above purity threshold (``default=95.0``)
    #. Below GC content standard deviation threshold (``default=5.0``)
    #. Below coverage standard deviation threshold (``default=25.0``)
#. Unbinned contigs will be re-binned until no more acceptable genome bins are yielded

If you include a taxonomy table Autometa will attempt to further partition the data based
on ascending taxonomic specificity (i.e. in the order superkingdom, phylum, class, order,
family, genus, species) when binning unclustered contigs from a previous attempt. We found
that this is mainly useful if you have a highly complex metagenome (lots of species), or
you have several related species at similar coverage level.

Use the following command to perform binning:

.. code-block:: bash

    autometa-binning \
        --kmers $HOME/tutorial/78mbp_metagenome.bacteria.kmers.embedded.tsv \
        --coverages $HOME/tutorial/78mbp_metagenome.coverages.tsv \
        --gc-content $HOME/tutorial/78mbp_metagenome.gc_content.tsv \
        --markers $HOME/tutorial/78mbp_metagenome.markers.tsv \
        --clustering-method dbscan \
        --completeness 20 \
        --purity 90 \
        --cov-stddev-limit 25 \
        --gc-stddev-limit 5 \
        --taxonomy $HOME/tutorial/78mbp_metagenome.taxonomy.tsv \
        --output-binning $HOME/tutorial/78mbp_metagenome.binning.tsv \
        --output-main $HOME/tutorial/78mbp_metagenome.main.tsv \
        --starting-rank superkingdom \
        --rank-filter superkingdom
        --rank-name-filter bacteria

Let us dissect the above command:

+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| Flag                    | Function                                                                                | Requirement |
+=========================+=========================================================================================+=============+
| ``--kmers``             | Path to embedded k-mer frequencies table                                                | Required    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--coverages``         | Path to metagenome coverages table                                                      | Required    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--gc-content``        | Path to metagenome GC contents table                                                    | Required    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--markers``           | Path to Autometa annotated markers table                                                | Required    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--output-binning``    | Path to write Autometa binning results                                                  | Required    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--output-main``       | Path to write Autometa main table                                                       | Required    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--clustering-method`` | Clustering algorithm to use for recursive binning. Choices dbscan (default) and hdbscan | Optional    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--completeness``      | completeness cutoff to retain cluster (default 20)                                      | Optional    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--purity``            | purity cutoff to retain cluster (default 95)                                            | Optional    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--cov-stddev-limit``  | coverage standard deviation limit to retain cluster (default 25)                        | Optional    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--gc-stddev-limit``   | GC content standard deviation limit to retain cluster (default 5)                       | Optional    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--taxonomy``          | Path to Autometa assigned taxonomies table                                              | Required    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--starting-rank``     | Canonical rank at which to begin subsetting taxonomy (default: superkingdom)            | Optional    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+
| ``--domain``            | Kingdom to consider. Choices bacteria (default) and archaea                             | Optional    |
+-------------------------+-----------------------------------------------------------------------------------------+-------------+

You can view the complete command-line options using ``autometa-binning -h``

The above command generates the following files:

#. ``78mbp_metagenome.binning.tsv`` contains the final binning results along with a few more metrics regarding each genome bin.
#. ``78mbp_metagenome.main.tsv`` which contains the feature table that was utilized during the genome binning process as well as the corresponding output predictions.

The following table describes each column for the resulting binning outputs. We'll start with the columns present in ``78mbp_metagenome.binning.tsv``
then describe the additional columns that are present in ``78mbp_metagenome.main.tsv``.

+-------------------+------------------------------------------------------------------------------------------------------------------------+
| Column            | Description                                                                                                            |
+===================+========================================================================================================================+
| Contig            | Name of the contig in the input fasta file                                                                             |
+-------------------+------------------------------------------------------------------------------------------------------------------------+
| Cluster           | Genome bin assigned by autometa to the contig                                                                          |
+-------------------+------------------------------------------------------------------------------------------------------------------------+
| Completeness      | Estimated completeness of the Genome bin, based on single-copy marker genes                                            |
+-------------------+------------------------------------------------------------------------------------------------------------------------+
| Purity            | Estimated purity of the Genome bin, based on the number of single-copy marker genes that are duplicated in the cluster |
+-------------------+------------------------------------------------------------------------------------------------------------------------+
| coverage_stddev   | Coverage standard deviation of the Genome bin                                                                          |
+-------------------+------------------------------------------------------------------------------------------------------------------------+
| gc_content_stddev | GC content standard deviation of the Genome bin                                                                        |
+-------------------+------------------------------------------------------------------------------------------------------------------------+

In addition to the above columns ``78mbp_metagenome.main.tsv`` file has the following additional columns:

+--------------+-------------------------------------------------+
| Column       | Description                                     |
+==============+=================================================+
| Coverage     | Estimated coverage of the contig                |
+--------------+-------------------------------------------------+
| gc_content   | Estimated GC content of the contig              |
+--------------+-------------------------------------------------+
| length       | Estimated length of the contig                  |
+--------------+-------------------------------------------------+
| species      | Assigned taxonomic species for the contig       |
+--------------+-------------------------------------------------+
| genus        | Assigned taxonomic genus for the contig         |
+--------------+-------------------------------------------------+
| family       | Assigned taxonomic family for the contig        |
+--------------+-------------------------------------------------+
| order        | Assigned taxonomic order for the contig         |
+--------------+-------------------------------------------------+
| class        | Assigned taxonomic class for the contig         |
+--------------+-------------------------------------------------+
| phylum       | Assigned taxonomic phylum for the contig        |
+--------------+-------------------------------------------------+
| superkingdom | Assigned taxonomic superkingdom for the contig  |
+--------------+-------------------------------------------------+
| taxid        | Assigned NCBI taxonomy ID number for the contig |
+--------------+-------------------------------------------------+
| x_1          | The first coordinate after dimension reduction  |
+--------------+-------------------------------------------------+
| x_2          | The second coordinate after dimension reduction |
+--------------+-------------------------------------------------+

You can attempt to improve your genome bins with an unclustered recruitment step which uses features from existing genome bins to recruit unbinned contigs.
Alternatively you can use these initial genome bin predictions and continue to the :ref:`Examining Results` section.

.. _advanced-usage-binning:

Advanced Usage
^^^^^^^^^^^^^^

.. code-block::

    Completeness = Number of single copy marker genes present just once / Total number of single copy marker genes

    Purity = Number of single copy marker genes present more than once / Total number of single copy marker genes

These are default parameters that autometa uses to accept clusters are 20% complete, 95% pure, below 25% coverage standard deviation
and below 5% GC content standard deviation. These parameters can be altered using the flags, ``--completeness``, ``--purity``, ``--cov-stddev-limit`` and ``--gc-stddev-limit``.

There are two binning algorithms to choose from Density-Based Spatial Clustering of Applications with Noise (`DBSCAN <https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html>`_)
and Hierarchical Density-Based Spatial Clustering of Applications with Noise (`HDBSCAN <https://hdbscan.readthedocs.io/en/latest/index.html>`_). The default is DBSCAN.

It is important to note that if recursively binning with taxonomy, only contigs at the specific taxonomic rank are analyzed and once the binning algorithm has moved on to the next rank, these
are not considered until they fall under another taxonomic rank under consideration. I.e. Iterate through phyla. Contig of one phylum is only considered for that phylum then not
for the rest of the phyla. If it is still unbinned at the Class rank, then it will be considered only during its respective Class's iteration. The canonical rank from which to start
binning can be changed using the ``--starting-rank`` flag. The default is ``superkingdom``.

8. Unclustered recruitment (Optional)
-------------------------------------

An unclustered recruitment step which uses features from existing genome bins is used to classify the unbinned contigs to the genome bins that were produced in the previous step.
This step is optional and the results should be verified before proceeding with these results.

.. note::

    The machine learning step has been observed to bin contigs that do not necessarily belong to the predicted genome. Careful inspection of coverage and taxonomy should be done before proceeding with these results.

Use the following command to run the unclustered recruitment step:

.. code-block:: bash

    autometa-unclustered-recruitment \
        --kmers $HOME/tutorial/78mbp_metagenome.bacteria.kmers.normalized.tsv \
        --coverage $HOME/tutorial/78mbp_metagenome.coverages.tsv \
        --binning $HOME/tutorial/78mbp_metagenome.binning.tsv \
        --markers $HOME/tutorial/78mbp_metagenome.markers.tsv \
        --taxonomy $HOME/tutorial/78mbp_metagenome.taxonomy.tsv \
        --output-binning $HOME/tutorial/78mbp_metagenome.recruitment.binning.tsv \
        --output-features $HOME/tutorial/78mbp_metagenome.recruitment.features.tsv \
        --output-main $HOME/tutorial/78mbp_metagenome.recruitment.main.tsv \
        --classifier decision_tree \
        --seed 42

Let us dissect the above command:

+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| Flag                  | Function                                                                                        | Required (Y/N) |
+=======================+=================================================================================================+================+
| ``--kmers``           | Path to normalized k-mer frequencies table                                                      |        Y       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| ``--coverages``       | Path to metagenome coverages table                                                              |        Y       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| ``--binning``         | Path to autometa binning output                                                                 |        Y       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| ``--markers``         | Path to Autometa annotated markers table                                                        |        Y       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| ``--output-binning``  | Path to write Autometa unclustered recruitment table                                            |        Y       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| ``--taxonomy``        | Path to taxonomy table                                                                          |        N       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| ``--output-features`` | Path to write Autometa main table used during/after unclustered recruitment                     |        N       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| ``--output-main``     | Path to write Autometa main table used during/after unclustered recruitment                     |        N       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| ``--classifier``      | classifier to use for recruitment of contigs. Choices decision_tree (default) and random_forest |        N       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+
| ``--seed``            | Seed to use for RandomState when initializing classifiers (default: 42)                         |        N       |
+-----------------------+-------------------------------------------------------------------------------------------------+----------------+

You can view the complete command-line options using ``autometa-unclustered-recruitment -h``

The above command would generate ``78mbp_metagenome.recruitment.binning.tsv`` and ``78mbp_metagenome.recruitment.main.tsv``.

``78mbp_metagenome.recruitment.binning.tsv`` contains the final predictions of ``autometa-unclustered-recruitment``. ``78mbp_metagenome.recruitment.features.tsv``
is the feature table utilized during/after the unclustered recruitment algorithm. This represents unbinned contigs with their respective annotations and output predictions of their recruitment into a genome bin.
The taxonomic features have been encoded using “one-hot encoding” or a presence/absence matrix where each column is a canonical taxonomic rank and its respective value for each row represents its presence or absence.
Presence and absence are denoted with 1 and 0, respectively. Hence "one-hot" encoding being an encoding of presence and absence of the respective annotation type. In our case taxonomic designation.

The ``78mbp_metagenome.recruitment.binning.tsv`` file contains the following columns:

+-------------------+----------------------------------------------------------------------------------+
| Column            | Description                                                                      |
+===================+==================================================================================+
| contig            | Name of the contig in the input fasta file                                       |
+-------------------+----------------------------------------------------------------------------------+
| cluster           | Genome bin assigned by autometa to the contig                                    |
+-------------------+----------------------------------------------------------------------------------+
| recruited_cluster | Genome bin assigned by autometa to the contig after unclustered recruitment step |
+-------------------+----------------------------------------------------------------------------------+

.. _advanced-usage-unclustered-recruitment:

Advanced Usage
^^^^^^^^^^^^^^

The clustering method for the unclustered recruitment step can be performed either using a decision tree classifier (default) or using a random forst algorithm. The choice of method can be selected using the  ``--classifier`` flag.

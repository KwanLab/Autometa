************
Benchmarking
************

.. note::

    The most recent benchmarking results are hosted on our `KwanLab/metaBenchmarks <https://github.com/KwanLab/metaBenchmarks>`_ Github repository
    and provide a range of analyses covering multiple stages and parameter sets. These benchmarks are available with their own respective
    modules so that the community may easily assess how their novel (``taxon-profiling``, ``clustering``, ``binning``, ``refinement``) algorithms
    perform compared to current state-of-the-art methods. Tools were selected for benchmarking based on their relevance
    to environmental, single-assembly, reference-free binning pipelines.

Example benchmarking with simulated communities
===============================================

Downloading Test Datasets
^^^^^^^^^^^^^^^^^^^^^^^^^

The first step in benchmarking is to download the test data files. We will be benchmarking the taxon-profiling step and binning steps and thus the files for the same will be downloaded. reference_assignments are needed as something to compare out output with. Same reference_assignments file will be used for both benchmarking.

Autometa is packaged with a built-in module that allows any user to download any of the available test datasets.
To use these utilities simply run the command ``autometa-download-dataset``.

For example, to download all of the simulated communities reference binning/taxonomy assignments as well as the Autometa
v2.0 binning/taxonomy predictions:

.. code:: bash

    # Note: community is the test dataset that was used for clustering or classification. e.g.
    # choices: 78Mbp,156Mbp,312Mbp,625Mbp,1250Mbp,2500Mbp,5000Mbp,10000Mbp
    community_sizes=(78Mbp 156Mbp 312Mbp 625Mbp 1250Mbp 2500Mbp 5000Mbp 10000Mbp)

    autometa-download-dataset \
    --community-type simulated \
    --community-sizes ${community_sizes[@]} \
    --file-names reference_assignments.tsv.gz binning.tsv.gz taxonomy.tsv.gz \
    --dir-path simulated

Benchmark clustering
^^^^^^^^^^^^^^^^^^^^

Here we are benchmarking the clustering step. All the simulated communities are being benchmarked.

.. code:: bash

    for community_size in ${community_sizes[@]};do
        autometa-benchmark \
            --benchmark clustering \
            --predictions simulated/${community_size}/binning.tsv.gz \
            --reference simulated/${community_size}/reference_assignments.tsv.gz \
            --output-wide ${community_size}.clustering_benchmarks.wide.tsv.gz \
            --output-long ${community_size}.clustering_benchmarks.long.tsv.gz
    done

Aggregate resulst across simulated communities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we have provided some convinient python scripts for you to aggregate multiple binning results (whether they have unique dataset values or not).

When dataset index is unique
----------------------------

When all the dataset files are unique the index column is the dataset column.

.. code:: python

    import pandas as pd
    import glob
    df = pd.concat([
        pd.read_csv(fp, sep="\t", index_col="dataset")
        for fp in glob.glob("*.clustering_benchmarks.long.tsv.gz")
    ])
    df.to_csv("benchmarks.tsv", sep='\t', index=True, header=True)

When dataset index is `not` unique
----------------------------------

When the dataset files are not unique, then the index (dataset value) will be renamed to the filename of the respective file. 

.. code:: python

    import pandas as pd
    import os
    import glob
    dfs = []
    for fp in glob.glob("*.clustering_benchmarks.long.tsv.gz"):
        df = pd.read_csv(fp, sep="\t", index_col="dataset")
        df.index = df.index.map(lambda fpath: os.path.basename(fpath))
        dfs.append(df)
    df = pd.concat(dfs)
    df.to_csv("benchmarks.tsv", sep='\t', index=True, header=True)

Benchmark classification
^^^^^^^^^^^^^^^^^^^^^^^^

Previously we benchmarked the clustering step. Now we are going to benchmark the taxon profiling step. The only difference is that the ``--classification`` flag is provide instead of the ``--clustering`` flag.

.. code:: bash

    autometa-benchmark \
        --benchmark classification \
        --predictions taxonomy.tsv.gz \
        --reference simulated/${community_size}/reference_assignments.tsv.gz \
        --output-wide ${community_size}.classification_benchmarks.wide.tsv.gz \
        --output-classification-reports taxa_reports \
        --ncbi /ncbi

Benchmark clustering-classification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This will benchmark the final binning results using the classification metrics, 
precision, recall and F1-score.

.. note::
    
    First, F1-scores are determined for each reference genome across all recovered genome bins (specifically for clusters containing any of the reference genome's contigs). Then, each reference genome's F1-score is selected corresponding to the respective cluster with the largest set of the respective reference genome's contigs (by length in bp).

.. code:: bash
    
    autometa-download-dataset \
        --community-type simulated \
        --community-sizes 78Mbp \
        --file-names binning.tsv.gz\
        --dir-path simulated

Now running the command to benchmark:

.. code:: bash

    autometa-benchmark \
        --benchmark binning-classification \
        --predictions 78Mbp/binning.tsv \
        --reference 78Mbp/reference_assignments.tsv.gz \
        --output-wide binning_classification_wide.tsv.gz

Specifying multiple results
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``autometa-benchmark`` provides the flexibility of specifying multiple output files of the same dataset and benchmarking them simultaneously. For example, in case you have done the clustering for a sample using three different binner, you can benchmark all of them simulateously. This is really handy when benchmarking multiple tools on the same dataset.

For clustering
--------------

.. code:: bash

    autometa-benchmark \
       --benchmark clustering \
       --predictions binner1_78Mbp_output.tsv.gz binner2_78Mbp_output.tsv.gz binner3_78Mbp_output.tsv.gz \
       --reference simulated/78Mbp/reference_assignments.tsv.gz \
       --output-wide 78Mbp.clustering_benchmarks.wide.tsv.gz \
       --output-long 78Mbp.clustering_benchmarks.long.tsv.gz

For classification
------------------

.. code:: bash

    autometa-benchmark \
        --benchmark classification \
        --predictions taxonomic_profiler1_78Mbp.tsv.gz taxonomic_profiler2_78Mbp.tsv.gz taxonomic_profiler3_78Mbp.tsv.gz  \
        --reference simulated/78Mbp/reference_assignments.tsv.gz \
        --output-wide 78Mbp.classification_benchmarks.wide.tsv.gz \
        --output-classification-reports taxa_reports \
        --ncbi /ncbi

Autometa Test Datasets
======================

Simulated Communities
^^^^^^^^^^^^^^^^^^^^^

.. csv-table:: Autometa Simulated Communities
    :file: simulated_community.csv
    :header-rows: 1

You can download all the Simulated communities using this `link <https://drive.google.com/drive/folders/1JFjVb-pfQTv4GXqvqRuTOZTfKdT0MwhN?usp=sharing>`__.
Individual communities can be downloaded using the links in the above table.

For more information on simulated communities,
check the `README.md <https://drive.google.com/file/d/1Ti05Qp13FleuMQdnp3C5L-sXnIM25EZE/view?usp=sharing>`__
located in the ``simulated_communities`` directory.

Generating New Simulated Communities
------------------------------------

Communities were simulated using `ART <https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm>`__,
a sequencing read simulator, with a collection of 3000 bacteria randomly retrieved.
Genomes were retrieved until the provided total length was reached.

e.g. ``-l 1250`` would translate to 1250Mbp as the sum of total lengths for all bacterial genomes retrieved.

.. code:: bash

    # Work out coverage level for art_illumina
    # C = [(LN)/G]/2
    # C = coverage
    # L = read length (total of paired reads)
    # G = genome size in bp
    # -p  : indicate a paired-end read simulation or to generate reads from both ends of amplicons
    # -ss : HS25 -> HiSeq 2500 (125bp, 150bp)
    # -f  : fold of read coverage simulated or number of reads/read pairs generated for each amplicon
    # -m  : the mean size of DNA/RNA fragments for paired-end simulations
    # -s  : the standard deviation of DNA/RNA fragment size for paired-end simulations.
    # -l  : the length of reads to be simulated
    $ coverage = ((250 * reads) / (length * 1000000))
    $ art_illumina -p -ss HS25 -l 125 -f $coverage -o simulated_reads -m 275 -s 90 -i asm_path

Synthetic Communities
^^^^^^^^^^^^^^^^^^^^^

51 bacterial isolates were assembled into synthetic communities which we've titled ``MIX51``.

The initial synthetic community was prepared using a mixture of fifty-one bacterial isolates.
The synthetic community's DNA was extracted for sequencing, assembly and binning.

You can download the MIX51 community using this `link <https://drive.google.com/drive/folders/1x8d0o6HO5N72j7p_D_YxrSurBfpi9zmK?usp=sharing>`__.

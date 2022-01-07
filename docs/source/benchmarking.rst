============
Benchmarking
============

.. headers/sections formatting
    See: https://docs.typo3.org/m/typo3/docs-how-to-document/main/en-us/WritingReST/HeadlinesAndSection.html

.. note::

    The most recent Autometa benchmarking results covering multiple modules and input parameters are hosted on our 
    `KwanLab/metaBenchmarks <https://github.com/KwanLab/metaBenchmarks>`_ Github repository and provide a range of 
    analyses covering multiple stages and parameter sets. These benchmarks are available with their own respective
    modules so that the community may easily assess how Autometa's novel (``taxon-profiling``, ``clustering``, 
    ``binning``, ``refinement``) algorithms perform compared to current state-of-the-art methods. Tools were selected for 
    benchmarking based on their relevance to environmental, single-assembly, reference-free binning pipelines.

Benchmarking with the ``autometa-benchmark`` module
===================================================

Autometa includes the ``autometa-benchmark`` entrypoint, a script to benchmark Autometa taxon-profiling, clustering
and binning-classification prediction results using clustering and classification evaluation metrics. To select the
appropriate benchmarking method, supply the ``--benchmark`` parameter with the respective choice. The three benchmarking
methods are detailed below.

.. note::
    If you'd like to follow along with the benchmarking commands, you may download the test datasets
    using:

    .. code:: bash

        autometa-download-dataset \
            --community-type simulated \
            --community-sizes 78Mbp \
            --file-names reference_assignments.tsv.gz binning.tsv.gz taxonomy.tsv.gz \
            --dir-path $HOME/Autometa/autometa/datasets/simulated

    This will download three files:

    - ``reference_assignments``: tab-delimited file containing contigs with their reference genome assignments. ``cols: [contig, reference_genome, taxid, organism_name, ftp_path, length]``
    - ``binning.tsv.gz``: tab-delimited file containing contigs with Autometa binning predictions, ``cols: [contig, cluster]``
    - ``taxonomy.tsv.gz``: tab-delimited file containing contigs with Autometa taxon-profiling predictions ``cols: [contig, kingdom, phylum, class, order, family, genus, species, taxid]``


Taxon-profiling
---------------

Example benchmarking with simulated communities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    # Set community size (see above for selection/download of other community types)
    community_size=78Mbp
    
    # Inputs
    ## NOTE: predictions and reference were downloaded using autometa-download-dataset
    predictions="$HOME/Autometa/autometa/datasets/simulated/${community_size}/taxonomy.tsv.gz" # required columns -> contig, taxid
    reference="$HOME/Autometa/autometa/datasets/simulated/${community_size}/reference_assignments.tsv.gz"
    ncbi=$HOME/Autometa/autometa/databases/ncbi

    # Outputs
    output_wide="${community_size}.taxon_profiling_benchmarks.wide.tsv.gz" # file path
    output_long="${community_size}.taxon_profiling_benchmarks.long.tsv.gz" # file path
    reports="${community_size}_taxon_profiling_reports" # directory path

    autometa-benchmark \
        --benchmark classification \
        --predictions $predictions \
        --reference $reference \
        --ncbi $ncbi \
        --output-wide $output_wide \
        --output-long $output_long \
        --output-classification-reports $reports

.. note::
    Using ``--benchmark=classification`` requires the path to a directory containing files (nodes.dmp, names.dmp, merged.dmp) 
    from NCBI's taxdump tarball. This should be supplied using the ``--ncbi`` parameter.

Clustering
----------

Example benchmarking with simulated communities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    # Set community size (see above for selection/download of other community types)
    community_size=78Mbp

    # Inputs
    ## NOTE: predictions and reference were downloaded using autometa-download-dataset
    predictions="$HOME/Autometa/autometa/datasets/simulated/${community_size}/binning.tsv.gz" # required columns -> contig, cluster
    reference="$HOME/Autometa/autometa/datasets/simulated/${community_size}/reference_assignments.tsv.gz"

    # Outputs
    output_wide="${community_size}.clustering_benchmarks.wide.tsv.gz"
    output_long="${community_size}.clustering_benchmarks.long.tsv.gz"
    
    autometa-benchmark \
        --benchmark clustering \
        --predictions $predictions \
        --reference $reference \
        --output-wide $output_wide \
        --output-long $output_long


Binning
-------

Example benchmarking with simulated communities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    # Set community size (see above for selection/download of other community types)
    community_size=78Mbp
    
    # Inputs
    ## NOTE: predictions and reference were downloaded using autometa-download-dataset
    predictions="$HOME/Autometa/autometa/datasets/simulated/${community_size}/binning.tsv.gz" # required columns -> contig, cluster
    reference="$HOME/Autometa/autometa/datasets/simulated/${community_size}/reference_assignments.tsv.gz"
    
    # Outputs
    output_wide="${community_size}.binning_benchmarks.wide.tsv.gz"
    output_long="${community_size}.binning_benchmarks.long.tsv.gz"
    
    autometa-benchmark \
        --benchmark binning-classification \
        --predictions $predictions \
        --reference $reference \
        --output-wide $output_wide \
        --output-long $output_long


Autometa Test Datasets
======================

Descriptions
------------

Simulated Communities
~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Autometa Simulated Communities
    :file: simulated_community.csv
    :header-rows: 1

You can download all the Simulated communities using this `link <https://drive.google.com/drive/folders/1JFjVb-pfQTv4GXqvqRuTOZTfKdT0MwhN?usp=sharing>`__.
Individual communities can be downloaded using the links in the above table.

For more information on simulated communities,
check the `README.md <https://drive.google.com/file/d/1Ti05Qp13FleuMQdnp3C5L-sXnIM25EZE/view?usp=sharing>`__
located in the ``simulated_communities`` directory.

Synthetic Communities
~~~~~~~~~~~~~~~~~~~~~

51 bacterial isolates were assembled into synthetic communities which we've titled ``MIX51``.

The initial synthetic community was prepared using a mixture of fifty-one bacterial isolates.
The synthetic community's DNA was extracted for sequencing, assembly and binning.

You can download the MIX51 community using this `link <https://drive.google.com/drive/folders/1x8d0o6HO5N72j7p_D_YxrSurBfpi9zmK?usp=sharing>`__.

Download
--------

Using ``autometa-download-dataset``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Autometa is packaged with a built-in module that allows any user to download any of the available test datasets.
To use retrieve these datasets one simply needs to run the ``autometa-download-dataset`` command.

For example, to download the reference assignments for a simulated community as well as the most recent Autometa
binning and taxon-profiling predictions for this community, provide the following parameters: 

.. code:: bash

    # choices for simulated: 78Mbp,156Mbp,312Mbp,625Mbp,1250Mbp,2500Mbp,5000Mbp,10000Mbp
    autometa-download-dataset \
        --community-type simulated \
        --community-sizes 78Mbp \
        --file-names reference_assignments.tsv.gz binning.tsv.gz taxonomy.tsv.gz \
        --dir-path simulated


This will download ``reference_assignments.tsv.gz``, ``binning.tsv.gz``, ``taxonomy.tsv.gz`` to the ``simulated/78Mbp`` directory.

- ``reference_assignments``: tab-delimited file containing contigs with their reference genome assignments. ``cols: [contig, reference_genome, taxid, organism_name, ftp_path, length]``
- ``binning.tsv.gz``: tab-delimited file containing contigs with Autometa binning predictions, ``cols: [contig, cluster]``
- ``taxonomy.tsv.gz``: tab-delimited file containing contigs with Autometa taxon-profiling predictions ``cols: [contig, kingdom, phylum, class, order, family, genus, species, taxid]``

Using ``gdrive``
----------------

You can download the individual assemblies of different datasests with the help of ``gdown`` using command line
(This is what ``autometa-download-dataset`` is using behind the scenes). If you have installed ``autometa`` using
``conda`` then ``gdown`` should already be installed. If not, you can install it using 
``conda install -c conda-forge gdown`` or ``pip install gdown``.

Example for the 78Mbp simulated community
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Navigate to the 78Mbp community dataset using the `link <https://drive.google.com/drive/u/2/folders/1McxKviIzkPyr8ovj8BG7n_IYk-QfHAgG>`_ mentioned above.
2. Get the file ID by navigating to any of the files and right clicking, then selecting the ``get link`` option. 
    This will have a ``copy link`` button that you should use. The link for the metagenome assembly 
    (ie. ``metagenome.fna.gz``) should look like this : ``https://drive.google.com/file/d/15CB8rmQaHTGy7gWtZedfBJkrwr51bb2y/view?usp=sharing``
3. The file ID is within the ``/`` forward slashes between ``file/d/`` and ``/``, e.g:

.. code:: bash

    # Pasted from copy link button:
    https://drive.google.com/file/d/15CB8rmQaHTGy7gWtZedfBJkrwr51bb2y/view?usp=sharing
    #                 begin file ID ^ ------------------------------^ end file ID

4. Copy the file ID
5. Now that we have the File ID, you can specify the ID or use the ``drive.google.com`` prefix. Both should work.

.. code:: bash

    file_id="15CB8rmQaHTGy7gWtZedfBJkrwr51bb2y"
    gdown --id ${file_id} -O metagenome.fna.gz
    # or
    gdown https://drive.google.com/uc?id=${file_id} -O metagenome.fna.gz

.. note::

    Unfortunately, at the moment ``gdown`` doesn't support downloading entire directories from Google drive.
    There is an open `Pull request <https://github.com/wkentaro/gdown/pull/90#issue-569060398>`_ on the ``gdown`` repository
    addressing this specific issue which we are keeping a close eye on and will update this documentation when it is merged.


Advanced
========

Data Handling
-------------

Aggregating benchmarking results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When dataset index is unique
""""""""""""""""""""""""""""

.. code:: python

    import pandas as pd
    import glob
    df = pd.concat([
        pd.read_csv(fp, sep="\t", index_col="dataset")
        for fp in glob.glob("*.clustering_benchmarks.long.tsv.gz")
    ])
    df.to_csv("benchmarks.tsv", sep='\t', index=True, header=True)

When dataset index is `not` unique
""""""""""""""""""""""""""""""""""

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


Downloading multiple test datasets at once
------------------------------------------

To download all of the simulated communities reference binning/taxonomy assignments as well as the Autometa
v2.0 binning/taxonomy predictions all at once, you can provide the multiple arguments to ``--community-sizes``.

e.g. ``--community-sizes 78Mbp 156Mbp 312Mbp 625Mbp 1250Mbp 2500Mbp 5000Mbp 10000Mbp``

An example of this is shown in the bash script below:

.. code:: bash

    # choices: 78Mbp,156Mbp,312Mbp,625Mbp,1250Mbp,2500Mbp,5000Mbp,10000Mbp
    community_sizes=(78Mbp 156Mbp 312Mbp 625Mbp 1250Mbp 2500Mbp 5000Mbp 10000Mbp)

    autometa-download-dataset \
        --community-type simulated \
        --community-sizes ${community_sizes[@]} \
        --file-names reference_assignments.tsv.gz binning.tsv.gz taxonomy.tsv.gz \
        --dir-path simulated

Generating new simulated communities
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
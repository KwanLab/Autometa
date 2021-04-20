============
Benchmarking
============

This page contains information regarding test datasets as well as benchmarking statistics. The current Autometa pipeline was compared against its previous version as well other binning pipelines.

Test datasets
=============

Simulated
---------

Communities were simulated using `ART <https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm>`__, a sequencing read simulator, with a collection of 3000 bacteria randomly retrieved. Genomes were retrieved until the provided total length was reached.

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

.. csv-table:: Simulated communities
    :file: simulated_community.csv
    :header-rows: 1

You can download all the Simulated communities using this `link <https://drive.google.com/drive/folders/1JFjVb-pfQTv4GXqvqRuTOZTfKdT0MwhN?usp=sharing>`__. Individual communities can be downloaded using the links in the above table.

For more information on simulated communities, check the `README.md <https://drive.google.com/file/d/1Ti05Qp13FleuMQdnp3C5L-sXnIM25EZE/view?usp=sharing>`__ located in the `simulated_communities` directory.

Synthetic
---------

51 bacterial isolates were assembled into synthetic communities which we've titled ``MIX51``.

The initial synthetic community was prepared using a mixture of fifty-one bacterial isolates. The synthetic community's DNA was extracted for sequencing, assembly and binning.

You can download the MIX51 community using this `link <https://drive.google.com/drive/folders/1x8d0o6HO5N72j7p_D_YxrSurBfpi9zmK?usp=sharing>`__.

Download datasests
==================

Using autometa built-in module
------------------------------

.. todo::
    Address `Issue #110 <https://github.com/KwanLab/Autometa/issues/110#issue-707373779>`_ and add steps here.


Using command line
------------------

You can download the individual assemblies of different datasests with the help of ``gdown`` using command line. If you have installed autometa using ``conda`` then ``gdown`` should already be installed. If not, you can install it using ``conda install -c conda-forge gdown`` or ``pip install gdown``.

`Example for the 78Mbp community`

1. Navigate to the 78Mbp community dataset using the `link <https://drive.google.com/drive/u/2/folders/1McxKviIzkPyr8ovj8BG7n_IYk-QfHAgG>`_ mentioned above.
2. Get the file ID by navigating to any of the files and right clicking, then selecting the ``get link`` option. This will have a ``copy link`` button that you should use. The link for the metagenome assembly (ie. metagenome.fna.gz) should look like this : ``https://drive.google.com/file/d/15CB8rmQaHTGy7gWtZedfBJkrwr51bb2y/view?usp=sharing``
3. The file ID is within the / forward slashes between file/d/ and /, e.g:

.. code:: bash

    # Pasted from copy link button:
    https://drive.google.com/file/d/15CB8rmQaHTGy7gWtZedfBJkrwr51bb2y/view?usp=sharing
    #                 begin file ID ^ ------------------------------^ end file ID

4. Copy the file ID
5. Now that we have the File ID, you can specify the ID or use the drive.google.com prefix. Both should work.

.. code:: bash

    file_id="15CB8rmQaHTGy7gWtZedfBJkrwr51bb2y" 
    gdown --id ${file_id} -O metagenome.fna.gz
    # or
    gdown https://drive.google.com/uc?id=${file_id} -O metagenome.fna.gz

.. note:: 

    Unfortunately, at the moment ``gdown`` doesn't support downloading entire directories from Google drive. There is an open `Pull request <https://github.com/wkentaro/gdown/pull/90#issue-569060398>`_ on the ``gdown`` repository addressing this specific issue which we are keeping a close eye on and will update this documentation when it is merged.

Benchmarks
==========

.. todo::
    Add the Benchmarking statistics

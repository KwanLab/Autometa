============
Benchmarking
============

This page contains information regarding the test datasests the users can use as well as benchmarking of the current autometa algorithm against its previous version as well other binning pipelines.

Test datasests
==============



Simulated
---------

Communities were simulated using ART, a sequencing read simulator, with a collection of 3000 bacteria randomly retrieved. Genomes were retrieved until the provided total length was reached.

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

Benchmarks
==========

.. todo::
    Add the Benchmarking statistics

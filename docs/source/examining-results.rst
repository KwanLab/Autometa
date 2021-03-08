=================
Examining Results
=================

Automappa
=========

An interactive interface for exploration and refinement of metagenomes Automappa is a tool built to interface with Autometa output to help you explore your binning results.

For details, see the Automappa page `here. <https://github.com/WiscEvan/Automappa>`__

Binning summary
===============

Create a master table with all the binning annotations including embedded kmers, taxonomy and coverage.

.. code-block:: python

    import pandas as pd
    kmers = pd.read_csv("kmers.am_clr.bhsne.tsv", sep='\t', index_col='contig')
    binning = pd.read_csv("binning.hdbscan.tsv", sep='\t', index_col='contig')
    binning = pd.read_csv("taxonomy.tsv", sep='\t', index_col='contig')
    coverage=pd.read_csv("coverage.tsv", sep="\t", index_col="contig")
    # how='right' will only take contigs in the binning table.
    df = pd.merge(kmers, binning, left_index=True, right_index=True, how='right')
    df = pd.merge(taxonomy, df, left_index=True, right_index=True, how='right')
    df = pd.merge(coverage, df, left_index=True, right_index=True, how='right')
    df.write_csv("master.tsv", sep="\t", index="contig")

You can now run the following R scripts (preferably in RStudio) to examine your results.

Visualize bins
--------------

.. code-block:: R

    library("ggplot2")
    fpath="master.tsv"
    data = read.table(fpath, header=TRUE, sep='\t')
    ggplot(data, aes(x=x, y=y, color=cluster, group=cluster)) + geom_point(size=(sqrt(data$length))/100, shape=20, alpha=0.5) + theme_classic() + xlab('BH-tSNE X') + ylab('BH-tSNE Y') 

.. todo::
    Add an outoput image and explain the results

.. code-block:: R

    library("ggplot2")
    fpath="master.tsv"
    ggplot(data, aes(x=x, y=y, color=phylum, group=phylum)) + geom_point(size=(sqrt(data$length))/100, shape=20, alpha=0.5) + theme_classic() + xlab('BH-tSNE X') + ylab('BH-tSNE Y') 

.. todo::
    Add an outoput image and explain the results
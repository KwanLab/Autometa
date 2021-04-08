=================
Examining Results
=================

Automappa
=========

An interactive interface for exploration and refinement of metagenomes Automappa is a tool built to interface with Autometa output to help you explore your binning results.

For details, see the `Automappa page <https://github.com/WiscEvan/Automappa>`__

Visualize bins
==============

You can now run the following R scripts (preferably in RStudio) to examine your results.

.. code-block:: R

    # Load packages
    library("ggplot2")

    # Read the main binning table
    fpath="/Users/sidd/Research/autometa_v2/78mbp_metagenome.main.tsv"
    data = read.table(fpath, header=TRUE, sep='\t')
    
    # Fill empty cells as unclustered
    data$cluster <- sub("^$", "Unclustered", data$cluster)

    ggplot(data, aes(x=x_1, y=x_2, color=cluster, group=cluster)) +
        geom_point(size=(sqrt(data$length))/100, shape=20, alpha=0.5) + 
        theme_classic() + xlab('BH-tSNE X') + ylab('BH-tSNE Y') +
        guides( color = guide_legend( title = 'Cluster/bin' ))

.. image:: ../img/78Mbp_col_contigs.svg
    :alt: 78Mbp_col_contigs

In the above chart the points represent contigs. They are plotted on the two dimensions that result from dimension-reduction by BH-tSNE (or any other dimension reduction method that you have chosen) - you can think of the distance between points as being roughly proportional to their differences in 5-mer frequency. The points are also scaled in size according to the contig length, and they are colored by the assigned cluster/bin. You can see that there are some bins which are well-separated from others, but there are other bins that are close together. Cases like these might be worth investigating manually if you think, for instance, that multiple Autometa bins close together could actually be different parts of the same genome. If clusters are close together, there is also the possibility that contigs in the region have been misassigned.

In addition to using nucleotide composition, Autometa uses coverage and can also use taxonomy to distinguish contigs with similar composition. We can also visualize these differences with R.

.. code-block:: R

    ggplot(data, aes(x=x_1, y=x_2, color=phylum, group=phylum)) +
        geom_point(size=(sqrt(data$length))/100, shape=20, alpha=0.5) + 
        theme_classic() + xlab('BH-tSNE X') + ylab('BH-tSNE Y') +
        guides( color = guide_legend( title = 'Phylum' )) 

.. image:: ../img/78Mbp_col_phylum.svg
    :alt: 78Mbp_col_phylum

In the above plot, we have now colored the points by taxonomic phylum, and this reveals that several clusters that are close together in BH-tSNE space are in fact quite divergent from one another (like bottom left). This is probably the basis for Autometa's assignment of separate bins in these cases. 

In some cases, the contigs in a bin may in fact look divergent. You may want to manually examine cases such as these, but they could well be real if, for example, some contigs have few protein coding genes, or the organism is highly divergent from known sequences (see our paper `here <https://www.nature.com/articles/srep34362>`__ for some examples).

In this particular dataset, the coverages of all genomes are fairly similar, as revealed in the next plot:

.. code-block:: R

    ggplot(data, aes(x=coverage, y=gc_content, color=cluster, group=cluster)) +
        geom_point(size=(sqrt(data$length))/100, shape=20, alpha=0.5) +
        theme_classic() + xlab('BH-tSNE X') + ylab('BH-tSNE Y') +
        guides( color = guide_legend( title = 'Cluster/bin' )) 


.. image:: ../img/78Mbp_col_coverage.svg
    :alt: 78Mbp_col_coverage

In the above plot, the points are colored by cluster/bin again, and you can see that in this case, coverage is not much of a distinguishing feature. In other datasets, you may see closely related genomes at different coverages, which will be separatable by Autometa.
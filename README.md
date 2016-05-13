AutoMeta
========

An automated pipeline to deconvolute single metagenomic assemblies to separate individual bacterial and archaeal genomes. Contig separation is based on k-mer frequency 
analysis and [nonlinear dimension reduction](http://www.nature.com/articles/srep04516), an algorithm designed by Paul Wilmes' group. Contigs are automatically separated into 
groups using the [DBSCAN R package](https://cran.r-project.org/web/packages/dbscan/index.html).

Just tell me how to use it!
---------------------------

USAGE: 
```bash
autometa.py asm.fasta output_directory
```

Optionally, you can carry out the pre-processing step of separating contigs based on taxonomic kingdom (Eukaryota, Bacteria, Archaea and unclassified).  
```bash
separate_kingdoms.py asm.fasta output_directory 16
```
The last specifies the number of threads to use when running [Diamond](https://github.com/bbuchfink/diamond/) blast searches.

Dependencies/Requirements
-------------------------

Tested exclusively in a linux environment

* [DBSCAN R package](https://cran.r-project.org/web/packages/dbscan/index.html)
* (Optional, see above) [Diamond](https://github.com/bbuchfink/diamond/)
* [ggplot2 R package](http://ggplot2.org)
* [HMMER](http://hmmer.org)
* [MEGAN 6](http://ab.inf.uni-tuebingen.de/software/megan6/) 
* Perl
* Python 2.x
* Java >=1.8
* [R](https://www.r-project.org)



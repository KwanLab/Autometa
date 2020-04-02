Autometa
=========

An automated binning pipeline for single metagenomes, in particular host-associated and highly complex ones. Autometa is copyright 2020 Ian J. Miller, Evan Rees, Izaak Miller and Jason C. Kwan, and is released under the GNU Affero General Public License v3 (see LICENSE.txt). If you find Autometa useful to your work, please cite:

Miller, I. J.; Rees, E. R.; Ross, J.; Miller, I.; Baxa, J.; Lopera, J.; Kerby, R. L.; Rey, F. E.; Kwan, J. C. Autometa: Automated extraction of microbial genomes from individual shotgun metagenomes. *Nucleic Acids Research*, **2019**. [DOI: https://doi.org/10.1093/nar/gkz148](https://doi.org/10.1093/nar/gkz148)

Pipeline
---------

1. Filter by length cutoff
2. Count k-mer frequencies
3. Get coverage profiles
4. Assign taxonomies
5. Split by kingdom (bacteria or archaea)
6. Annotate kingdom with its respective markers
7. Get metagenome-assembled genomes (MAGs) from kingdom

For more details on the above process, please see our [paper](https://doi.org/10.1093/nar/gkz148).

Analysis of results
-------------------



Installation
------------

```bash
conda create -n autometa "python>=3.7" --yes
conda install -n autometa -c bioconda -c conda-forge --yes \
    biopython \
    pandas \
    tqdm \
    numpy \
    scikit-learn \
    scipy \
    samtools \
    bedtools \
    bowtie2 \
    hmmer \
    prodigal \
    diamond \
    ipython \
    ndcctools \
    parallel \
    requests \
    hdbscan \
    umap-learn \
    && conda clean --all --yes
```

Usage
-----

### Data Preparation

Before you run Autometa, you need to have assembled your shotgun metagenome. We
recommend using MetaSPAdes (part of the [SPAdes](http://cab.spbu.ru/software/spades/)
 package) after removing Illumina adaptors with
 [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

Note that if you use SPAdes or something else that names contigs like this: `NODE_1_length_319818_cov_8.03695` then Autometa can use the coverage information
 in the contig names. If you have used another assembler, then you first have to
 make a coverage table.

Fortunately, Autometa can construct this table for you with: `python -m autometa.commmon.coverage`

### (Quick Start) Metagenome Job Submission(s)

An Autometa metagenome job submission file is provided in the `tests` directory.
you may also find it [here](https://github.com/WiscEvan/Autometa/blob/dev/tests/metagenome.config) and below.

After you have filled out the job submission form, you may run the job via the
command:

`python autometa.py --metagenomes-configs </path/to/metagenome.config>`

Notice *multiple* configs can be passed to run binning runs. In the future, this
 will allow pan-genome analyses.

### (Advanced) Running modules

Many of the Autometa modules may be run standalone.

Simply pass in the `-m` flag when calling a script to signify to python you are
running an Autometa *module*.

I.e. `python -m autometa.common.kmers -h`

Database Dependencies
------------

Autometa comes packaged with the necessary markers files and will download and
format NCBI files during the first run.

#### NCBI:

* non-redundant protein database (nr) - [link](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz)
* taxdump.tar.gz - [link](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
    - required files within taxdump tarchive are *nodes.dmp*, *names.dmp* and *merged.dmp*
* prot.accession2taxid.gz - [link](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz)

#### Markers:

- bacteria single-copy-markers - [link](https://github.com/WiscEvan/Autometa/raw/dev/databases/markers/bacteria.single_copy.hmm)
- bacteria single-copy-markers cutoffs - [link](https://raw.githubusercontent.com/WiscEvan/Autometa/dev/databases/markers/bacteria.single_copy.cutoffs?token=AGF3KQVL3J4STDT4TJQVDBS6GG5FE)
- archaea single-copy-markers - [link](https://github.com/WiscEvan/Autometa/raw/dev/databases/markers/archaea.single_copy.hmm)
- archaea single-copy-markers cutoffs - [link](https://raw.githubusercontent.com/WiscEvan/Autometa/dev/databases/markers/archaea.single_copy.cutoffs?token=AGF3KQXVUDFIH6ECVTYMZQS6GG5KO)

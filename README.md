Autometa 2.0
=========

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

Database Dependencies
------------

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

Running Autometa
==================

### Running modules

Many of the Autometa modules may be run standalone.

Simply pass in the `-m` flag when calling a script to signify to python you are running an Autometa *module*.

I.e. `python -m autometa.common.kmers -h`

### Metagenome Job Submission(s):

An Autometa metagenome job submission file is provided in the `projects/tests` directory.
you may also find it [here](https://github.com/WiscEvan/Autometa/blob/dev/projects/tests/test_metagenome.config) and below.

After you have filled out the job submission form, you may run the job via the command:

`python autometa.py --metagenomes-configs </path/to/metagenome.config>`

Notice *multiple* configs can be passed to run binning runs. In the future, this will allow
pan-genome aware algorithms to use all metagenomes found within a project config.

#### <metagenome.config>
```shell
########################################
###### Metagenome Submission Files #####
########################################

#Place file path in respective location. If unavailable. You may leave as None.
[files]
metagenome = None
fwd_reads = None
rev_reads = None
length_filtered = None
coverages = None
kmer_counts = None
kmer_normalized = None
kmer_embedded = None
nucleotide_orfs = None
amino_acid_orfs = None
blastp = None
blastp_hits = None
blastx = None
taxonomy = None
bacteria_hmmscan = None
bacteria_markers = None
archaea_hmmscan = None
archaea_markers = None
binning = None

########################################
### Metagenome Submission Parameters ###
########################################

# See autometa.py --help for details
[parameters]
projects = None
project = 0
resume = 0
add_metagenome = None
length_cutoff = 10000
cov_from_spades = False
kmer_size = 5
kmer_multiprocess = True
kmer_normalize = True
do_pca = True
pca_dims = 50
embedding_method = UMAP
taxon_method = majority_vote
reversed = True
binning_method = recursive_dbscan
completeness = 20.0
purity = 90.0
verbose = False
force = False
usepickle = True
noparallel = True
cpus = 1
```

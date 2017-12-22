Autometa
========

An automated pipeline to deconvolute single metagenomic assemblies to separate individual bacterial and archaeal genomes. Autometa is free to use for academic and non-commercial users (see LICENSE.md). For commercial use contact [Jason Kwan](mailto:jason.kwan@wisc.edu). If you find Autometa useful to your work, please cite:

Miller, I. J.; Rees, E.; Ross, J.; Miller, I.; Baxa, J.; Lopera, J.; Kerby, R. L.; Rey, F. E.; Kwan, J. C. Autometa: Automated extraction of microbial genomes from individual shotgun metagenomes. *bioRxiv*, **2017**, xx

Dependencies
------------

Third party programs

* [Prodigal](https://github.com/hyattpd/prodigal/releases/) (Tested with v2.6.2)
* [HMMER](http://hmmer.org) (Tested with v3.1b2)
* [DIAMOND](https://github.com/bbuchfink/diamond) (Tested with v0.7.9.58)
* [Anaconda Python](https://www.anaconda.com)

Databases (these will be automatically downloaded the first time you run make\_taxonomy\_table.py)

* [nr](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz) (Note: tested with newer NR versions without GI numbers)
* [taxdump.tar.gz](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
* [prot.accession2taxid.gz](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz)

Python packages (Note: this list assumes the use of Anaconda Python)

* [tqdm](https://pypi.python.org/pypi/tqdm/4.19.5)
* [BioPython](http://biopython.org)
* [tsne](https://pypi.python.org/pypi/tsne)
* [joblib](https://pypi.python.org/pypi/joblib)

Installation
------------

The following was tested on a new install of Ubuntu Server 16.04.3 LTS.

First we install Prodigal, HMMER and DIAMOND.

```
sudo apt-get install prodigal hmmer build-essential zlib1g-dev
mkdir diamond
cd diamond
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.14/diamond-linux64.tar.gz
tar xvf diamond-linux64.tar.gz
```

At this point you should add the diamond directory to your $PATH environmental variable.

Install Anaconda Python.

```
wget https://repo.continuum.io/archive/Anaconda2-5.0.1-Linux-x86_64.sh
chmod +x Anaconda2-5.0.1-Linux-x86_64.sh
./Anaconda2-5.0.1-Linux-x86_64.sh
```

Install Python modules.

```
conda install tqdm joblib biopython
git clone https://github.com/danielfrg/tsne.git
sudo apt-get install libatlas-base-dev
cd tsne
python setup.py install
```

Now download Autometa.

```
git clone https://bitbucket.org/jason_c_kwan/autometa
cd autometa/pipeline
python setup_lca_functions.py build_ext --inplace
```

You should now add the autometa/pipeline directory to your $PATH environmental variable


Usage
-----
Before you run Autometa, you need to have assembled your shotgun metagenome. We recommend using MetaSPAdes (part of the [SPAdes](http://cab.spbu.ru/software/spades/) package) after removing Illumina adaptors with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

[Need to include here the procedure for getting coverage from a sam file]

Here we demonstrate the use of Autometa with the test dataset:

```
autometa/test_data/scaffolds.fasta
```

This dataset is the simulated dataset referred to as "78.125 Mbp" in our paper.

### Step 1: Split data into kingdom bins [optional]

We found that in host-associated metagenomes, this step vastly improves the binning performance of Autometa (and other pipelines) because less eukaryotic contigs will be binned into bacterial bins. However, if you are confident that you do not have eukaryotic contamination, or a mixture of Archaea and Bacteria, then you can skip this stage because it is rather computationally intensive. 

```
make_taxonomy_table.py -a ~/autometa/test_data/scaffolds.fasta -p 16 \
	-t ~/autometa/databases -l 3000
```

In the above command, we give the script make\_taxonomy\_table.py the assembly fasta file (-a), specify the number of CPUs to use (-p), the directory to use for database files (-t) and that we will just consider contigs above 3,000 bp (-l). The first time you run this script, it will automatically download the database files listed above, and format the nr database for DIAMOND to use. This script will then do the following:

1. Genes are identified in each contig with Prodigal.
2. Gene protein sequences are searched against nr with DIAMOND.
3. From the hits for each gene, the lowest common ancestor (LCA) is determined.
4. The taxonomy of each contig is determined by examining the LCA of each component protein (see paper for details)

#### Output files produced by make\_taxonomy\_table.py

File                         | Description 
-----------------------------|------------
Bacteria.fasta               | Contigs classified as bacterial  
scaffolds_filtered.fasta     | All contigs above the length cutoff 
scaffolds_filtered.fasta.tab | Table describing the GC content, length and coverage of filtered contigs
scaffolds_filtered.orfs.daa  | The output from DIAMOND (binary format)
scaffolds_filtered.orfs.faa  | ORF translations obtained from Prodigal
scaffolds_filtered.orfs.lca  | Table describing the lowest common ancestor (LCA) for each ORF
scaffolds_filtered.orfs.tab  | The output from DIAMOND (tab delimited format)
scaffolds_filtered.txt       | The full output from Prodigal
taxonomy.tab                 | Table with information from scaffolds_filtered.fasta.tab, and also the full assigned taxonomy for each contig
unclassified.fasta           | Contigs unclassified on the Kingdom level

Note that in our test data, there are no non-bacterial contigs. For other datasets, make\_taxonomy\_table.py will produce other fasta files, for example Eukaryota.fasta and Viruses.fasta.

### Step 2: Bin bacterial contigs with BH-tSNE and DBSCAN

Note: This procedure can also be used on archaeal sequences, but will most likely not work well on other kingdoms, because the relevant single-copy gene markers are not included in Autometa. 

Here we are running Autometa on the Bacteria.fasta file made in step 1.

```
run_autometa.py -a Bacteria.fasta -p 16 -l 3000 -t taxonomy.tab
```

In the above command, we are supplying Bacteria.fasta to Autometa, and also the taxonomy table (taxonomy.tab) produced in step 1. If we supply a taxonomy table, then this information is used to help with clustering. Otherwise, Autometa clusters solely on 5-mer frequency and coverage. In the above command we are also specifying to use 16 CPUs (-p), and that our length cutoff is 3,000 bp (-l). We are using the default output directory of the current working directory (this can be set with the -o flag), and by default the pipeline assumes we are looking at bacterial contigs (use -k archaea otherwise). The script will do the following:

1. Find single-copy marker genes in the input contigs with HMMER
2. Generate two-dimensional BH-tSNE coordinates for each contig based on 5-mer frequencies
3. Cluster contigs based on BH-tSNE coordinates, coverage and (optionally) taxonomy
4. Accept clusters that are estimated to be over 20% complete and 90% pure based on single-copy marker genes
5. Unclustered contigs leftover will be re-clustered until no more acceptable clusters are yielded

For more details on the above process, please see our paper. If you include a taxonomy table in the run\_autometa.py command, Autometa will attempt to further partition the data based on ascending taxonomic specificity (i.e. in the order phylum, class, order, family, genus, species) when clustering unclustered contigs from a previous attempt. We found that this is mainly useful if you have a highly complex metagenome (lots of species), or you have several related species at similar coverage level.

#### Output files produced by run\_autometa.py

File | Description
-----|------------
Bacteria\_filtered.hmm.tbl | Output from HMMER
Bacteria\_filtered\_marker.tab | Table describing the marker genes found in each contig
k-mer\_matrix | Raw 5-mer frequencies for each contig
recursive\_dbscan\_output.tab | Output table containing the cluster (bin) for each contig


### Step 3: Recruit unclustered contigs to bins through supervised machine learning [optional]

In this step we use supervised machine learning to classify the unclustered contigs to the bins that we have produced (formally, the bins produced in step 2 are the training set). Depending on the size of your dataset, this step can be computationally intensive. 

```
ML_recruitment.py -t recursive_dbscan_output.tab -r -m k-mer_matrix -o ML_recruitment_output.tab
```

In the above command, we give ML\_recruitment.py the output table from step 2 (recursive\_dbscan\_output.tab, -t), as well as the k-mer\_matrix file produced in step 2 (-m), and specify the output file (ML\_recruitment\_output.tab, -o). We also use the -r flag, specifying that the script will run recursively, adding contigs it classifies to the training set and re-classifying until 0 more classifications are yielded. By default, classifications are only made if 10 out of 10 repeat classifications agree, and if the classification would not increase the apparent contamination estimated by the presence of single-copy marker genes. The specified output file is a table with the following columns:



### Running all steps in sequence

For convenience, it is possible to run all three of the above steps through run\_autometa.py, as shown below:

```
[Fill in later]
```


=============
metagenome.py
=============

.. code-block:: shell
 
	usage: temp_write.py [-h] --ncbi NCBI [--cutoff CUTOFF]
	                     [--kmer-size KMER_SIZE]
	                     [--kmer-normalized KMER_NORMALIZED]
	                     [--kmer-fpath KMER_FPATH]
	                     [--taxon-method {majority_vote}]
	                     [--vote-fname VOTE_FNAME] [--cpus CPUS] [--noparallel]
	                     [--force] [--verbose]
	                     assembly

	Metagenome class holding attributes and methods to manipulate metagenome
	assemblies.

	positional arguments:
	  assembly              </path/to/assembly.fasta>

	optional arguments:
	  -h, --help            show this help message and exit
	  --ncbi NCBI           </path/to/ncbi/dir>
	  --cutoff CUTOFF       length to filter sequences
	  --kmer-size KMER_SIZE
	  --kmer-normalized KMER_NORMALIZED
	                        Perform CLR transform on k-mer frequencies if
	                        provided. (</path/to/kmers.normalized.tsv>)
	  --kmer-fpath KMER_FPATH
	                        </path/to/kmers.tsv>
	  --taxon-method {majority_vote}
	  --vote-fname VOTE_FNAME
	                        <taxon-vote.tsv filename>
	  --cpus CPUS           num cpus to use
	  --noparallel          Do not use GNU parallel
	  --force               overwrite existing files
	  --verbose             add verbosity

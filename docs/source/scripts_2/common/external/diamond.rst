==========
diamond.py
==========

.. code-block:: shell
 
	usage: temp_write.py [-h] [--evalue EVALUE] [--maxtargetseqs MAXTARGETSEQS]
	                     [--cpus CPUS] [--tmpdir TMPDIR] [--top-pct TOP_PCT]
	                     [--force] [--verbose]
	                     fasta database acc2taxids outfile {blastp,blastx}

	Retrieves blastp hits with provided input assembly

	positional arguments:
	  fasta                 </path/to/faa/file>
	  database              </path/to/diamond/formatted/database>
	  acc2taxids            </path/to/ncbi/prot.accession2taxid.gz>
	  outfile               </path/to/diamond/output/file>
	  {blastp,blastx}       [blastp]: A.A -> A.A. [blastx]: Nucl. -> A.A.

	optional arguments:
	  -h, --help            show this help message and exit
	  --evalue EVALUE       diamond evalue threshold
	  --maxtargetseqs MAXTARGETSEQS
	                        max target sequences to retrieve per query
	  --cpus CPUS           num cpus to use
	  --tmpdir TMPDIR       </path/to/tmp/directory>
	  --top-pct TOP_PCT     top percentage of hits to retrieve
	  --force               force overwrite of diamond output table
	  --verbose             add verbosity

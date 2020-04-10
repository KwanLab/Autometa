================
majority_vote.py
================

.. code-block:: shell
 
	usage: modified majority vote [-h] [--dbdir DBDIR] [--lca-fname LCA_FNAME]
	                              [--blast-table BLAST_TABLE] [--nopickle]
	                              [--verbose]
	                              fasta outdir outfname

	positional arguments:
	  fasta                 </path/to/prot/orfs.faa>
	  outdir                <path/to/output/dir>
	  outfname              <output filename>

	optional arguments:
	  -h, --help            show this help message and exit
	  --dbdir DBDIR         <path/to/ncbi/dir>
	  --lca-fname LCA_FNAME
	                        <lca filename>
	  --blast-table BLAST_TABLE
	                        </path/to/blastp.tsv>
	  --nopickle            do not pickle objects to disk
	  --verbose             add verbosity

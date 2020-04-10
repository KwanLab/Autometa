======
lca.py
======

.. code-block:: shell
 
	usage: Script to determine lowest common ancestor [-h]
	                                                  [--blast-hits BLAST_HITS]
	                                                  [--dbdir DBDIR] [--nopickle]
	                                                  [--verbose] [--force]
	                                                  orfs blast outdir outfname

	positional arguments:
	  orfs                  <path/to/orfs.fasta>
	  blast                 <path/to/blast.tsv>
	  outdir                <path/to/output/dir>
	  outfname              <lca filename>

	optional arguments:
	  -h, --help            show this help message and exit
	  --blast-hits BLAST_HITS
	                        <path/to/blast.pkl.gz> (with taxids already added)
	  --dbdir DBDIR         <path/to/ncbi/dir>
	  --nopickle            do not pickle objects to disk
	  --verbose             add verbosity
	  --force               force overwrite if file already exists

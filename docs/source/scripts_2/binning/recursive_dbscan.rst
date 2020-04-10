===================
recursive_dbscan.py
===================

.. code-block:: shell
 
	usage: temp_write.py [-h] [--embedded-kmers EMBEDDED_KMERS]
	                     [--embedding-method {TSNE,UMAP}]
	                     [--clustering-method {DBSCAN,HDBSCAN}]
	                     [--completeness COMPLETENESS] [--purity PURITY]
	                     [--taxonomy TAXONOMY] [--domain {bacteria,archaea}]
	                     [--verbose]
	                     kmers coverage markers out

	Perform decomposition/embedding/clustering via PCA/[TSNE|UMAP]/DBSCAN.

	positional arguments:
	  kmers                 </path/to/kmers.tsv>
	  coverage              </path/to/coverages.tsv>
	  markers               </path/to/markers.tsv>
	  out                   </path/to/output.tsv>

	optional arguments:
	  -h, --help            show this help message and exit
	  --embedded-kmers EMBEDDED_KMERS
	                        </path/to/embedded_kmers.tsv>
	  --embedding-method {TSNE,UMAP}
	                        Embedding method to use
	  --clustering-method {DBSCAN,HDBSCAN}
	                        Clustering method to use
	  --completeness COMPLETENESS
	                        <completeness cutoff>
	  --purity PURITY       <purity cutoff>
	  --taxonomy TAXONOMY   </path/to/taxonomy.tsv>
	  --domain {bacteria,archaea}
	                        Kingdom to consider (archaea|bacteria)
	  --verbose             log debug information

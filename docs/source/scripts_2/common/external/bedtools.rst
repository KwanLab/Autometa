===========
bedtools.py
===========

.. code-block:: shell
 
	usage: temp_write.py [-h] [--coverage COVERAGE] [--force-bed] [--force-cov]
	                     ibam lengths bed

	positional arguments:
	  ibam                 </path/to/BAM/alignment.bam>
	  lengths              </path/to/genome/lengths.tsv> tab-delimited
	                       cols=[contig,length]
	  bed                  </path/to/alignment.bed> tab-delimited
	                       cols=[contig,length]

	optional arguments:
	  -h, --help           show this help message and exit
	  --coverage COVERAGE  </path/to/coverage.tsv>
	  --force-bed          force overwrite `bed`
	  --force-cov          force overwrite `--coverage`

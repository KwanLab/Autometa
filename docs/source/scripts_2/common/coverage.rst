===========
coverage.py
===========

.. code-block:: shell
 
	usage: coverage.py

	Construct contig coverage table given an input assembly and reads.

	optional arguments:
	  -h, --help            show this help message and exit
	  -f ASSEMBLY, --assembly ASSEMBLY
	                        </path/to/metagenome.fasta>
	  -1 [FWD_READS [FWD_READS ...]], --fwd-reads [FWD_READS [FWD_READS ...]]
	                        </path/to/forwards-reads.fastq>
	  -2 [REV_READS [REV_READS ...]], --rev-reads [REV_READS [REV_READS ...]]
	                        </path/to/reverse-reads.fastq>
	  -U [SE_READS [SE_READS ...]], --se-reads [SE_READS [SE_READS ...]]
	                        </path/to/single-end-reads.fastq>
	  --sam SAM             </path/to/alignments.sam>
	  --bam BAM             </path/to/alignments.bam>
	  --lengths LENGTHS     </path/to/lengths.tsv>
	  --bed BED             </path/to/alignments.bed>
	  --nproc NPROC         Num processors to use. (default: 4)
	  --from-spades         Extract k-mer coverages from contig IDs. (Input
	                        assembly is output from SPAdes)
	  --out OUT             </path/to/coverages.tsv>

=========
bowtie.py
=========

.. code-block:: shell
 
	usage: temp_write.py [-h] [-1 [FWD_READS [FWD_READS ...]]]
	                     [-2 [REV_READS [REV_READS ...]]]
	                     [-U [SE_READS [SE_READS ...]]] [--nproc NPROC]
	                     assembly database sam

	Align provided reads to metagenome `assembly` and write alignments to `sam`.
	NOTE: At least one reads file is required.

	positional arguments:
	  assembly              </path/to/assembly.fasta>
	  database              </path/to/alignment.database>. Will construct database
	                        at provided path if not found.
	  sam                   </path/to/alignment.sam>

	optional arguments:
	  -h, --help            show this help message and exit
	  -1 [FWD_READS [FWD_READS ...]], --fwd-reads [FWD_READS [FWD_READS ...]]
	                        </path/to/forward-reads.fastq>
	  -2 [REV_READS [REV_READS ...]], --rev-reads [REV_READS [REV_READS ...]]
	                        </path/to/reverse-reads.fastq>
	  -U [SE_READS [SE_READS ...]], --se-reads [SE_READS [SE_READS ...]]
	                        </path/to/single-end-reads.fastq>
	  --nproc NPROC         Num processors to use.

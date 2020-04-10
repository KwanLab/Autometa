========
hmmer.py
========

.. code-block:: shell
 
	usage: Retrieves markers with provided input assembly [-h] [--log LOG]
	                                                      [--force] [--cpus CPUS]
	                                                      [--parallel] [--verbose]
	                                                      orfs hmmdb cutoffs
	                                                      hmmscan markers

	positional arguments:
	  orfs         </path/to/assembly.orfs.faa>
	  hmmdb        </path/to/hmmpressed/hmmdb>
	  cutoffs      </path/to/hmm/cutoffs.tsv>
	  hmmscan      </path/to/hmmscan.out>
	  markers      </path/to/markers.tsv>

	optional arguments:
	  -h, --help   show this help message and exit
	  --log LOG    </path/to/parallel.log>
	  --force      force overwrite of out filepath
	  --cpus CPUS  num cpus to use
	  --parallel   Enable GNU parallel
	  --verbose    add verbosity

===========
prodigal.py
===========

.. code-block:: shell
 
	usage: prodigal.py

	Calls ORFs with provided input assembly

	positional arguments:
	  assembly     </path/to/assembly>
	  nucls_out    </path/to/nucls.out>
	  prots_out    </path/to/prots.out>

	optional arguments:
	  -h, --help   show this help message and exit
	  --force      force overwrite of ORFs out filepaths
	  --cpus CPUS  num cpus to use
	  --parallel   Enable GNU parallel
	  --verbose    add verbosity

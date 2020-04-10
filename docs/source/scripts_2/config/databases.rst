============
databases.py
============

.. code-block:: shell
 
	usage: databases config [-h] [--config CONFIG] [--dryrun] [--nproc NPROC]
	                        [--out OUT]

	optional arguments:
	  -h, --help       show this help message and exit
	  --config CONFIG  </path/to/input/database.config>
	  --dryrun         Log configuration actions but do not perform them.
	  --nproc NPROC    num. cpus to use for DB formatting. (default 4)
	  --out OUT        </path/to/output/database.config>

	By default, with no arguments, will download/format databases into default
	databases directory.

=============
work_queue.py
=============

.. code-block:: shell
 
	usage: temp_write.py [-h] [--stats STATS] [--transactions TRANSACTIONS]
	                     [--debug DEBUG]
	                     password project-name

	Uses WorkQueue as a task-manager to run different tasks within the Autometa
	pipeline on a scalable computing system.

	positional arguments:
	  password              </path/to/work-queue/master/password.txt>
	  project-name          <unique project identifier>

	optional arguments:
	  -h, --help            show this help message and exit
	  --stats STATS         </path/to/workqueue/stats.log
	  --transactions TRANSACTIONS
	                        </path/to/workqueue/transactions.log
	  --debug DEBUG         </path/to/workqueue/debug.log

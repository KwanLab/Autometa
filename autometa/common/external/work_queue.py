#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT

Work Queue API to scale Autometa job submissions allowing for dynamically
scalable and robust application across multiple machines, clusters, or
computing grids.

"""


import logging
import os
import sys

import work_queue as wq

logger = logging.getLogger(__name__)


def start_master(port, name, password, stats=None, debug=None, transactions=None):
    """Start a work-queue master to handle task submissions.

    Notes
    -----

    * Reference: https://cctools.readthedocs.io/en/latest/work_queue/

    Parameters
    ----------
    port : int
        port for master to operate on
    name : str
        Name of work-queue master to reference with workers
    password : str
        </path/to/file/password.txt>
        Following recommended security practices, a password allows only specific
        workers to be able to connect to work-queue master.
    debug : str, optional
        </path/to/output/debug.log>
    stats : str, optional
        </path/to/output/stats.log>

    Returns
    -------
    work_queue.WorkQueue
        work-queue master awaiting task submissions

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    q = wq.WorkQueue(
        port=port,
        name=name,
        stats_log=stats,
        transactions=transactions,
        debug_log=debug)
    q.specify_password_file(password)
    logger.debug(f'({q.name}) wq-master started')
    return q

def main(args):

    q = start_master(
        port=args.port,
        name=args.project_name,
        password=args.password,
        debug_log=args.debug,
        stats_log=args.stats)

    # Create a task for the queue
    # task = wq.Task("exe_name infname args outfname")

    # Specify resource requirements: (Default will consume entire worker)
    # t.specify_cores(required_cores)    # Integer
    # t.specify_memory(required_memory) # Integer specified in MB
    # t.specify_disk(required_disk)  # Integer specified in MB
    raise NotImplementedError
    # Add required executables
    for exe_fp,exe_name in executables:
        task.specify_input_file(exe_fp, exe_name, cache=True)
    # Add required input filepaths
    for infpath,infname in infpaths:
        task.specify_input_file(infpath, infname, cache=False)
    # Add required output filepaths
    for outfpath,outfname in outfpaths:
        task.specify_output_file(outfpath, outfname, cache=False)
    # Submit task to the queue
    task_id = q.submit(task)

    # Retrieve results from the queue
    while not q.empty():
        t = q.wait(5)
        if t:
            print("Task {} has returned!".format(t.id))

            if t.return_status == 0:
                print("command exit code:\n{}".format(t.exit_code))
                print("stdout:\n{}".format(t.output))
            else:
                print("There was a problem executing the task.")

    # NOTE: The size of output is limited to 1 GB so the output should be file paths!
    # ... Not the contents of the file.

    # Pipelined Submission.
    # If you have a very large number of tasks to run, it may not be possible to submit all
    # of the tasks, and then wait for all of them. Instead, submit a small number of tasks,
    # then alternate waiting and submiting to keep a constant number in the queue. The
    # hungry will tell you if more submissions are warranted:

    if q.hungry():
        # submit more tasks...
        pass

    # Watching Output Files:

    task.specify_output_file("my-file", flags = wq.WORK_QUEUE_WATCH)

    # Kill workers that are executing tasks twice as slow as compared to the
    # average.
    q.activate_fast_abort(2)

    # create task as usual and tag it with an arbitrary string.
    t = wq.Task(...)
    t.specify_tag("my-tag")

    taskid = q.submit(t)

    # cancel task by id. Return the canceled task.
    t = q.cancel_by_taskid(taskid)

    # or cancel task by tag. Return the canceled task.
    t = q.cancel_by_tasktag("my-tag")


    # if task fails given a worker misconfiguration:
    q.blacklist(t.hostname)

    # Retrieve all stats from the queue
    stats = q.stats
    print(stats.workers_busy)

# Main program
if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='[%(asctime)s %(levelname)s] %(name)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)

    parser = argparse.ArgumentParser(description="""
    Uses WorkQueue as a task-manager to run different tasks within the Autometa
    pipeline on a scalable computing system.
    """)
    parser.add_argument('password', help='</path/to/work-queue/master/password.txt>')
    parser.add_argument('project-name', help='<unique project identifier>')
    # TODO: Allow argparse to handle range of ints for port
    parser.add_argument('--port', help="""
    The port number to listen on. If zero, then a random port is chosen. A range
     of possible ports (low, hight) can be also specified instead of a single integer.
    """,
        default=wq.WORK_QUEUE_DEFAULT_PORT, type=int)
    parser.add_argument('--stats', help='</path/to/workqueue/stats.log')
    parser.add_argument('--transactions', help='</path/to/workqueue/transactions.log')
    parser.add_argument('--debug', help='</path/to/workqueue/debug.log')

    args = parser.parse_args()
    main(args)

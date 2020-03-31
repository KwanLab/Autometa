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

Uses Work Queue as a task-manager to run different tasks within the autometa pipeline on a scalable computing system.
"""

import work_queue as wq

import os
import sys


def init_queue(name='autometa', password, port=0, debug_log=None, stats_log=None):
    # create a new queue
    # We create the tasks queue using the default port. If this port is already
    # been used by another program, you can try setting port = 0 to use an
    # available port.
    try:
        q = wq.WorkQueue(
            name=name,
            port=port,
            debug_log=debug_log,
            stats_log=stats_log)
    except:
        print('Instantiation of Work Queue failed!')
        sys.exit(1)

    q.specify_password_file(password)
    return q

def check_executable(executable):
    executable_path = os.path.realpath(executable)
    if not os.path.exists(executable_path):
        executable_path = '/usr/bin/sleep'
        if not os.path.exists(executable_path):
            print(f'{executable_path} was not found. Please modify the executable_path')
            sys.exit(1)
    return executable_path

def main(args):
    q = init_queue(args.project_name, port=args.port, args.debug, args.stats)
    # We have to specify precisely which files need to be transmitted to the workers.
    # We record the location of executable in 'executable_path'

    print(f'listening on port {q.port}...')
    tasks = {
        'run_autometa':'autometa.py',
        'call_orfs':'prodigal.py',
        'get_markers':'hmmer.py',
    }
    # We create and dispatch a task for each filename given in the argument list
    for task,executable in tasks.items():
        executable_path = check_executable(executable)

        # Note that we write ./gzip here, to guarantee that the gzip version we
        # are using is the one being sent to the workers.
        command = "./gzip < %s > %s" % (infile, outfile)

        t = Task(command)

        # gzip is the same across all tasks, so we can cache it in the workers.
        # Note that when specifying a file, we have to name its local name
        # (e.g. executable_path), and its remote name (e.g. "gzip"). Unlike the
        # following line, more often than not these are the same.
        t.specify_file(executable_path, executable, wq.WORK_QUEUE_INPUT, cache=True)
        t.specify_cores(2)    #needs 2 cores
        t.specify_memory(100) #needs 100 MB memory
        t.specify_disk(1000)  #needs 1 GB disk
        # files to be compressed are different across all tasks, so we do not
        # cache them. This is, of course, application specific. Sometimes you may
        # want to cache an output file if is the input of a later task.
        t.specify_file(infile, infile, wq.WORK_QUEUE_INPUT, cache=False)
        t.specify_file(outfile, outfile, wq.WORK_QUEUE_OUTPUT, cache=False)

        # Once all files has been specified, we are ready to submit the task to the queue.
        taskid = q.submit(t)
        print("submitted task (id# %d): %s" % (taskid, t.command))
    stats = q.stats
    print(stats.workers_busy)
    print("waiting for tasks to complete...")
    while not q.empty():
        t = q.wait(5)
        if t:
            print("task (id# %d) complete: %s (return code %d)" % (t.id, t.command, t.return_status))
            if t.return_status != 0:
              # The task failed. Error handling (e.g., resubmit with new parameters, examine logs, etc.) here
              None
        #task object will be garbage collected by Python automatically when it goes out of scope

    print("all tasks complete!")

    #work queue object will be garbage collected by Python automatically when it goes out of scope
    sys.exit(0)

# Main program
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('password', help='</path/to/work-queue/master/password/file>')
    parser.add_argument('--port', default=wq.WORK_QUEUE_DEFAULT_PORT)
    parser.add_argument('--project-name', default='myautometaproject')
    parser.add_argument('--debug', default='my.debug.log')
    parser.add_argument('--stats', default='my.stats.log')
    args = parser.parse_args()
    main(args)

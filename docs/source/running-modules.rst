Running Autometa
================

Running modules
---------------

Many of the Autometa modules may be run standalone.

Simply pass in the ``-m`` flag when calling a script to signify to python you are running an Autometa *module*.

I.e. ``python -m autometa.common.kmers -h``

Metagenome Job Submission(s)
----------------------------

An Autometa metagenome job submission file is provided in the `projects/tests` directory.
you may also find it [here](https://github.com/WiscEvan/Autometa/blob/dev/projects/tests/test_metagenome.config) and below.

After you have filled out the job submission form, you may run the job via the command:

``python autometa.py --metagenomes-configs </path/to/metagenome.config>``

Notice *multiple* configs can be passed to run binning runs. In the future, this will allow
pan-genome aware algorithms to use all metagenomes found within a project config.

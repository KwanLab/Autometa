================
Running Autometa
================

Before running anything make sure you

1. Running python3
2. Have activated the conda environment using:
    ``conda activate autometa``
3. Are in the Autometa directory
4. If you are running the main autometa.py make sure that the Autometa directory is in your $PATH environmental variable

Running modules
===============

Many of the Autometa modules may be run standalone.

Simply pass in the ``-m`` flag when calling a script to signify to python you are 
running an Autometa *module*.

I.e. ``python -m autometa.common.kmers -h``

Running functions
=================

Many of the Autometa functions may be run standalone as well. It is same as importing any other python
function.

.. code-block:: python
    :linenos:

    python3
    import autometa.common.external.samtools
    autometa.common.common.external.samtools.sort(sam = <path/to/sam/file>, out=<path/to/output/file>, npoc = 4)


Data preparation
================

Before you run Autometa, you need to have assembled your shotgun metagenome. 
We recommend using MetaSPAdes (part of the SPAdes_ package) after removing Illumina 
adaptors with Trimmomatic_.

Note that if you use SPAdes or something else that names contigs like 
this: ``NODE_1_length_319818_cov_8.03695`` then Autometa can use the coverage 
information in the contig names. If you have used another assembler, then 
you first have to make a coverage table.

Fortunately, Autometa can construct this table for you with: ``python -m autometa.commmon.coverage``


Metagenome Job Submission(s)
============================

An Autometa metagenome job submission file is provided in the `projects/tests` directory.
You may also find it `here <https://github.com/WiscEvan/Autometa/blob/dev/tests/metagenome.config>`_.

After you have filled out the job submission form, you may run the job via the command:

``python autometa.py --metagenomes-configs </path/to/metagenome.config>``

Notice *multiple* configs can be passed to run binning runs. In the future, this will allow
pan-genome aware algorithms to use all metagenomes found within a project config.

.. _SPAdes: http://cab.spbu.ru/software/spades/
.. _Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

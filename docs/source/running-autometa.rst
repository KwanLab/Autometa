================
Running Autometa
================


Before running anything make sure you have activated the conda environment using
``conda activate autometa`` and are in the Autometa directory (``cd </path/to/Autometa>``).

See the install page for details on setting up your conda environment.


Testing Autometa
================

.. code-block:: shell

    # Check dependencies to see if Autometa environment is appropriately configured
    python autometa.py --check-dependencies

    # If any of the checks return False, you can check which failed using
    python autometa.py --check-dependencies --debug


Once you have resolved the dependencies and your checks are returning ``True``, you
can move on to running the pipeline.

Data preparation
================

Before you run Autometa, you need to have assembled your shotgun metagenome.
We recommend using MetaSPAdes (part of the SPAdes_ package) after removing Illumina
adaptors with Trimmomatic_.

Note that if you use SPAdes or something else that names contigs like
this: ``NODE_235_length_161117_cov_116.709`` then Autometa can use the coverage
information in the contig names. If you have used another assembler, then
you first have to make a coverage table.

Fortunately, Autometa can construct this table for you with: ``python -m autometa.commmon.coverage``

Metagenome Job Submission(s)
============================

An Autometa metagenome job submission file is provided in the `tests` directory.
You may find it here: `Autometa/tests/metagenome.config <https://github.com/KwanLab/Autometa/blob/dev/tests/metagenome.config>`_.

If you are not providing a file, under the ``[files]`` section, leave the value as None,
and Autometa will write update this value with the file it will generate.

The required value(s) for running the full pipeline are ``metagenome`` and at least one file
for calculating contig coverage. Note, if you have assembled your metagenome
with SPAdes, you can specify ``cov_from_spades = True`` in the ``[parameters]``
section and *only* your metagenome assembly is required.

Otherwise, you will need to provide *at least one* file in the metagenome.config ``[files]`` section:

.. code-block:: python

    [files]
    ...
    # Multiple Reads files of respective format may be provided using a comma-delimiter
    # Paired-end forward reads
    fwd_reads = None
    # Paired-end reverse reads
    rev_reads = None
    # Single-end reads
    se_reads = None
    # Reads aligned to metagenome assembly (in sam format)
    sam = None
    # Reads sorted from alignment to metagenome assembly (in bam format)
    bam = None
    # bedtools genomecov output file (in bed format)
    bed = None
    # Tab-delimited table of (at least) contig and coverage columns
    coverages = None
    ...

After you have filled out the job submission form, you may run the job via the command:

``python autometa.py </path/to/metagenome.config>``

Notice *multiple* configs can be passed to run binning runs. In the future, this will allow
pan-genome aware algorithms to use all metagenomes found within a project config.

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

    from autometa.common.external import samtools

    samtools.sort(sam=<path/to/sam/file>, out=<path/to/output/file>, nproc=4)


.. _SPAdes: http://cab.spbu.ru/software/spades/
.. _Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

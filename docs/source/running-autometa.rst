================
Running Autometa
================

Data preparation
================

Before you run Autometa, you need to have assembled your shotgun metagenome. The following workflow is recommended:

#. Trim adapter sequences from the reads. We prefer to use Trimmomatic_, but you can go ahead and use any tool of your preference.
#. Quality check of reads to make sure that the adapters have been removed, we use FastQC_ for this.
#. Assemble the trimmed reads. We recommend using MetaSPAdes which is a part of the SPAdes_ package to assemble the trimmed reads but you can use any other assembler as well.
#. An optional thing to do here would be to check the quality of your assembly as well. This would give you a variety of assembly statistics one of which is N50 which will be useful in selecting the cutoff value during the Autometa length-filter step. We tend to use metaQuast_ for this (use ``--min-contig 1`` option to get an accurate N50).

.. note::

    If you use SPAdes then Autometa can use the k-mer coverage information in the contig names. If you have used any other assembler, then you first have to make a coverage table.

    Fortunately, Autometa can construct this table for you with: ``autometa-coverage``. Use ``--help`` to get the complete usage.

Nextflow walkthrough
====================

Why nextflow
------------

Nextflow helps Autometa produce reproducible results while allowing the pipeline to scale across different platforms and hardware.


System Requirements
-------------------

Nextflow 

Currently the nextflow pipeline only works with Docker so Docker must be installed on your system `Get Docker <https://docs.docker.com/get-docker>`_. We do plan on removing this Docker dependency.

Nextflow runs on any Linux compatible system or MacOS with Java installed. 


Quick Start
------------


Installation
^^^^^^^^^^^^
Using `conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ install nf-core and nextflow into an environment (here called 'autometa-nf')

.. code-block:: bash

    conda create -c conda-forge -c bioconda --name autometa-nf python=3 nf-core nextflow -y


Once it finishes installing be sure to active the environment:

.. code-block:: bash

    conda activate autometa-nf

Lastly, download the pipeline from GitHub using nextflow

.. code-block:: bash

    nextflow pull https://github.com/KwanLab/Autometa -r main



Launching Autometa
^^^^^^^^^^^^^^^^^^

Run the pipeline using the command below. 

Note: Unless specified using the parameters (:code:`--interim_dir`, :code:`--outdir`, and :code:`--tracedir`), intermediate and temporary files will be created relative to the directory this is executed from. 

.. code-block:: bash

    nf-core launch KwanLab/Autometa

You will then be provided two "launch method" options: "Web based" or "Command line". While it is possible to use the command line version, it is preferred and easier to use the Web based GUI.
Use the arrow keys to select one or the other and then press return/enter.


Note: You can use `tmux <https://github.com/tmux/tmux/wiki>`_ or `screen <https://www.gnu.org/software/screen/>`_ in case you want to exit the window or disconnect from the server.

Set autometa parameters with web based GUI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The GUI will present all available parameters, though some extra parameters may be hidden. Toggle "Show hidden params" on the right side of the page to reveal hidden parameters.

The only mandatory parameters are :code:`--input` which is the path to your input metagenome's nucleotide FASTA file, and :code:`-profile`.

Currently the availble options for :code:`-profile` are 

* **standard**: runs all process jobs locally. If you use slurm, don't use this.
* **basic_slurm**: submits all process jobs into the slurm queue. See :ref:`using-slurm:` before using
* **docker**: is currently required

An example input for locally-executed jobs would be 

:code:`-profile`: standard,docker

An example input for slurm-executed jobs would be 

:code:`-profile`: basic_slurm,docker


Running the pipeline
^^^^^^^^^^^^^^^^^^^^

After you are finished double-checking your parameter settings, click "Launch" at the top right of web based GUI page, or "Launch workflow" at the bottom of the page.




Advanced Nextflow
-----------------


Multiple Inputs
^^^^^^^^^^^^^^^

You can also input multiple asseblies at once with the help of wildcards. In the below example all the files with extension "fasta" would be taken as input by nextflow_.

.. code-block:: bash

    // Find this section of code in parameters.config
    params.metagenome = "$HOME/tutorial/test_data/*.fna"
    params.interim = "$HOME/tutorial/interim/" 
    params.processed = "$HOME/tutorial/processed/"

.. note::
    1. Wildcard characters will only be interpreted when "double quotes" are used
    2. If the interim and processed directories doens't exist, nextflow would generate them 

Database directory
^^^^^^^^^^^^^^^^^^

.. todo::

Autometa uses the following NCBI databses throughout its pipeline:

- Non-redundant `nr database <ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz>`_
- `prot.accession2taxid.gz <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_
- *nodes.dmp*, *names.dmp* and *merged.dmp* from `taxdump tarball <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_ 

If you are running autometa for the first time you'll have to download these databases. Use ``autometa-update-databases --update-ncbi``. This will download the databases to the default path. You can check the default paths using ``autometa-config --print``. If you need to change the default download directory you can use ``autometa-config --section databases --option ncbi --value <path/to/new/ncbi_database_directory>``. See ``autometa-update-databases -h`` and ``autometa-config-h`` for full list of options.

In your ``parameters.config`` file you also need to specify the directory where the different databases are present. Make sure that the directory path contains the following databases:

- Diamond formatted nr file => nr.dmnd
- Extracted files from tarball taxdump.tar.gz
- prot.accession2taxid.gz

.. code-block:: bash

    // Find this section of code in parameters.config
    // Update this path to folder with all NCBI databases
    params.single_db_dir = "/Autometa/autometa/databases/ncbi"

CPUs, Memory, Disk
^^^^^^^^^^^^^^^^^^

Like nf-core pipelines, we have set some automatic defaults for Autometa's processes. These are dynamic and each process will try a second attempt using more resources if the first fails due to resources. Resources are always capped by the parameters (show with defaults):
 - :code:`-max_cpus = 2` 
 - :code:`-max_memory = 6.GB`
 - :code:`-max_time = 48.h`

The best practice to change the resources is to create a new config file and point to it at runtime by adding the flag :code:`-c path/to/config_file`


For example, to give all resource-intensive jobs more memory, create a file called :code:`overwrite_config.config` and insert

.. code-block:: bash
    
    process {
      withLabel:process_high {
        memory = 200.GB
      }
    }

Then your command to run the pipeline (assuming you've already run :code:`nf-core launch KwanLab/Autometa` which created a :code:`nf-params.json` file) would look something like:

.. code-block:: bash
    
    nextflow run KwanLab/Autometa -params-file nf-params.json -c overwrite_config.config



For addtional information and examples see "Tuning workflow resources" `here <https://nf-co.re/usage/configuration#running-nextflow-on-your-system>`_



Additional autometa parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Up to date descriptions and default values of Autometa's nextflow parameters can be viewed using the following command: 

.. code-block:: bash

    nextflow run KwanLab/Autometa -r main --help






You can also adjust other pipeline parameters that ultimately control how the binning is performed.

*params.length_cutoff* : Smallest contig you want binned (default is 3000bp)

*params.kmer_size* : kmer size to use

*params.norm_method* : Which normalization method to use. See :ref:`advanced-usage-kmers` section for deails

*params.pca_dimensions* : Number of dimensions of which to reduce the initial k-mer frequencies matrix (default is 50). See :ref:`advanced-usage-kmers` section for deails

*params.embedding_method* :  Choices are "sksne", "bhsne", "umap" (default is bhsne) See :ref:`advanced-usage-kmers` section for deails

*params.embedding_dimensions* : Final dimensions of the kmer frequencies matrix (default is 2). See :ref:`advanced-usage-kmers` section for deails

*params.kingdom* : Bin contigs belonging to this kingdom. Choices are "bacteria" and "archaea" (default is bacteria). 

*params.clustering_method* : Cluster contigs using which clustering method. Choices are "dbscan" and "hdbscan" (default is "dbscan"). See :ref:`advanced-usage-binning` section for deails

*params.binning_starting_rank* : Which taxonomic rank to start the binning from. Choices are "superkingdom", "phylum", "class", "order", "family", "genus", "species" (default is "superkingdom"). See :ref:`advanced-usage-binning` section for deails

*params.classification_method* : Which clustering method to use for unclustered recruitment step. Choices are "decision_tree" and "random_forest" (default is "decision_tree"). See :ref:`advanced-usage-unclustered-recruitment` section for deails

*params.completeness* :  Minimum completeness needed to keep a cluster (default is atleast 20% complete). See :ref:`advanced-usage-binning` section for deails

*params.purity* : Minimum purity needed to keep a cluster (default is atleast 95% pure). See :ref:`advanced-usage-binning` section for deails

*params.cov_stddev_limit* : Which clusters to keep depending on the covergae std.dev (default is 25%). See :ref:`advanced-usage-binning` section for deails

*params.gc_stddev_limit* : Which clusters to keep depending on the GC% std.dev (default is 5%). See :ref:`advanced-usage-binning` section for deails


Customizing Autometa's Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


In case you want to tweak some of the scripts, run on your own scheduling system or modify the pipeline you can clone the repository and then run nextflow directly from the scripts as below:
.. code-block:: bash

    # Clone the autometa repository into current directory
    git clone git@github.com:KwanLab/Autometa.git 
    # Modify some code
    # Then run nextflow
    nextflow run $HOME/Autometa/nextflow

Without docker
^^^^^^^^^^^^^^

.. todo::

Useful options
^^^^^^^^^^^^^^

``-c`` : In case you have configured nextflow_ with your executor (see :ref:`Configure nextflow with your 'executor'`) and have made other modifications on how to run nextflow_ using your ``nexflow.config`` file, you can specify that file using the ``-c`` flag

To see all of the command line options available you can refer to `nexflow CLI documentation <https://www.nextflow.io/docs/latest/cli.html#command-line-interface-cli>`_

Resuming the workflow
^^^^^^^^^^^^^^^^^^^^^

One of the most powerful features of nextflow_ is resuming the workflow from the last completed process. If your pipeline was interrupted for some reason you can resume it from the last completed process using the resume flag (``-resume``). Eg, ``nextflow run KwanLab/Autometa -params-file nf-params.json -c my_other_parameters.config -resume``

Execution Report
^^^^^^^^^^^^^^^^

After running nextflow you can see the execution statistics of your autometa run, including the time taken, CPUs used, RAM used, etc separately for each process. Nextflow would generate a summary report, a timeline report and a trace report automatically for you in the ``"${params.tracedir}/pipeline_info`` directory (``"${params.tracedir}`` defaults to ``autometa_tracedir``). You can read more about these execution reports `here <https://www.nextflow.io/docs/latest/tracing.html#execution-report>`_. 

Workflow Visualized
^^^^^^^^^^^^^^^^^^^

You can also visualize the entire workflow ie. create the DAG from the written DOT file. Install `Graphviz <https://graphviz.org/>`_ and do ``dot -Tpng < pipeline_info/autometa-dot > autometa-dag.png`` to get the in the ``png`` format.

Configure nextflow with your 'executor'
---------------------------------------

.. todo::

For nextflow_ to run the Autometa pipeline through a job scheduler you will need to update the respective 'profile' section in nextflow's config file. Each 'profile' may be configured with any available scheduler as noted in the `nextflow executors docs <https://www.nextflow.io/docs/latest/executor.html>`_. By default nextflow_ will use your local computer as the 'executor'. The next section briefly walks through nextflow_ executor configuration to run with the slurm job scheduler.

We have prepared a template for ``nextflow.config`` which you can access from our GitHub repository using this `link <https://github.com/WiscEvan/Autometa/blob/4b4e3c60e076706e28deae4ae4d45f26b5df7dee/nextflow.config>`_. Go ahead and copy this file to your desired location and open it in your favorite text editor (eg. Vim, nano, VSCode, etc).


.. _using-slurm:

SLURM
^^^^^

This allows you to run the pipeline using the SLURM resource manager. To do this you'll first needed to identify the slurm partition to use. You can find the available slurm partitions by running ``sinfo``. Example: On running ``sinfo`` on our cluster we get the following:

.. image:: ../img/slurm_partitions.png
    :alt: Screen shot of ``sinfo`` output showing ``queue`` listed under partition  

The slurm partition available on our cluster is queue.  You'll need to update this in ``nextflow.config``. 

.. todo::
    Change the path to ``nextflow.config`` after the merge.

.. code-block:: groovy

    // Find this section of code in nextflow.config
    }
    cluster {
    process.executor = "slurm"
    // queue is the slurm partition to use in our case
    // Set SLURM partition with queue directive.
    process.queue = "queue" // <<-- change this to whatever your partition is called
    // See https://www.nextflow.io/docs/latest/executor.html#slurm for more details.
    }

More parameters that are available for the slurm executor are listed in the nextflow `executor docs for slurm <https://www.nextflow.io/docs/latest/executor.html#slurm>`_.

HTCondor
^^^^^^^^

This allows you to run the pipeline using the HTCondor resource manager. To do this you'll need to enable the HTCondor executor to condor value in the ``nextflow.config``.

.. code-block:: groovy

    // Find this section of code in nextflow.config
    }
    chtc {
        process.executor = "condor"
        // See https://www.nextflow.io/docs/latest/executor.html#htcondor for more configuration options.
    }

More parameters that are available for the htcondor executor are listed in the nextflow executor `docs for HTCondor <https://www.nextflow.io/docs/latest/executor.html#htcondor>`_.

.. note::
    1. The pipeline must be launched from a node where the ``condor_submit`` command is available, that is, in a common usage scenario, the cluster head node.
    2. The HTCondor executor for Nextflow_ does not support at this time the HTCondor ability to transfer input/output data to the corresponding job computing node. Therefore the data needs to be made accessible to the computing nodes using a shared file system directory from where the Nextflow_ workflow has to be executed (or specified via the -w option).











.. todo:: Below python specific maybe there should be two "running..." files, one for nextflow and one for python?



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


.. _nextflow: https://www.nextflow.io/
.. _Docker: https://www.docker.com/
.. _SPAdes: http://cab.spbu.ru/software/spades/
.. _Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _metaQuast: http://quast.sourceforge.net/metaquast
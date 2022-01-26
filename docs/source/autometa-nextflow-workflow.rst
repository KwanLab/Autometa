================================
🍏 Autometa Nextflow Workflow 🍏
================================

.. _autometa-nextflow-workflow:

Why nextflow?
#############

Nextflow helps Autometa produce reproducible results while allowing the pipeline to scale across different platforms and hardware.

System Requirements
###################

Currently the nextflow pipeline requires Docker 🐳 so it must be installed on your system.
If you don't have Docker installed you can install it from `docs.docker.com/get-docker <https://docs.docker.com/get-docker>`_.
We plan on removing this dependency in future versions, so that other dependency managers
(e.g. Conda, Singularity, etc) can be used.

Nextflow runs on any Posix compatible system. Detailed system requirements
can be found in the `nextflow documentation <https://www.nextflow.io/docs/latest/getstarted.html#requirements>`_

Nextflow (required) and nf-core tools (optional but highly recommended) installation will be discussed in :ref:`install-nextflow-nfcore-with-conda`.

🍏 Data Preparation 🍏
######################

#. [Metagenome](:ref:`metagenome-assembly`)
#. [Sample Sheet](:ref:`sample-sheet-preparation`)

.. _metagenome-assembly:

Metagenome Assembly
*******************

You will first need to assemble your shotgun metagenome, to provide to Autometa as input.

The following is a typical workflow for metagenome assembly:

#. Trim adapter sequences from the reads

   * We usually use Trimmomatic_.

#. Quality check the trimmed reads to ensure the adapters have been removed

   * We usually use FastQC_.

#. Assemble the trimmed reads

   * We usually use MetaSPAdes which is a part of the SPAdes_ package.

#. Check the quality of your assembly (Optional)

   * We usually use metaQuast_ for this (use ``--min-contig 1`` option to get an accurate N50).
    This tool can compute a variety of assembly statistics one of which is N50.
    This can often be useful for selecting an appropriate length cutoff value for pre-processing the metagenome.

.. _sample-sheet-preparation:

Preparing a Sample Sheet
************************

An example sample sheet for three possible ways to provide a sample as an input is provided below. The first example
provides a metagenome with paired-end read information, such that contig coverages may be determined using a read-based alignment
sub-workflow. The second example uses pre-calculated coverage information by providing a coverage table *with* the input metagenome assembly.
The third example retrieves coverage information from the assembly contig headers (Currently, this is only available with metagenomes assembled using SPAdes)

.. attention::
    If you have paired-end read information, you can supply these file paths within the sample sheet and the coverage
    table will be computed for you (See ``example_1`` in the example sheet below).

    If you have used any other assembler, then you may also provide a coverage table (See ``example_2`` in the example sheet below).
    Fortunately, Autometa can construct this table for you with: ``autometa-coverage``.
    Use ``--help`` to get the complete usage or for a few examples see :ref:`coverage-calculation`.

    If you use SPAdes then Autometa can use the k-mer coverage information in the contig names (``example_3`` in the example sample sheet below).

+-----------+--------------------------------------+----------------------------------------+----------------------------------------+-----------------------+-------------------------+
| sample    | assembly                             | fastq_1                                | fastq_2                                | coverage_tab          | cov_from_assembly       |
+===========+======================================+========================================+========================================+=======================+=========================+
| example_1 | /path/to/example/1/metagenome.fna.gz | /path/to/paired-end/fwd_reads.fastq.gz | /path/to/paired-end/rev_reads.fastq.gz |                       | 0                       |
+-----------+--------------------------------------+----------------------------------------+----------------------------------------+-----------------------+-------------------------+
| example_2 | /path/to/example/2/metagenome.fna.gz |                                        |                                        | /path/to/coverage.tsv | 0                       |
+-----------+--------------------------------------+----------------------------------------+----------------------------------------+-----------------------+-------------------------+
| example_3 | /path/to/example/3/metagenome.fna.gz |                                        |                                        |                       | spades                  |
+-----------+--------------------------------------+----------------------------------------+----------------------------------------+-----------------------+-------------------------+

.. note::

   To retrieve coverage information from a sample's contig headers, provide the ``assembler`` used for the sample value in the row under the ``cov_from_assembly`` column.
   Using a ``0`` will designate to the workflow to try to retrieve coverage information from the coverage table (if it is provided)
   or coverage information will be calculated by read alignments using the provided paired-end reads. If both paired-end reads and a coverage table are provided,
   the pipeline will prioritize the coverage table.

   If you are providing a coverage table to ``coverage_tab`` with your input metagenome, it must be tab-delimited and contain *at least* the header columns, ``contig`` and ``coverage``.

Supported Assemblers for ``cov_from_assembly``
----------------------------------------------

+--------------+-----------------+-----------------------+
|    Assembler | Supported (Y/N) | ``cov_from_assembly`` |
+==============+=================+=======================+
| [meta]SPAdes |        Y        |       ``spades``      |
+--------------+-----------------+-----------------------+
|    Unicycler |        N        |     ``unicycler``     |
+--------------+-----------------+-----------------------+
|      Megahit |        N        |      ``megahit``      |
+--------------+-----------------+-----------------------+


You may copy the below table as a csv and paste it into a file to begin your sample sheet. You will need to update your input parameters, accordingly.

Example ``sample_sheet.csv``
----------------------------

.. code-block:: bash

    sample,assembly,fastq_1,fastq_2,coverage_tab,cov_from_assembly
    example_1,/path/to/example/1/metagenome.fna.gz,/path/to/paired-end/fwd_reads.fastq.gz,/path/to/paired-end/rev_reads.fastq.gz,,0
    example_2,/path/to/example/2/metagenome.fna.gz,,,/path/to/coverage.tsv,0
    example_3,/path/to/example/3/metagenome.fna.gz,,,,spades

.. caution::

    Paths to any of the file inputs **must be absolute file paths**.

    e.g.

    #. Replacing any instance of the ``$HOME`` variable with the real path

        (``$HOME/Autometa/tests/data/metagenome.fna.gz`` --> ``/home/user/Autometa/tests/data/metagenome.fna.gz``)

    #. Using the entire file path of the input

        (``tests/data/metagenome.fna.gz`` --> ``/home/user/Autometa/tests/data/metagenome.fna.gz``)

Basic
#####

While the Autometa Nextflow pipeline can be run using Nextflow directly, we designed
it using nf-core standards and templating to provide an easier user experience through
use of the nf-core "tools" python library. The directions below demonstrate using a minimal
Conda environment to install Nextflow and nf-core tools and then running the Autometa pipeline.

.. _install-nextflow-nfcore-with-conda:

Installing Nextflow and nf-core tools with Conda
************************************************

If you have not previously installed/used Conda, you can get it using the
Miniconda installer appropriate to your system, here: `<https://docs.conda.io/en/latest/miniconda.html>`_

After installing conda, running the following command will create a minimal
Conda environment named "autometa-nf", and install Nextflow and nf-core tools.

.. code-block:: bash

    conda env create --file=https://raw.githubusercontent.com/KwanLab/Autometa/main/environment.yml

If you receive the message...

.. code-block:: bash

    CondaValueError: prefix already exists:

...it means you have already created the environment. If you want to overwrite/update
the environment then add the :code:`--force` flag to the end of the command.

.. code-block:: bash

    conda env create --file=https://raw.githubusercontent.com/KwanLab/Autometa/main/environment.yml --force

Once Conda has finished creating the environment be sure to activate it:

.. code-block:: bash

    conda activate autometa-nf


Using nf-core
*************

Download/Launch the Autometa Nextflow pipeline using nf-core tools.
The stable version of Autometa will always be the "main" git branch.
To use an in-development git branch switch "main" in the command with
the name of the desired branch. After the pipeline downloads, nf-core will
start the pipeline launch process.

.. code-block:: bash

    nf-core launch KwanLab/Autometa -r main

You will then be asked to choose "Web based" or "Command line" for selecting/providing options.
While it is possible to use the command line version, it is preferred and easier to use the web-based GUI.
Use the arrow keys to select one or the other and then press return/enter.


Setting parameters with a web-based GUI
***************************************

The GUI will present all available parameters, though some extra
parameters may be hidden (these can be revealed by selecting
"Show hidden params" on the right side of the page).

Required parameters
*******************

1. :code:`--input`: the path to your input sample sheet
2. :code:`-profile`: this sets options specified within the "profiles" section in the pipeline's nextflow.config file
    - :code:`standard` (default): runs all process jobs locally, (currently this requires Docker).
    - :code:`slurm`: submits all process jobs into the slurm queue. See :ref:`using-slurm` before using

.. caution::

    Notice the number of hyphens used between ``--input`` and ``-profile``. ``--input`` is an `Autometa` workflow parameter
    where as ``-profile`` is a `nextflow` argument. This difference in hyphens is true for passing in all arguments to the `Autometa`
    workflow and `nextflow`, respectively.

Running the pipeline
********************

After you are finished double-checking your parameter settings, click "Launch"
at the top right of web based GUI page, or "Launch workflow" at the bottom of
the page. After returning to the terminal you should be provided the option
:code:`Do you want to run this command now?  [y/n]`  enter :code:`y` to begin the pipeline.

.. note::

    This process will lead to nf-core tools creating a file named :code:`nf-params.json`.
    This file contains your specified parameters that differed from the pipeline's defaults.
    This file can also be manually modified and/or shared to allow reproducible configuration
    of settings (e.g. among members within a lab sharing the same server).

    Additionally all Autometa specific pipeline parameters can be used as command line arguments
    using the :code:`nextflow run ...` command by prepending the parameter name with two hyphens
    (e.g. :code:`--outdir /path/to/output/workflow/results`)


Advanced
########

Parallel computing and computer resource allotment
**************************************************

While you might want to provide Autometa all the compute resources available in order to get results
faster, that may or may not actually achieve the fastest run time.

Within the Autometa pipeline, parallelization happens by providing all the assemblies at once
to software that internally handles parallelization.

The Autometa pipeline will try and use all resources available to individual
pipeline modules. Each module/process has been pre-assigned resource allotments via a low/medium/high tag.
This means that even if you don't select for the pipeline to run in parallel some modules (e.g. DIAMOND BLAST)
may still use multiple cores.

* The maximum number of CPUs that any single module can use is defined with the :code:`--max_cpus` option (default: 4).
* You can also set :code:`--max_memory` (default: 16GB)
* :code:`--max_time` (default: 240h). :code:`--max_time` refers to the maximum time *each process* is allowed to run, *not* the execution time for the the entire pipeline.

Databases
*********

Autometa uses the following NCBI databases throughout its pipeline:

- Non-redundant nr database
    - `ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz <https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz>`_
- prot.accession2taxid.gz
    - `ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz <https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz>`_
- nodes.dmp, names.dmp and merged.dmp - Found within
    - `ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz <ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_

If you are running autometa for the first time you'll have to download these databases.
You may use ``autometa-update-databases --update-ncbi``. This will download the databases to the default path. You can check
the default paths using ``autometa-config --print``. If you need to change the default download directory you can use
``autometa-config --section databases --option ncbi --value <path/to/new/ncbi_database_directory>``.
See ``autometa-update-databases -h`` and ``autometa-config -h`` for full list of options.

In your ``nf-params.json`` file you also need to specify the directory where the different databases are present.
Make sure that the directory path contains the following databases:

- Diamond formatted nr file => nr.dmnd
- Extracted files from tarball taxdump.tar.gz
- prot.accession2taxid.gz

.. code-block::
    {
        "single_db_dir" = "$HOME/Autometa/autometa/databases/ncbi"
    }

.. note::
    Find the above section of code in ``nf-params.json`` and update this path to the folder
    with all of the downloaded/formatted NCBI databases.

CPUs, Memory, Disk
******************

.. note::

    Like nf-core pipelines, we have set some automatic defaults for Autometa's processes. These are dynamic and each
    process will try a second attempt using more resources if the first fails due to resources. Resources are always
    capped by the parameters (show with defaults):

    - :code:`--max_cpus = 2`
    - :code:`--max_memory = 6.GB`
    - :code:`--max_time = 48.h`

The best practice to change the resources is to create a new config file and point to it at runtime by adding the
flag :code:`-c path/to/config_file`


For example, to give all resource-intensive (i.e. having ``label process_high``) jobs additional memory and cpus, create a file called :code:`process_high_mem.config` and insert

.. code-block:: groovy

    process {
        withLabel:process_high {
            memory = 200.GB
            cpus = 32
        }
    }

Then your command to run the pipeline (assuming you've already run :code:`nf-core launch KwanLab/Autometa` which created
a :code:`nf-params.json` file) would look something like:

.. code-block:: bash

    nextflow run KwanLab/Autometa -params-file nf-params.json -c process_high_mem.config

.. caution::

    If you are restarting from a previous run, **DO NOT FORGET** to also add the ``-resume`` flag to the nextflow run command.

    **Notice only 1 hyphen is used with the :code:`-resume` nextflow parameter**


For additional information and examples see `Tuning workflow resources <https://nf-co.re/usage/configuration#running-nextflow-on-your-system>`_

Additional Autometa parameters
******************************

Up to date descriptions and default values of Autometa's nextflow parameters can be viewed using the following command:

.. code-block:: bash

    nextflow run KwanLab/Autometa -r main --help


You can also adjust other pipeline parameters that ultimately control how binning is performed.

``params.length_cutoff`` : Smallest contig you want binned (default is 3000bp)

``params.kmer_size`` : kmer size to use

``params.norm_method`` : Which kmer frequency normalization method to use. See
:ref:`advanced-usage-kmers` section for details

``params.pca_dimensions`` : Number of dimensions of which to reduce the initial k-mer frequencies
matrix (default is ``50``). See :ref:`advanced-usage-kmers` section for details

``params.embedding_method`` :  Choices are ``sksne``, ``bhsne``, ``umap``, ``densmap``, ``trimap``
(default is ``bhsne``) See :ref:`advanced-usage-kmers` section for details

``params.embedding_dimensions`` : Final dimensions of the kmer frequencies matrix (default is ``2``).
See :ref:`advanced-usage-kmers` section for details

``params.kingdom`` : Bin contigs belonging to this kingdom. Choices are ``bacteria`` and ``archaea``
(default is ``bacteria``).

``params.clustering_method`` : Cluster contigs using which clustering method. Choices are "dbscan" and "hdbscan"
(default is "dbscan"). See :ref:`advanced-usage-binning` section for details

``params.binning_starting_rank`` : Which taxonomic rank to start the binning from. Choices are ``superkingdom``, ``phylum``,
``class``, ``order``, ``family``, ``genus``, ``species`` (default is ``superkingdom``). See :ref:`advanced-usage-binning` section for details

``params.classification_method`` : Which clustering method to use for unclustered recruitment step.
Choices are ``decision_tree`` and ``random_forest`` (default is ``decision_tree``). See :ref:`advanced-usage-unclustered-recruitment` section for details

``params.completeness`` :  Minimum completeness needed to keep a cluster (default is at least 20% complete, e.g. ``20``).
See :ref:`advanced-usage-binning` section for details

``params.purity`` : Minimum purity needed to keep a cluster (default is atleast 95% pure, e.g. ``95``).
See :ref:`advanced-usage-binning` section for details

``params.cov_stddev_limit`` : Which clusters to keep depending on the covergae std.dev (default is 25%, e.g. ``25``).
See :ref:`advanced-usage-binning` section for details

``params.gc_stddev_limit`` : Which clusters to keep depending on the GC% std.dev (default is 5%, e.g. ``5``).
See :ref:`advanced-usage-binning` section for details


Customizing Autometa's Scripts
******************************

In case you want to tweak some of the scripts, run on your own scheduling system or modify the pipeline you can clone
the repository and then run nextflow directly from the scripts as below:

.. code-block:: bash

    # Clone the autometa repository into current directory
    git clone git@github.com:KwanLab/Autometa.git
    # Modify some code
    # e.g. one of the local modules
    code /$HOME/Autometa/modules/local/align_reads.nf
    # Then run nextflow
    nextflow run $HOME/Autometa

Useful options
**************

``-c`` : In case you have configured nextflow with your executor (see :ref:`Configuring your nextflow Executor`)
and have made other modifications on how to run nextflow using your ``nexflow.config`` file, you can specify that file
using the ``-c`` flag

To see all of the command line options available you can refer to
`nexflow CLI documentation <https://www.nextflow.io/docs/latest/cli.html#command-line-interface-cli>`_

Resuming the workflow
*********************

One of the most powerful features of nextflow is resuming the workflow from the last completed process. If your pipeline
was interrupted for some reason you can resume it from the last completed process using the resume flag (``-resume``).
Eg, ``nextflow run KwanLab/Autometa -params-file nf-params.json -c my_other_parameters.config -resume``

Execution Report
****************

After running nextflow you can see the execution statistics of your autometa run, including the time taken, CPUs used,
RAM used, etc separately for each process. Nextflow will generate summary, timeline and trace reports automatically for
you in the ``${params.outdir}/trace`` directory. You can read more about this in the
`nextflow docs on execution reports <https://www.nextflow.io/docs/latest/tracing.html#execution-report>`_.

Visualizing the Workflow
------------------------

You can visualize the entire workflow ie. create the directed acyclic graph (DAG) of processes from the written DOT file. First install
`Graphviz <https://graphviz.org/>`_ (``conda install -c anaconda graphviz``) then do ``dot -Tpng < pipeline_info/autometa-dot > autometa-dag.png`` to get the
in the ``png`` format.

Configuring your nextflow Executor
**********************************

For nextflow to run the Autometa pipeline through a job scheduler you will need to update the respective ``profile``
section in nextflow's config file. Each ``profile`` may be configured with any available scheduler as noted in the
`nextflow executors docs <https://www.nextflow.io/docs/latest/executor.html>`_. By default nextflow will use your
local computer as the 'executor'. The next section briefly walks through nextflow executor configuration to run
with the slurm job scheduler.

We have prepared a template for ``nextflow.config`` which you can access from the KwanLab/Autometa GitHub repository using this
`nextflow.config template <https://raw.githubusercontent.com/KwanLab/Autometa/main/nextflow.config>`_. Go ahead
and copy this file to your desired location and open it in your favorite text editor (eg. Vim, nano, VSCode, etc).


.. _using-slurm:

SLURM
-----

This allows you to run the pipeline using the SLURM resource manager. To do this you'll first needed to identify the
slurm partition to use. You can find the available slurm partitions by running ``sinfo``. Example: On running ``sinfo``
on our cluster we get the following:

.. image:: ../img/slurm_partitions.png
    :alt: Screen shot of ``sinfo`` output showing ``queue`` listed under partition

The slurm partition available on our cluster is ``queue``.  You'll need to update this in ``nextflow.config``.

.. code-block:: groovy

    profiles {
        // Find this section of code in nextflow.config
        slurm {
            process.executer       = "slurm"
            // NOTE: You can determine your slurm partition (e.g. process.queue) with the `sinfo` command
            // Set SLURM partition with queue directive.
            process.queue = "queue" // <<-- change this to whatever your partition is called
            // queue is the slurm partition to use in our case
            docker.enabled         = true
            docker.userEmulation   = true
            singularity.enabled    = false
            podman.enabled         = false
            shifter.enabled        = false
            charliecloud.enabled   = false
            executor {
                queueSize = 8
            }
        }
    }

More parameters that are available for the slurm executor are listed in the nextflow
`executor docs for slurm <https://www.nextflow.io/docs/latest/executor.html#slurm>`_.


Docker image selection
######################


Especially when developing new features it may be necessary to run the pipeline with a custom docker image.
Create a new image by navigating to the top Autometa directory and running ``make image``. This will create a new
Autometa Docker image, tagged with the name of the current Git branch.

To use this tagged version (or any other Autometa image tag) add the argument ``--autometa_image tag_name`` to the nextflow run command


Running modules
###############

Many of the Autometa modules may be run standalone.

Simply pass in the ``-m`` flag when calling a script to signify to python you are
running the script as an Autometa *module*.

I.e. ``python -m autometa.common.kmers -h``

.. note::
    Autometa has many *entrypoints* available that are utilized by the nextflow workflow. If you have installed autometa,
    all of these entrypoints will be available to you.

Using Autometa's Python API
###########################

Autometa's classes and functions are available after installation.
To access these, do the same as importing any other python library.

Examples
********


Samtools wrapper
----------------

To incorporate a call to ``samtools sort`` inside of your python code using the Autometa samtools wrapper.

.. code-block:: python

    from autometa.common.external import samtools

    # To see samtools.sort parameters try the commented command below:
    # samtools.sort?

    # Run samtools sort command in ipython interpreter
    samtools.sort(sam="<path/to/alignment.sam>", out="<path/to/output/alignment.bam>", cpus=4)


Metagenome Description
----------------------

Here is an example to easily assess your metagenome's characteristics using Autometa's Metagenome class

.. code-block:: python

    from autometa.common.metagenome import Metagenome

    # To see input parameters, instance attributes and methods
    # Metagenome?

    # Create a metagenome instance
    mg = Metagenome(assembly="/path/to/metagenome.fasta")

    # To see available methods (ignore any elements in the list with a double underscore)
    dir(mg)

    # Get pandas dataframe of metagenome details.
    metagenome_df = mg.describe()

    metagenome_df.to_csv("path/to/metagenome_description.tsv", sep='\t', index=True, header=True)


k-mer frequency counting, normalization, embedding
--------------------------------------------------

To quickly perform a k-mer frequency counting, normalization and embedding pipeline...

.. code-block:: python

    from autometa.common import kmers

    # Count kmers
    counts = kmers.count(
        assembly="/path/to/metagenome.fasta",
        size=5
    )

    # Normalize kmers
    norm_df = kmers.normalize(
        df=counts,
        method="ilr"
    )

    # Embed kmers
    embed_df = kmers.embed(
        norm_df,
        pca_dimensions=50,
        embed_dimensions=3,
        method="densmap"
    )


.. _nextflow: https://www.nextflow.io/
.. _Docker: https://www.docker.com/
.. _SPAdes: http://cab.spbu.ru/software/spades/
.. _Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _metaQuast: http://quast.sourceforge.net/metaquast

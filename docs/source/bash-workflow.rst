.. _autometa-bash-workflow:

===================
🐚 Bash Workflow 🐚
===================

Getting Started
###############

#. :ref:`Compute Environment Setup`
#. :ref:`Download Workflow Template`
#. :ref:`Configure Required Inputs`

Compute Environment Setup
*************************

If you have not previously installed/used Conda, you can get it using the
Miniconda installer appropriate to your system, here: `<https://docs.conda.io/en/latest/miniconda.html>`_

You may either create a new Conda environment named "autometa"...

.. code-block:: bash

    conda create -n autometa -c bioconda autometa
    # Then, once Conda has finished creating the environment
    # you may activate it:
    conda activate autometa

\.\.\. or install Autometa into any of your existing environments.

This installs Autometa in your current active environment:

.. code-block:: bash

    conda install -c bioconda autometa

The next command installs Autometa in the provided environment:

.. code-block:: bash

    conda install -n <your-env-name> -c bioconda autometa

Download Workflow Template
**************************

To run Autometa using the bash workflow you will simply need to download and configure the workflow template to your
metagenomes specifications.

* `autometa.sh <https://github.com/KwanLab/Autometa/blob/main/workflows/autometa.sh>`_
* `autometa-large-data-mode.sh <https://github.com/KwanLab/Autometa/blob/main/workflows/autometa-large-data-mode.sh>`_

Here are a few download commands if you do not want to navigate to the workflow on GitHub

via curl
--------

.. code-block:: bash

    curl -o autometa.sh https://raw.githubusercontent.com/KwanLab/Autometa/main/workflows/autometa.sh

via wget
--------

.. code-block:: bash

    wget https://raw.githubusercontent.com/KwanLab/Autometa/main/workflows/autometa.sh

.. note::

    The ``autometa-large-data-mode`` workflow is also available and is configured similarly to the ``autometa`` bash workflow.

Configure Required Inputs
*************************

The Autometa bash workflow requires the following input file and directory paths. To see how to prepare each input, see :ref:`bash-workflow-data-preparation`

#. Assembly (``assembly``)
#. Alignments (``bam``)
#. ORFs (``orfs``)
#. Diamond blastp results table (``blast``)
#. NCBI database directory (``ncbi``)
#. Input sample name (``simpleName``)
#. Output directory (``outdir``)

.. _bash-workflow-data-preparation:

Data preparation
################

#. :ref:`metagenome-preparation` (``assembly``)
#. :ref:`alignments-preparation` (``bam``)
#. :ref:`orfs-preparation` (``orfs``)
#. :ref:`blastp-preparation` (``blast``)
#. :ref:`ncbi-preparation` (``ncbi``)

.. _metagenome-preparation:

Metagenome Assembly
*******************

You will first need to assemble your shotgun metagenome, to provide to Autometa as input.

The following is a typical workflow for metagenome assembly:

#. Trim adapter sequences from the reads

    We usually use Trimmomatic_.

#. Quality check the trimmed reads to ensure the adapters have been removed

    We usually use FastQC_.

#. Assemble the trimmed reads

    We usually use MetaSPAdes which is a part of the SPAdes_ package.

#. Check the quality of your assembly (Optional)

    We usually use metaQuast_ for this (use ``--min-contig 1`` option to get an accurate N50).
    This tool can compute a variety of assembly statistics one of which is N50.
    This can often be useful for selecting an appropriate length cutoff value for pre-processing the metagenome.

.. _alignments-preparation:

Alignments Preparation
**********************

.. note::
    The following example requires ``bwa``, ``kart`` and ``samtools``

    ``conda install -c bioconda bwa kart samtools``

.. code-block:: bash

    # First index metagenome assembly
    bwa index \
        -b 550000000000 \ # block size for the bwtsw algorithm (effective with -a bwtsw) [default=10000000]
        metagenome.fna     # Path to input metagenome

    # Now perform alignments (we are using kart, but you can use another alignment tool if you'd like)
    kart \
        -i metagenome.fna                   \ # Path to input metagenome
        -t 20                               \ # Number of cpus to use
        -f /path/to/forward_reads.fastq.gz  \ # Path to forward paired-end reads
        -f2 /path/to/reverse_reads.fastq.gz \ # Path to reverse paired-end reads
        -o alignments.sam                      # Path to alignments output

    # Now sort alignments and convert to bam format
    samtools sort \
        -@ 40              \ # Number of cpus to use
        -m 10G             \ # Amount of memory to use
        alignments.sam     \ # Input alignments file path
        -o alignments.bam     # Output alignments file path

.. _orfs-preparation:

ORFs
****

.. note::
    The following example requires ``prodigal``. e.g. ``conda install -c bioconda prodigal``

.. code-block:: bash

    prodigal -i metagenome.fna \
        -f "gbk" \
        -d "metagenome.orfs.fna" \
        -o "metagenome.orfs.gbk" \
        -a "metagenome.orfs.faa" \ # This generated file is required as input to the bash workflow
        -s "metagenome.all_orfs.txt"

.. _blastp-preparation:

Diamond blastp Preparation
**************************

.. note::
    The following example requires ``diamond``. e.g. ``conda install -c bioconda diamond``

.. code-block:: bash

    diamond blastp \
        --query "metagenome.orfs.faa" \ # See prodigal output from above
        --db /path/to/nr.dmnd         \ # See NCBI section
        --threads <num cpus to use>   \
        --out blastp.tsv # This generated file is required as input to the bash workflow

.. _ncbi-preparation:

NCBI Preparation
****************

If you are running Autometa for the first time you'll have to download the NCBI databases.

.. code-block:: bash

    # First configure where you want to download the NCBI databases
    autometa-config \
        --section databases \
        --option ncbi \
        --value <path/to/your/ncbi/database/directory>

    # Now download and format the NCBI databases
    autometa-update-databases --update-ncbi

.. note::
    You can check the default config paths using ``autometa-config --print``.

    See ``autometa-update-databases -h`` and ``autometa-config -h`` for full list of options.

The previous command will download the following NCBI databases:

- Non-redundant nr database
    - `ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz <https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz>`_
- prot.accession2taxid.gz
    - `ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz <https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz>`_
- nodes.dmp, names.dmp and merged.dmp - Found within
    - `ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz <ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_

Input Sample Name
*****************

A crucial step prior to running the Autometa bash workflow is specifying the metagenome sample name and where to output
Autometa's results.

.. code-block:: bash

    # Default
    simpleName="TemplateAssemblyName"
    # Replace with your sample name
    simpleName="MySample"

.. note::

    The ``simpleName`` that is provided will be used as a prefix to all of the resulting autometa output files.

Output directory
****************

Immediately following the ``simpleName`` parameter, you will need to specify where to write all results.

.. code-block:: bash

    # Default
    outdir="AutometaOutdir"
    # Replace with your output directory...
    outdir="MySampleAutometaResults"

Running the pipeline
####################

After you are finished configuring/double-checking your parameter settings..

You may run the pipeline via bash:

.. code-block:: bash

    bash autometa.sh

or submit the pipeline into a queue:

For example, with slurm:

.. code-block:: bash

    sbatch autometa.sh

.. caution::

    Make sure your conda autometa environment is activated or the autometa entrypoints will not be available.

Additional parameters
#####################

You can also adjust other pipeline parameters that ultimately control how binning is performed.
These are located at the top of the workflow just under the required inputs.

``length_cutoff`` : Smallest contig you want binned (default is 3000bp)

``kmer_size`` : kmer size to use

``norm_method`` : Which kmer frequency normalization method to use. See
:ref:`advanced-usage-kmers` section for details

``pca_dimensions`` : Number of dimensions of which to reduce the initial k-mer frequencies
matrix (default is ``50``). See :ref:`advanced-usage-kmers` section for details

``embed_method`` :  Choices are ``sksne``, ``bhsne``, ``umap``, ``densmap``, ``trimap``
(default is ``bhsne``) See :ref:`advanced-usage-kmers` section for details

``embed_dimensions`` : Final dimensions of the kmer frequencies matrix (default is ``2``).
See :ref:`advanced-usage-kmers` section for details

``cluster_method`` : Cluster contigs using which clustering method. Choices are "dbscan" and "hdbscan"
(default is "dbscan"). See :ref:`advanced-usage-binning` section for details

``binning_starting_rank`` : Which taxonomic rank to start the binning from. Choices are ``superkingdom``, ``phylum``,
``class``, ``order``, ``family``, ``genus``, ``species`` (default is ``superkingdom``). See :ref:`advanced-usage-binning` section for details

``classification_method`` : Which clustering method to use for unclustered recruitment step.
Choices are ``decision_tree`` and ``random_forest`` (default is ``decision_tree``). See :ref:`advanced-usage-unclustered-recruitment` section for details

``completeness`` :  Minimum completeness needed to keep a cluster (default is at least 20% complete, e.g. ``20``).
See :ref:`advanced-usage-binning` section for details

``purity`` : Minimum purity needed to keep a cluster (default is at least 95% pure, e.g. ``95``).
See :ref:`advanced-usage-binning` section for details

``cov_stddev_limit`` : Which clusters to keep depending on the coverage std.dev (default is 25%, e.g. ``25``).
See :ref:`advanced-usage-binning` section for details

``gc_stddev_limit`` : Which clusters to keep depending on the GC% std.dev (default is 5%, e.g. ``5``).
See :ref:`advanced-usage-binning` section for details

.. note::

    If you are configuring an autometa job using the ``autometa-large-data-mode.sh`` template,
    there will be an additional parameter called, ``max_partition_size`` (default=10,000). This is the maximum size
    partition the Autometa clustering algorithm will consider. Any taxon partitions larger than this setting
    will be skipped.

.. _SPAdes: http://cab.spbu.ru/software/spades/
.. _Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _metaQuast: http://quast.sourceforge.net/metaquast

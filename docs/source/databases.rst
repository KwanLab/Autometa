=========
Databases
=========

If you are running Autometa for the first time you will need to download and format a few databases.
You may do this manually or using a few Autometa helper scripts. If you would like to use Autometa's
scripts for this, you will first need to install Autometa (See :ref:`Installation`).

The following sections use a pair of commands to configure autometa such that the database is updated
according to its respective path.

Markers
#######

.. code-block:: bash
    
    # Point Autometa to where you would like your markers database directory
    autometa-config \
        --section databases --option markers \
        --value <path/to/your/markers/database/directory>
    
    # Update your markers database directory
    autometa-update-databases --update-markers

.. alert::
    
    Do NOT use a trailing slash, e.g. NO ``/`` for the database directory paths!

Links to these markers files and their associated cutoff values are below:

- bacteria single-copy-markers - `link <https://raw.githubusercontent.com/KwanLab/Autometa/main/autometa/databases/markers/bacteria.single_copy.hmm>`__
- bacteria single-copy-markers cutoffs - `link <https://raw.githubusercontent.com/KwanLab/Autometa/main/autometa/databases/markers/bacteria.single_copy.cutoffs>`__
- archaea single-copy-markers - `link <https://raw.githubusercontent.com/KwanLab/Autometa/main/autometa/databases/markers/archaea.single_copy.hmm>`__
- archaea single-copy-markers cutoffs - `link <https://raw.githubusercontent.com/KwanLab/Autometa/main/autometa/databases/markers/archaea.single_copy.cutoffs>`__

NCBI
####

.. code-block:: bash

    # First configure where you want to download the NCBI databases
    autometa-config \
        --section databases --option ncbi \
        --value <path/to/your/ncbi/database/directory>

    # Now download and format the NCBI databases
    autometa-update-databases --update-ncbi

.. note::

    You can check the config paths using ``autometa-config --print``.

    See ``autometa-update-databases -h`` and ``autometa-config -h`` for full list of options.

The previous command will download the following NCBI databases:

- Non-redundant nr database
    - `ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz <https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz>`_
- prot.accession2taxid.gz
    - `ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz <https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz>`_
- nodes.dmp, names.dmp, merged.dmp and delnodes.dmp - Found within
    - `ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz <ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_

After these files are downloaded, the ``taxdump.tar.gz`` tarball's files are extracted and the non-redundant protein database (``nr.gz``)
is formatted as a diamond database (i.e. ``nr.dmnd``). This will significantly speed-up the ``diamond blastp`` searches.

Genome Taxonomy Database (GTDB)
###############################

If you would like to incorporate the benefits of using the Genome Taxonomy Database,
you can either run the following script or manually download the respective databases.

.. code-block:: bash

    # First configure where you want to download the GTDB databases
    autometa-config \
        --section databases --option gtdb \
        --value <path/to/your/gtdb/database/directory>

    # To use a specific GTDB release
    autometa-config \
        --section gtdb --option release \
        --value latest
        # Or --value r207 or --value r202, etc.

    # Download and format the configured GTDB databases release
    autometa-update-databases --update-gtdb


.. note::

    You can check the default config paths using ``autometa-config --print``.

    See ``autometa-update-databases -h`` and ``autometa-config -h`` for full list of options.

The previous command will download the following GTDB databases and format the `gtdb_proteins_aa_reps.tar.gz` to generate `gtdb.dmnd` to be used by Autometa:

- Amino acid sequences of representative genome
    - `gtdb_proteins_aa_reps.tar.gz <https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz>`_
- gtdb-taxdump.tar.gz from `shenwei356/gtdb-taxdump <https://github.com/shenwei356/gtdb-taxdump/releases>`_
    - `gtdb-taxdump.tar.gz <https://github.com/shenwei356/gtdb-taxdump/releases/latest/download/gtdb-taxdump.tar.gz>`_


Once unzipped `gtdb-taxdump.tar.gz` will have the taxdump files of all the respective GTDB releases. 
Make sure that the release you use is in line with the `gtdb_proteins_aa_reps.tar.gz` release version. 
It's better to always use the latest version. 

All the taxonomy files for a specific taxonomy database should be in a single directory. 
You can now copy the taxdump files of the desired release version in the sample directory as `gtdb.dmnd`

Alternatively if you have manually downloaded `gtdb_proteins_aa_reps.tar.gz` and `gtdb-taxdump.tar.gz` you can run the 
following command to format the `gtdb_proteins_aa_reps.tar.gz` to generate `gtdb.dmnd` and make it ready for Autometa.

.. code-block:: bash

    autometa-setup-gtdb --reps-faa <path/to/gtdb_proteins_aa_reps.tar.gz> --dbdir <path/to/output_directory> --cpus 20

.. note::

    Again Make sure that the formatted `gtdb_proteins_aa_reps.tar.gz` database and gtdb taxdump files are in the same directory. 

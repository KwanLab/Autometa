===========
bedtools.py
===========

File containing wrapper functions for bedtools.

**Running the script as stand alone module**
While in the Autometa directory</path/to/output/file/prots.out> and with the conda environment activated, run: 

.. code-block:: bash

    python -m autometa.common.external.prodigal </path/to/BAM/alignment.bam> 
    </path/to/genome/output/lengths.tsv> </path/to/output/alignment.bed> --[options]

**Inputs** : ibam file

**Returns** : A tab-delimited table with index = contig and column = coverage

.. rubric:: Usage and command line options

.. program-output:: cd ../.. ; python -m autometa.common.external.bedtools -h
    :shell:

.. rubric:: Functions

.. autosummary::

    bedtools.genomecov
    bedtools.parse

.. automodule:: bedtools
    :members:
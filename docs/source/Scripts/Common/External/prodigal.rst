===========
prodigal.py
===========

File containing functions to retrieve orfs from provided assembly using prodigal.

**Running the script as stand alone module**
While in the Autometa directory and with the conda environment activated, run: 

.. code-block:: bash

    python -m autometa.common.external.prodigal </path/to/assembly/input/file> 
    </path/to/output/file/nucls.out> </path/to/output/file/prots.out> --[options]

**Inputs** : Assembly file (fasta format)

**Returns** : Nucleotide file (ffn format) and protein file (faa format) containing orfs

.. rubric:: Usage and command line options

.. program-output:: cd ../.. ; python -m autometa.common.external.prodigal -h
    :shell:

.. rubric:: Functions

.. autosummary::
   
    prodigal.contigs_from_headers
    prodigal.orf_records_from_contigs
    prodigal.run
   

.. automodule:: prodigal
    :members:

==========
diamond.py
==========

File containing class and functions related to running diamond on metagenome sequences

**Running the script as stand alone module**
While in the Autometa directory and with the conda environment activated, run:

.. code-block:: bash

    python -m autometa.common.external.diamond <path/to/input/file.faa> 
    </path/to/diamond/formatted/database> </path/to/ncbi/prot.accession2taxid.gz> 
    </path/to/diamond/output/file> <blastp or blastx> --[options]

**Inputs** : Protein file (faa format)

**Returns** : A dictionary of Blastp hits, with key = qseqid and value = Diamond result

.. rubric:: Usage and command line options

.. program-output:: cd ../.. ; python -m autometa.common.external.diamond -h
    :shell:

.. rubric:: Classes

.. autosummary::
   
      diamond.DiamondResult

.. automodule:: diamond
    :members:
=========
kmers.py
=========

File containing functions to count and retrieve k-mers of a given sequences.

**Running the script as stand alone module**
While in the Autometa directory and with the conda environment activated, run: 

.. code-block:: bash

  python -m autometa.common.kmers --fasta <path/to/input/fasta/file.fasta> 
  --kmers <path/to/output/file.tsv> --[options]

**Inputs** : Assembly file (in fna format), k-mer size

**Returns** : A tab-delimited matrix of k-mer frequency counts with an index column corresponding
to each contig in the input fasta file

.. rubric:: Usage and command line options

.. program-output:: cd ../.. ; python -m autometa.common.kmers -h
    :shell:

.. rubric:: Functions

.. autosummary::
   
  kmers.count
  kmers.embed
  kmers.init_kmers
  kmers.load
  kmers.mp_counter
  kmers.normalize
  kmers.record_counter
  kmers.revcomp
  kmers.seq_counter

.. automodule:: kmers
    :members:


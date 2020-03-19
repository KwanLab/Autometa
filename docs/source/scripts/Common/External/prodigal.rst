Prodigal
========

File containing functions to retrieve orfs from provided assembly using prodigal.

**Running the script as stand alone module**
While in the Autometa directory and with the conda environment activated, run: 
``python3 -m autometa.common.external.prodigal </path/to/assembly/input/file> </path/to/output/file/nucls.out> </path/to/output/file/prots.out> --[options]``

**Inputs** : Assembly file (fasta format)

**Returns** : Nucleotide file (ffn format) and protein file (faa format) containing orfs

**Usage and Command line options:**

.. code-block:: bash

    python -m autometa.common.external.prodigal -h

    Usage: Calls ORFs with provided input assembly [-h] [--force] [--cpus CPUS]
                                                        [--parallel] [--verbose]
                                                        assembly nucls_out prots_out

    Positional arguments:
    assembly     </path/to/assembly/input/file>
    nucls_out    </path/to/output/file/nucls.out>
    prots_out    </path/to/output/file/prots.out>

    Optional arguments:
    -h, --help   show this help message and exit
    --force      force overwrite of ORFs out filepaths
    --cpus       number of cpus to use
    --parallel   enable GNU parallel
    --verbose    add verbosity


.. automodule:: prodigal
    :members:
    :show-inheritance:
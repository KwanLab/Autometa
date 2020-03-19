Bedtools
========

File containing wrapper functions for bedtools.

**Running the script as stand alone module**
While in the Autometa directory</path/to/output/file/prots.out> and with the conda environment activated, run: 
``python3 -m autometa.common.external.prodigal </path/to/BAM/alignment.bam> </path/to/genome/output/lengths.tsv> </path/to/output/alignment.bed> --[options]``

**Inputs** : ibam file

**Returns** : A tab-delimited table with index = contig and column = coverage

**Usage and Command line options:**

.. code-block:: bash

    python -m autometa.common.external.bedtools -h

    Usage: Retrieves genome coverage from the input ibam file [-h] [--coverage COVERAGE] 
                                                              [--force-bed] [--force-cov]
                                                              ibam lengths bed

    positional arguments:
    ibam                    </path/to/BAM/alignment.bam>
    lengths                 </path/to/genome/output/lengths.tsv> tab-delimited format
                            columns = [contig,length]
    bed                     </path/to/output/alignment.bed> tab-delimited format
                            columns = [contig,length]

    optional arguments:
    -h, --help              show this help message and exit
    --coverage COVERAGE     </path/to/coverage.tsv>
    --force-bed             force overwrite `bed` file
    --force-cov             force overwrite `--coverage` file 



.. automodule:: bedtools
    :members:
    :show-inheritance:
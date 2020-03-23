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

**Usage and Command line options:**

.. code-block:: bash
 
  Usage: Count the k-mers in the given sequence   [-h]  --fasta FASTA --kmers KMERS [--size SIZE]
                                                  [--normalized NORMALIZED] [--embedded EMBEDDED]
                                                  [--method {TSNE,UMAP}] [--n-components N_COMPONENTS]
                                                  [--do-pca] [--pca-dimensions PCA_DIMENSIONS]
                                                  [--multiprocess] [--nproc NPROC]

  Optional arguments:
    -h, --help                      show this help message and exit
    --fasta FASTA                   </path/to/input/sequences.fna>
    --kmers KMERS                   </path/to/output/kmers.tsv>  
                                    (will skip if file already exists)
    --size SIZE                     k-mer size
    --normalized NORMALIZED         </path/to/output/kmers.normalized.tsv> 
                                    (will skip if file already exists)
    --embedded EMBEDDED             </path/to/kmers.embedded.tsv>
    --method {TSNE,UMAP}            embedding method to use
    --n-components N_COMPONENTS     num components of lower dimension manifold
    --do-pca                        Whether to perform PCA prior to manifold learning
    --pca-dimensions PCA_DIMENSIONS number components to reduce to PCA feature space to
    --multiprocess                  count k-mers using multiprocessing
    --nproc NPROC                   number of processors to use if multiprocess is selected
                                    (default = 4)

.. automodule:: kmers
    :members:
    :show-inheritance:


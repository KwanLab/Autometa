kmers.py
=========


.. code-block:: python

    python -m autometa.common.kmers -h
    usage: Count k-mers [-h] --fasta FASTA --kmers KMERS [--size SIZE]
                        [--normalized NORMALIZED] [--embedded EMBEDDED]
                        [--method {TSNE,UMAP}] [--n-components N_COMPONENTS]
                        [--do-pca] [--pca-dimensions PCA_DIMENSIONS]
                        [--multiprocess] [--nproc NPROC]

    optional arguments:
      -h, --help            show this help message and exit
      --fasta FASTA         </path/to/sequences.fna>
      --kmers KMERS         </path/to/output/kmers.tsv> (will skip if file exists)
      --size SIZE           k-mer size
      --normalized NORMALIZED
                            </path/to/output/kmers.normalized.tsv> (will skip if
                            file exists)
      --embedded EMBEDDED   </path/to/kmers.embedded.tsv>
      --method {TSNE,UMAP}  embedding method
      --n-components N_COMPONENTS
                            <num components of lower dimension manifold>
      --do-pca              Whether to perform PCA prior to manifold learning
      --pca-dimensions PCA_DIMENSIONS
                            <num components to reduce to PCA feature space
      --multiprocess        count k-mers using multiprocessing
      --nproc NPROC         num. processors to use if multiprocess is selected.
                            (default = num. cpus available)

Will output tab-delimited matrix of k-mer frequency counts with an index column corresponding
to each contig in the input `--fasta`.

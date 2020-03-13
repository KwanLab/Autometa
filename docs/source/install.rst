Install
=========

Currently installation is supported by constructing a conda_ environment

conda environment creation

.. code-block:: bash

    conda create -n autometa "python>=3.7" --yes
    conda install -n autometa -c bioconda -c conda-forge --yes \
        biopython \
        pandas \
        tqdm \
        numpy \
        scikit-learn \
        samtools \
        bedtools \
        bowtie2 \
        hmmer \
        prodigal \
        diamond \
        ipython \
        ndcctools \
        parallel \
        requests \
        hdbscan \
        umap-learn \
        && conda clean --all --yes


.. _conda: https://docs.conda.io/en/latest/

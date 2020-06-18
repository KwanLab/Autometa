=======
Install
=======


Currently installation is supported by constructing a conda_ environment. You need to be running
a conda_ package (eg. miniconda, anaconda, etc) to install autometa.

Conda environment creation:

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
        ndcctools \
        parallel \
        requests \
        umap-learn \
        && conda clean --all --yes


.. _conda: https://docs.conda.io/en/latest/

You can now activate the conda environment using

    ``conda activate autometa``

Now download autometa to your desired directory using:

    ``git clone https://github.com/KwanLab/Autometa.git``

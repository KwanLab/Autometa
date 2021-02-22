=======
Install
=======


Currently installation is supported by constructing a conda_ environment. You need to be running
a conda_ package (eg. miniconda, anaconda, etc) to install autometa.

Conda environment creation:

.. code-block:: bash

    conda create -n autometa "python>=3.7"
    conda install -n autometa -c bioconda -c conda-forge \
        biopython \
        pandas \
        tqdm \
        tsne \
        numpy \
        scikit-learn \
        scikit-bio \
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
        hdbscan


.. _conda: https://docs.conda.io/en/latest/

You can now activate the conda environment using

    ``conda activate autometa``

Now download autometa to your desired directory using:

    ``git clone https://github.com/KwanLab/Autometa.git``


If you are wanting to help develop autometa, you will need these additional dependencies:

.. code-block:: bash

    conda install -n autometa \
        black pre_commit pytest pytest-cov pytest-html pytest-repeat pytest-variables

    # Navigate to your autometa conda environment
    conda activate autometa

    # If you'd like some easier visuals when running tests...
    pip install pytest-emoji pytest-md

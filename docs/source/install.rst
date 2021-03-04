=======
Install
=======


Currently installation is supported by constructing a conda_ environment. You need to be running
a conda_ package (eg. miniconda, anaconda, etc) to install autometa.

Conda installation
==================

1. Install miniconda_
1. Create a new environment ``conda create -n autometa "python>=3.7"``
1. Install autometa ``conda install -c conda-forge -c bioconda autometa --yes``
1. Actiavate autometa envoirnonment ``conda activate autometa``


.. _conda: https://docs.conda.io/en/latest/

.. _miniconda: https://docs.conda.io/en/latest/miniconda.html 


If you are wanting to help develop autometa, you will need these additional dependencies:

.. code-block:: bash

    conda install -n autometa \
        black pre_commit pytest pytest-cov pytest-html pytest-repeat pytest-variables

    # Navigate to your autometa conda environment
    conda activate autometa

    # If you'd like some easier visuals when running tests...
    pip install pytest-emoji pytest-md

=======
Install
=======


Currently installation is supported by constructing a conda_ environment. You need to be running
a conda_ package (eg. miniconda, anaconda, etc) to install autometa.

Conda installation
==================

#. Install miniconda_
#. Create a new environment ``conda create -n autometa "python>=3.7"``
#. Install autometa ``conda install -c conda-forge -c bioconda autometa --yes``
#. Actiavate autometa envoirnonment ``conda activate autometa``

Docker image
============

You can also run Autometa using a prebuild Docker image. 

#. Install Docker_
#. Run the following commands

    .. code-block:: bash

        git clone https://github.com/KwanLab/Autometa.git
        docker pull jasonkwan/autometa:latest

Testing Autometa
================

Dependencies
------------

.. code-block:: shell

    # Check dependencies to see if Autometa environment is appropriately configured
    python autometa.py --check-dependencies

    # If any of the checks return False, you can check which failed using
    python autometa.py --check-dependencies --debug


Once you have resolved the dependencies your checks should be returning ``True``.

Unit tests
----------

You can also check the installation using autometa's built in unit tests.This is not at all necessary and is primarily meant for development and debuggging purposes. To run the tests however, you'll first need to install the following packages and download the test dataset.

.. code-block:: shell

    # Install packages for testing
    conda install -n autometa -c conda-forge \
        black pre_commit pytest pytest-cov pytest-html pytest-repeat pytest-variables gdown --yes

    # Download test data
    gdown https://drive.google.com/uc\?\id=1bSlPldaq3C6Cf9Y5Rm7iwtUDcjxAaeEk -O tests/data/test_data.json

Run the tests using ``make test``. This will run the tests and make sure that everything is working well. 

Additional unit tests are provided in the test directory. These are designed to aid in future development of autometa.

.. _conda: https://docs.conda.io/en/latest/
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html 
.. _Docker: https://www.docker.com/


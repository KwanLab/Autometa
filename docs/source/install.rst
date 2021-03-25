=======
Install
=======

Currently installation is supported by constructing a conda_ environment. You need to be running
a conda_ package (eg. miniconda, anaconda, etc) to install autometa. Before you proceed to installation download and install a conda_ package. You can start with miniconda_.

Direct installation (Quickest)
==============================

#. Install miniconda_
#. Create a new environment and install autometa using ``conda create -c bioconda -n autometa autometa``
#. Actiavate autometa envoirnonment ``conda activate autometa``

Install from source (using make)
================================

Download and install miniconda_. Now run the following commands:

.. code-block:: bash

    # Navigate to the directory where you need to clone Autometa
    cd $HOME
    # Clone the autimeta git repository
    git clone https://github.com/KwanLab/Autometa.git
    # Navigate into the cloned repository
    cd $HOME/Autometa
    # List all make commands
    make
    # create autometa conda environment
    make create_environment
    # activate autometa conda environment
    conda activate autometa
    # install autometa source code
    make install

Install from source (full commands)
===================================

Download and install miniconda_. Now run the following commands:

.. code-block:: bash

    # Navigate to the directory where you need to clone Autometa
    cd $HOME
    # Clone the autimeta git repository
    git clone https://github.com/KwanLab/Autometa.git
    # Navigate into the cloned repository
    cd $HOME/Autometa
    # Construct the autometa environment from requirements.txt
    conda create -n autometa --file=requirements.txt
    # Install the autometa code base from source
    python setup.py install

Docker image
============

You can build a docker image for your clone of Autometa repository. 

#. Install Docker_
#. Run the following commands

.. code-block:: bash

    # Navigate to the directory where you need to clone Autometa
    cd $HOME
    # Clone the autimeta git repository
    git clone https://github.com/KwanLab/Autometa.git
    # Navigate into the cloned repository
    cd $HOME/Autometa
    # This will tag the image as jason-c-kwan/autometa:<your current branch>
    make image

Testing Autometa
================

You can also check the installation using autometa's built in unit tests.This is not at all necessary and is primarily meant for development and debuggging purposes. To run the tests however, you'll first need to install the following packages and download the test dataset.

.. code-block:: bash

    # Activate your autometa conda environment
    conda activate autometa
    # List all make options
    make
    # Install dependencies for test environment
    make test_environment
    # Download test_data.json for unit testing to tests/data/
    make unit_test_data_download

You can now run different unit tests using the following commands:

.. code-block:: bash

    # Run all unit tests
    make unit_test
    # Run unit tests marked with entrypoint
    make unit_test_entrypoints
    # Run unit tests marked with WIP
    make unit_test_wip

For more information about the above commands see the :ref:`Contributing Guidelines` page. Additional unit tests are provided in the test directory. These are designed to aid in future development of autometa. 

.. _conda: https://docs.conda.io/en/latest/
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html 
.. _Docker: https://www.docker.com/
.. _anaconda: https://www.anaconda.com/
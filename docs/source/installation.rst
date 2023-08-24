.. _installation-page:

============
Installation
============

Currently Autometa package installation is supported by mamba_, and docker_.
For installation using mamba, download mamba from Mambaforge_.

.. attention::

    If you are only trying to run the Autometa workflow, you should start at :ref:`Getting Started` before proceeding.

Direct installation (Quickest)
==============================

#. Install mamba_

    .. code-block:: bash

        wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
        bash Mambaforge-$(uname)-$(uname -m).sh

    Follow the installation prompts and when you get to this:

    .. code-block:: bash

        Do you wish the installer to initialize Mambaforge
        by running conda init? [yes|no]
        [no] >>> yes

    This will require restarting the terminal, or resetting
    the terminal with the source command

    .. code-block:: bash

        # To resolve the comment:
        # ==> For changes to take effect, close and re-open your current shell. <==
        # type:
        source ~/.bashrc

    .. note::

        If you already have conda installed, you can install mamba as a drop-in replacement.

        .. code-block:: bash

            conda -n base -c conda-forge mamba -y


#. Create a new environment with ``autometa`` installed:

    .. code-block:: bash

        mamba create -c conda-forge -c bioconda -n autometa autometa

    .. note::

            You may add the ``bioconda`` and ``conda-forge`` channels to your mamba
            config to simplify the command.

            .. code-block:: bash

                mamba config --append channels bioconda
                mamba config --append channels conda-forge

            Now mamba will search the ``bioconda`` and ``conda-forge``
            channels alongside the defaults channel.

            .. code-block:: bash

                mamba create -n autometa autometa


#. Activate ``autometa`` environment:

    .. code-block::

        mamba activate autometa

Install from source (using make)
================================

Download and install mamba_. Now run the following commands:

.. code-block:: bash

    # Navigate to the directory where you would like to clone Autometa
    cd $HOME

    # Clone the Autometa repository
    git clone https://github.com/KwanLab/Autometa.git

    # Navigate into the cloned repository
    cd Autometa

    # create autometa mamba environment
    make create_environment

    # activate autometa mamba environment
    mamba activate autometa

    # install autometa source code in autometa environment
    make install

.. note::

    You can see a list of all available make commands by running ``make`` without any other arguments.

Install from source (full commands)
===================================

Download and install mamba_. Now run the following commands:

.. code-block:: bash

    # Navigate to the directory where you would like to clone Autometa
    cd $HOME

    # Clone the Autometa repository
    git clone https://github.com/KwanLab/Autometa.git

    # Navigate into the cloned repository
    cd Autometa

    # Construct the autometa environment from autometa-env.yml
    mamba env create --file=autometa-env.yml

    # Activate environment
    mamba activate autometa

    # Install the autometa code base from source
    python -m pip install . --ignore-installed --no-deps -vv

Building the Docker image
=========================

You can build a docker image for your clone of the Autometa repository.

#. Install Docker_
#. Run the following commands

.. code-block:: bash

    # Navigate to the directory where you need to clone Autometa
    cd $HOME

    # Clone the Autometa repository
    git clone https://github.com/KwanLab/Autometa.git

    # Navigate into the cloned repository
    cd Autometa

    # This will tag the image as jasonkwan/autometa:<your current branch>
    make image

    # (or the full command from within the Autometa repo)
    docker build . -t jasonkwan/autometa:`git branch --show-current`

Testing Autometa
================

You can also check the installation using autometa's built-in unit tests.
This is not at all necessary and is primarily meant for development and debugging purposes.
To run the tests, however, you'll first need to install the following packages and download the test dataset.

.. code-block:: bash

    # Activate your autometa mamba environment
    mamba activate autometa

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

.. note::

    As a shortcut you can also create the test environment and run **all** the unit tests using ``make unit_test`` command.

For more information about the above commands see the :ref:`Contributing Guidelines` page.
Additional unit tests are provided in the test directory. These are designed to aid in future development of autometa.

.. _mamba: https://mamba.readthedocs.io/en/latest/index.html
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _Docker: https://www.docker.com/

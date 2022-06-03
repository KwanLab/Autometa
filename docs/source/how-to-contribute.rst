=================
How to contribute
=================

Autometa is an open-source project developed on GitHub. If you would like to help develop
Autometa or have ideas for new features please see our `contributing guidelines <https://github.com/KwanLab/Autometa/blob/dev/.github/CONTRIBUTING.md>`__

Some good first issues are available on the KwanLab Autometa GitHub repository `good first issues 🥇💡 <https://github.com/KwanLab/Autometa/contribute>`__

If you are wanting to help develop Autometa, you will need these additional dependencies:

Documentation
=============

Autometa builds documentation using `readthedocs <https://readthedocs.org/>`__. You have to install the following to be able to build the docs

.. code-block:: bash

    # Activate your autometa conda environment
    conda activate autometa
    # Install dependencies
    conda install -n autometa -c conda-forge \
        sphinx sphinx_rtd_theme
    # List all make options
    make
    # Build documentation for autometa.readthedocs.io
    make docs

``make docs``
-------------

This command runs sphinx and generates autometa documentation for autometa.readthedocs.io.

Unit tests
==========

You will have to install certain dependencies as well as test data to be able to run and develop unit tests.

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

``make test_environment``
-------------------------

This command installs all the dependencies that you need to successfully run the unit tests.

``make unit_test_data_download``
--------------------------------

This command downloads the test_data.json object that you need to run the Unit tests. This is a necessary step when wanting to run unit tests as the test_data.json file will hold many of the variables necessary to conduct these tests.

``make unit_test_data_build``
-----------------------------

This is used to create your own ``test_data.json`` object locally. This step is NOT required for running unit tests, you can directly download the ``test_data.json`` object using the previous command. This command is needed in case you are changing file formats or adding more objects into the test suite. To do this you first need to download all the files from `here <https://drive.google.com/open?id=189C6do0Xw-X813gspsafR9r8m-YfbhTS>`__ in ``tests/data/`` and then run ``make unit_test_data_build``. This would generate a similar ``test_data.json`` object that you get by running the previous command.

The above command is used to manually build the test_data.json file for unit testing. I.e. it will run the script make_test_data.py which will aggregate all of the files in the tests/data folder that have been downloaded from `here <https://drive.google.com/open?id=189C6do0Xw-X813gspsafR9r8m-YfbhTS>`__. This is the first or perhaps 0th step when it comes to running the tests without an already generated ``test_data.json`` object as it generates the test_data.json file that is parsed to retrieve all of the pre-generated variables used for intermediate stages of the pipeline. This is done to reduce the test time and computational workload when running through the test suite.

``make unit_test``
------------------

This command runs all unit tests under the tests directory. This includes all tests marked as WIP or as entrypoints. However this will skip tests marked with the following decorator:

.. code-block:: python

    @pytest.mark.skip
    def test_some_function(...):
        ...

``make unit_test_entrypoints``
------------------------------

This command runs the tests marked as entrypoints. This is denoted in pytest with the decorator:

.. code-block:: python

    @pytest.mark.entrypoint
    def test_some_function_that_is_an_entrypoint(...):
    ...

Entrypoints correspond to the entry point functions listed out by 'console scripts' in setup.py. These entry point functions are aliased to provide more intuitive commands for the end user. These are important and sometimes referred to as "happy" tests because if one of these fail for the end-user, they will probably be quite unhappy and likely distrust the functionality of the rest of the codebase.

``make unit_test_wip``
----------------------

This command runs the tests marked as work-in-progress (WIP). This is denoted in pytest with the decorator:

.. code-block:: python

    @pytest.mark.wip
    def test_some_function_that_is_wip(...):
    ...

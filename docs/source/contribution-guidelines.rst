=======================
Contributing Guidelines
=======================

"Autometa is an open-source project developed on
GitHub. If you would like to help develop Autometa or
have ideas for new features please see our `contributing guidelines <https://github.com/KwanLab/Autometa/blob/master/.github/CONTRIBUTING.md>`__

For developers
==============

If you are wanting to help develop autometa, you will need these additional dependencies:

.. code-block:: bash

    conda install -n autometa -c conda-forge\
        black pre_commit pytest pytest-cov pytest-html pytest-repeat pytest-variables sphinx sphinx_rtd_theme

    # Navigate to your autometa conda environment
    conda activate autometa

    # If you'd like some easier visuals when running tests...
    pip install pytest-emoji pytest-md
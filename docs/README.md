Autometa documentation
========================

Building the documentation
--------------------------

This file will detail how to download the docs and contribute to them.
What files to install and how to run

Installation
------------

Each script in the source will be run with the help flag to automatically place,
usage information for the respective script. Therefore, the packages necessary to
run the script need to be available. In order to perform the build, the `autometa`
conda environment needs to be active and additional packages need to be installed.

```bash
# Note: This assumes you've already created the autometa environment.
# Install additional dependencies to autometa environment...
conda install -n autometa -c conda-forge sphinxcontrib-programoutput sphinx sphinx_rtd_theme
```

Running make
------------

You can build your own documentation for inspection (and to check for any
formatting issues in your newly written docstrings) using the command:
 `make clean html` within the docs directory. If you would like to build the docs
from the base Autometa directory, you can run the command: `make clean html -C docs`

```shell
# Note: You will need to have already installed dependencies described above.
# First activate env.
conda activate autometa

# Then build the docs either from base Autometa directory.
make clean html -C docs && open docs/build/html/index.html
# clean will remove any files from previous build
# -C flag will first navigate to specified directory and look for Makefile there.

# Or build from `docs` directory.
cd docs && make clean html && open build/html/index.html
```


Organization
------------

:TODO: What files are here and description of each.

Development
-----------

:TODO: Contributing to the documentation

* Read old issues before submitting new ones
* Guidelines
* Fork the repo and send in pull requests

References
----------

* [Matplotlib documentation](https://github.com/matplotlib/matplotlib/blob/master/doc/devel/documenting_mpl.rst)
* [Writing good commit messages](https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)

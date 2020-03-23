# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
from sphinx.ext.autodoc import between


sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('../../autometa/common'))
sys.path.insert(0, os.path.abspath('../../autometa/config'))
sys.path.insert(0, os.path.abspath('../../autometa/binning'))
sys.path.insert(0, os.path.abspath('../../autometa/datasets'))
sys.path.insert(0, os.path.abspath('../../autometa/validation'))
sys.path.insert(0, os.path.abspath('../../autometa/common/external'))

""" sys.path.insert(0, sys.path.append('../..'))
sys.path.insert(0, sys.path.append('../autometa'))
sys.path.insert(0, sys.path.append('../autometa/common'))
sys.path.insert(0, sys.path.append('../autometa/config'))
sys.path.insert(0, sys.path.append('../autometa/binning'))
sys.path.insert(0, sys.path.append('../autometa/datasets'))
sys.path.insert(0, sys.path.append('../autometa/validation'))
sys.path.insert(0, sys.path.append('../autometa/common/external'))
sys.path.insert(0, os.path.abspath('../../autometa/common/external'))
 """
# -- Project information -----------------------------------------------------

project = 'Autometa'
copyright = '2020, Evan R. Rees, Shaurya Chanana, Siddharth Uppal, Kyle Wolf, Jason C. Kwan'
author = 'Evan R. Rees, Shaurya Chanana, Siddharth Uppal, Kyle Wolf, Jason C. Kwan'

# The full version, including alpha/beta/rc tags
release = '2.0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage',
'sphinx_rtd_theme', 'sphinx.ext.napoleon'] 

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_static_path = []


#Function to not include the copyright in the autodic extension
def setup(app):
    # Register a sphinx.ext.autodoc.between listener to ignore everything
    # between lines that contain the word COPYRIGHT
    app.connect('autodoc-process-docstring', between('^.*COPYRIGHT.*$', exclude=True))
    return app
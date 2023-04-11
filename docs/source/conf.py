# Autometa documentation build configuration file.
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import os
import subprocess
import sys

from sphinx.ext.autodoc import between
from datetime import datetime


sys.path.append(os.path.abspath("./_ext"))

for dirpath, dirnames, filenames in os.walk("../../", topdown=True):
    sys.path.insert(0, os.path.abspath(dirpath))

autodoc_mock_imports = [
    "Bio",
    "hdbscan",
    "tsne",
    "sklearn",
    "umap",
    "tqdm",
    "pandas",
    "numpy",
    "attr",
]

# fmt: off
import parse_argparse

# -- Project information -----------------------------------------------------

project = "Autometa"
copyright = f"2016 - {datetime.now().year}, Ian J. Miller, Evan R. Rees, Izaak Miller, Shaurya Chanana, Siddharth Uppal, Kyle Wolf, Jason C. Kwan"
author = "Ian J. Miller, Evan R. Rees, Izaak Miller, Shaurya Chanana, Siddharth Uppal, Kyle Wolf, Jason C. Kwan"

# The short X.Y version.
version = "2.0"

# The full version, including alpha/beta/rc tags
release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.todo",
    "sphinx_rtd_theme",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosectionlabel",
]

todo_include_todos = True
# Includes doctrings of functions that begin with double underscore
napoleon_include_special_with_doc = True
autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The suffix of source filenames.
source_suffix = ".rst"

# The encoding of source files.
source_encoding = "utf-8-sig"

# The master toctree document.
master_doc = "index"

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

highlight_language = "shell"

# -- Options for HTML output -------------------------------------------------
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

if not os.path.isdir("_static"):
    os.mkdir("_static")

html_static_path = ["_static"]

# If not None, a 'Last updated on:' timestamp is inserted at every page
# bottom, using the given strftime format.
# The empty string is equivalent to '%b %d, %Y'.

html_last_updated_fmt = "%b %d, %Y"


def setup(app):
    # Register a sphinx.ext.autodoc.between listener to ignore everything
    # between lines that contain the word COPYRIGHT
    # Exclude COPYRIGHT block in scripts using the autodoc between function
    app.connect("autodoc-process-docstring", between("^.*COPYRIGHT.*$", exclude=True))
    return app


def sphinx_apidoc(header="Autometa modules"):
    """
    Runs `spinx-apidoc`
    Example usage: sphinx-apidoc[OPTIONS] - o < OUTPUT_PATH > <MODULE_PATH > [EXCLUDE_PATTERN, ...]

    Parameters
    ----------
    header : str, optional
        Project name to be displayed on its tab in build docs, by default "Autometa modules"
    """

    source_dir = os.path.dirname(__file__)
    pkg_dir = os.path.join(os.path.dirname(os.path.dirname(source_dir)), "autometa")
    dest_dirname = header.replace(" ", "_")
    dest_dir = os.path.join(source_dir, dest_dirname)
    cmd = f"sphinx-apidoc --force --output-dir {dest_dir} {pkg_dir} --doc-project '{header}'"
    subprocess.call(cmd, shell=True)


parse_argparse.main()
sphinx_apidoc()

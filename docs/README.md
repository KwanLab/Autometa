Autometa documentation
========================

Building the documentation
--------------------------

This file will detail how to download the docs and contribute to them.
What files to install and how to run

Installation
------------

```bash
conda install -c anaconda sphinx_rtd_theme
```

Organization
------------

This is the top level build directory for the Autometa documentation.  All of the documentation is written using sphinx, a python documentation system built on top of ReST. This docs folder contains:

* api - placeholders to automatically generate the api documentation

* index.rst - the top level include document for Autometa docs

* conf.py - the sphinx configuration

* Makefile and make.bat - entry points for building the docs

* _static - used by the sphinx build system

* _templates - used by the sphinx build system

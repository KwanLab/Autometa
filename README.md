Autometa
=========

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/KwanLab/Autometa?label=Latest%20Version&sort=semver&style=flat-square)

An automated binning pipeline for single metagenomes, in particular host-associated and highly complex ones. Autometa is copyright 2020 Ian J. Miller, Evan Rees, Izaak Miller and Jason C. Kwan, and is released under the GNU Affero General Public License v3 (see LICENSE.txt). If you find Autometa useful to your work, please cite:

Miller, I. J.; Rees, E. R.; Ross, J.; Miller, I.; Baxa, J.; Lopera, J.; Kerby, R. L.; Rey, F. E.; Kwan, J. C. Autometa: Automated extraction of microbial genomes from individual shotgun metagenomes. *Nucleic Acids Research*, **2019**. [DOI: https://doi.org/10.1093/nar/gkz148](https://doi.org/10.1093/nar/gkz148)

Documentation
-------------

<a href="https://autometa.readthedocs.io"><img src="https://img.shields.io/badge/readthedocs-Autometa-blue"/></a>
[![Documentation Status](https://readthedocs.org/projects/autometa/badge/?version=latest)](https://autometa.readthedocs.io/en/latest/?badge=latest)
<a href="https://github.com/KwanLab/Autometa/discussions"><img src="https://img.shields.io/github/discussions/KwanLab/Autometa"/></a>

Quickstart
-------------------------------------------------------------------------------------------------------------------------

### :shell: Bash workflow and :snake: Autometa package

[![Install with Conda Badge](https://anaconda.org/bioconda/autometa/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Platforms Badge](https://anaconda.org/bioconda/autometa/badges/platforms.svg)](https://anaconda.org/bioconda/autometa)
[![Downloads Badge](https://anaconda.org/bioconda/autometa/badges/downloads.svg)](https://anaconda.org/bioconda/autometa)

#### Install into current env

```bash
conda install -c bioconda autometa
```

#### Create new env

```bash
conda create -n autometa -c bioconda autometa
```

### :whale: Run with Docker

<a href="https://hub.docker.com/r/jasonkwan/autometa"><img src="https://img.shields.io/docker/image-size/jasonkwan/autometa/main?style=flat-square" alt="Docker Image Size (tag)"/></a>

```bash
docker run -it --rm jasonkwan/autometa:main
```

### :green_apple: Nextflow Workflow

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000?style=flat-square)](https://www.nextflow.io/)

#### Install into current env

```bash
conda env update -n <your-env> --file=https://raw.githubusercontent.com/KwanLab/Autometa/main/nextflow-env.yml
```

#### Create new env

```bash
conda env create --file=https://raw.githubusercontent.com/KwanLab/Autometa/main/nextflow-env.yml
```

For developers
--------------

<a href="https://github.com/KwanLab/Autometa/blob/main/.github/CONTRIBUTING.md"><img src="https://img.shields.io/badge/Contributing-Guidelines-blue" alt="Contributing Guidelines"/></a>
<a href="https://github.com/KwanLab/Autometa/contribute"><img alt="GitHub good first issues" src="https://img.shields.io/github/issues-search/KwanLab/Autometa?color=purple&label=Good%20First%20Issues&query=label%3A%22good%20first%20issue%22"></a>

![GitHub language count](https://img.shields.io/github/languages/count/KwanLab/Autometa)
![GitHub top language](https://img.shields.io/github/languages/top/KwanLab/Autometa)
[![codecov](https://codecov.io/gh/KwanLab/Autometa/branch/main/graph/badge.svg?token=N2X4F6P6M5)](https://codecov.io/gh/KwanLab/Autometa)

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
In addition, references of tools and data used in this pipeline are as follows:

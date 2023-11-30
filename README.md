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

Full documentation is hosted on [autometa.readthedocs.io](https://autometa.readthedocs.io "Autometa documentation")

Quickstart
-------------------------------------------------------------------------------------------------------------------------

### :shell: Bash workflow

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/autometa/README.html) [![Conda](https://img.shields.io/conda/dn/bioconda/autometa.svg)](https://anaconda.org/bioconda/autometa/files)

#### 1. Setup env

##### Install into your current env...

```bash
mamba install -c conda-forge -c bioconda autometa
```

##### ... or create a new env

```bash
mamba create -n autometa -c conda-forge -c bioconda autometa
```

#### 2. Download the bash workflow template, [autometa.sh](https://raw.githubusercontent.com/KwanLab/Autometa/main/workflows/autometa.sh "autometa.sh template")

#### 3. Edit the input parameters

#### 4. Run the workflow

### :green_apple: Nextflow Workflow

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000?style=flat-square)](https://www.nextflow.io/)

#### 1. Setup env

> :student: This workflow requires only nextflow and nf-core be installed.

##### Install into your current env...

```bash
mamba env update -n <your-env> --file=https://raw.githubusercontent.com/KwanLab/Autometa/main/nextflow-env.yml
```

##### ... or create a new env

```bash
mamba env create --file=https://raw.githubusercontent.com/KwanLab/Autometa/main/nextflow-env.yml
# Activate the env after creation
mamba activate autometa-nf
```

#### 2. Launch and run the workflow

```bash
nf-core launch KwanLab/Autometa
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

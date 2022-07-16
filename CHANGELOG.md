# Autometa: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.0] - 2022-05-11

### `Added`

- Addition of `autometa` entrypoint to retrieve info on citation, version and other commands
- :green_apple: Addition of autometa version written within each autometa-related nextflow process
- :shell: The autometa version is emitted at the beginning of the autometa bash workflow.

### `Fixed`

- :shell::bug: Renamed incorrect variables `$coverage` to `$coverages` and `$clustering_method` to `$cluster_method` in `workflows/autometa.sh`
- :snake::bug: Fix use of `kwargs` in kmers `embed(..., *kwargs)` function
- :snake: Fix passing `n_jobs` to kmers `embed(..., n_jobs=1)` function

### `Dependencies`

- :art: Pinned black version to 22.3.0
- ⬆️:green_heart: upgrade diamond to version 2.0.14

## [2.0.3] - 2022-04-07

### What's Changed

- :memo: Add CITATION file by @shaneroesemann in <https://github.com/KwanLab/Autometa/pull/242>
- :bug::art::green_apple: Allow utf-8-sig sample sheets #250
- :bug::memo::green_apple: Replace incorrect runtime defaults specified in nf-core launch command #252
- :art::green_apple: Continue nextflow workflow (without termination) if bins are not recovered from a particular dataset #253
- :art::bug: `autometa-benchmark` entrypoint, when checking taxid classifications, now converts any taxids of `0` to `1` #254
- :bug: TAXA_DB connection is checked prior to TAXA_DB database files download rather than pinging google DNS #258
- :bug: Fix `autometa-length-filter` bug where specifying a directory when writing to `--output-fasta` was required (now optional) #256
- 💚 ⬆️:whale::green_apple:  Bump version to 2.0.3 in `VERSION` and `manifest.version` in `nextflow.config` s.t. autometa docker images used in nextflow workflow correspond to most recent release

### New Contributors

- @shaneroesemann made their first contribution in <https://github.com/KwanLab/Autometa/pull/242>

**Full Changelog**: <https://github.com/KwanLab/Autometa/compare/2.0.2...2.0.3>

## [2.0.2] - 2022-03-15

### What's Changed

- :bug: Fix the NoneType error during initating the LCA object (#246) by @chtsai0105 in <https://github.com/KwanLab/Autometa/pull/247>

### New Contributors

- @chtsai0105 made their first contribution in <https://github.com/KwanLab/Autometa/pull/247>

**Full Changelog**: <https://github.com/KwanLab/Autometa/compare/2.0.1...2.0.2>

## [2.0.1] - 2022-02-24

### What's Changed

- Environment and directory structure updates by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/8>

- resolved #10 Contributors added and copyright year updated to 2020. by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/15>
- Resolving issues #16, #17, #18, #21  and update to Autometa API and Logger by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/25>
- Resolved #19 added docstring, fixed nproc and removed depth function by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/29>
- Documentation by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/34>
- Issue #5 Working conda recipe by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/38>
- 🐛found in coverage by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/49>
- fixes #2  by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/47>
- Contributing Guidelines by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/50>
- Add Markers class documentation. by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/62>
- Add **main**.py documentation (fixes #60) by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/63>
- Add functionality to bin without taxonomy. Update docstrings by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/65>
- Documentation by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/45>
- Remove merge conflict resolution lines (Fixes #68) by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/69>
- :art::bug: Add mock import of modules and link to contribution guidelines (fixes #22) by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/70>
- Resolves #55 Environ by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/76>
- fixes-#54 Metagenome by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/66>
- Update MAG class to MetaBin by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/67>
- hmmer by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/72>
- LCA by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/78>
- Fix writing by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/82>
- Update majority_vote by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/81>
- pre-commit hooks by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/92>
- verbose bug by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/90>
- Recursive DBSCAN by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/84>
- databases and utilities by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/77>
- Rank-specific binning by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/96>
- diamond.py by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/87>
- Add support request issue template. by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/97>
- ncbi.py by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/83>
- Samtools by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/103>
- Binning stats/taxonomy summary by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/99>
- decision tree classifier by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/100>
- Fix config and setup of user project by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/104>
- Update project docstrings by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/108>
- CI/CD by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/101>
- :bug: Change > to >= when calculating N50 by @chasemc in <https://github.com/KwanLab/Autometa/pull/119>
- Fix Dockerfile by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/123>
- Add support for gzipped assemblies by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/129>
- Update bug report template by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/130>
- Remove --multiprocess from autometa-kmers entrypoint by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/127>
- Add GC content std.dev. limit and coverage std. dev. limit Binning metrics by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/120>
- Nextflow implementation template by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/118>
- Update documentation by @Sidduppal in <https://github.com/KwanLab/Autometa/pull/121>
- Add feature to download google drive datasets by @ajlail98 in <https://github.com/KwanLab/Autometa/pull/138>
- Add densmap embed method and fix binning-summary cluster column bug by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/176>
- Classification and Clustering Benchmarking by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/141>
- Nfcore and structuring modules for collaboration by @chasemc in <https://github.com/KwanLab/Autometa/pull/157>
- Delete .gitattribute - there is a .gitattributes by @chasemc in <https://github.com/KwanLab/Autometa/pull/190>
- :fire::green_apple: Remove duplicate standard slurm profiles by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/195>
- Fix import error in databases.py by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/194>
- Fix/Create mock data subworkflow by @chasemc in <https://github.com/KwanLab/Autometa/pull/206>
- 🐳:bug: Docker fix :whale: by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/213>
- :memo: Update Documentation by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/212>
- :art: Add typehints and update kmers docstring by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/208>
- :art::snake: Add specific parsers for domtblout and tblout for hmmscan output formats by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/201>
- :snake::art::bug: Update metagenome.length_filter(...) by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/210>
- Fix bedtools genomecov deprecation (coverage calculation) by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/209>
- :art:🐚  Add bash-implementations of Autometa workflows by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/202>
- Nextflow documentation by @chasemc in <https://github.com/KwanLab/Autometa/pull/184>
- :fire::memo: Reformat benchmarking docs by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/215>
- Refactor autometa-taxonomy-lca by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/211>
- :art::whale: Replace jason-c-kwan with jasonkwan for docker images by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/217>
- Update mock_data_reporter.Dockerfile by @chasemc in <https://github.com/KwanLab/Autometa/pull/220>
- :whale::green_heart::memo: Add docker CI and update links by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/216>
- Simplify Licensing by @chasemc in <https://github.com/KwanLab/Autometa/pull/222>
- Add check for nr.dmnd and nr.gz by @chasemc in <https://github.com/KwanLab/Autometa/pull/221>
- 🍏 Change Nextflow I/O behavior by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/218>
- :green_apple::art: add/update coverage handling by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/223>
- :snake:🐎 Large data mode by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/207>
- Merge main into dev by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/224>
- pytest & codecov CI by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/227>
- Refactor samplesheet by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/228>
- 🐛 🎨 🍏  Fix kingdom-handling and mounting TAXA_DB databases into docker container by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/229>
- Add error handling strategies for nextflow processes by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/231>
- Release 2.0.0 by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/233>
- :memo::art::arrow_up::fire: update build files and respective docs by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/234>
- Add badges and links to README.md by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/235>
- GH action: Add dynamic docker tags by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/236>
- Update nextflow-workflow.rst by @samche42 in <https://github.com/KwanLab/Autometa/pull/238>
- :bug::fire: Fix bedtools coverage calculation bug by @WiscEvan in <https://github.com/KwanLab/Autometa/pull/243>

### New Contributors

- @chasemc made their first contribution in <https://github.com/KwanLab/Autometa/pull/119>

- @samche42 made their first contribution in <https://github.com/KwanLab/Autometa/pull/238>

**Full Changelog**: <https://github.com/KwanLab/Autometa/compare/1.0.3...2.0.1>

## [2.0.0] - 2022-02-09

Second release of Autometa

### `Added`

- 🍏 Created autometa nextflow workflow with the [nf-core](https://nf-co.re/) template.
- :shell: Created autometa bash workflow
- Refactored entire codebase for easy installation using pip and/or bioconda
- Addition of various k-mer embedding parameters:
    - k-mer size
    - normalization method
    - PCA dimensions
    - embedding method
- Addition of HDBSCAN clustering method
- Addition of autometa entrypoints
- Easy install with `make` commands in `Makefile`
- Easy installation using `pip`
- Easy installation using `conda install -c bioconda autometa`
- Autometa library of functions available within the python environment
- Granularity between autometa tasks to reduce re-computations on failed tasks
- Configuration to specify paths to required databases
- Added tests for CI/CD
- :whale: Update Dockerfile to conform to new autometa environment & commands
- :memo: Addition of documentation hosted on [autometa.readthedocs.io](https://autometa.readthedocs.io "Autometa documentation")
- Dependencies specified in `autometa-env.yml`
- Dependencies specified in `requirements.txt`

### `Fixed`

- :art::snake: Restructured source code directory to allow installation of autometa console scripts
- :art::snake: Restructured source code scripts to prevent calling code on import

### `Removed`

- python 2.7
- *All* autometa version `1.*` commands

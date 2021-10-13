.config files are prefixed by numbers (e.g. "01_optional_autometa.config")

This is the order in which the files are loaded from `~/Autometa/nextflow.config.

If a parameters is defined more than once, the one loaded last will be used.


Four config files are not numbered and are nf-core based files

  - `base.config`
  - `modules.config`
  - `test_full.config`
  - `test.config`

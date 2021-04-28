Work in progress...


modules under `/nextflow/modules/core` are required for autometa v2.0

Modules providing new functionality should be created in a new directory with a short, descriptive name, within `/nextflow/modules`

`process{}` statements should be placed in a subfolder named `process`
workflow statemnt(s) go in the top-level of your module's directory
For a good example of a subworkflow with processes, see: `/nextflow/modules/core/taxonomy`


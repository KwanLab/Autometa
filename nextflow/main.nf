#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

  params.config_profile_description = false 
  params.config_profile_contact = false 
  params.config_profile_url = false 
/*
========================================================================================
                         Autometa   
========================================================================================
 Autometa's Nextflow Analysis Pipeline
 #### Homepage
 https://github.com/KwanLab/Autometa
 #### Documentation
 https://autometa.readthedocs.io/en/latest/
----------------------------------------------------------------------------------------
*/

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run autometa --input 'input/path/sample.fna'" //TODO
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

// Path of fna files used as input:
if ( !params.input || params.input instanceof Boolean )
error """
You must supply the `input` parameter in the config or on the command line!
For a single input use the path
    nextflow run main.nf -c parameters.config --input "</path/to/your/metagenome.fna>"
For multiple input files, glob patterns may be used:
    nextflow run main.nf -c parameters.config --input "</path/to/your/metagenome_*.fna>"
"""

// Where to store final results:
if ( !params.outdir || params.outdir instanceof Boolean )
error """
You must supply the `--outdir` parameter in the config or on the command line!
e.g.
nextflow run main.nf -c parameters.config --outdir "</directory/path/to/final/results>"
"""

// Where to store intermediate results:
if ( !params.interim_dir )
error """
You must supply the `--interim_dir` parameter in the config or on the command line!
e.g.
nextflow run main.nf -c parameters.config --interim_dir "</directory/path/to/store/interimediate/results>"
"""
if ( params.interim_dir = true ) {
    params.interim_dir = "${params.outdir}/interim_outputs"
    println "Intermediate results will be saved to:\n${params.outdir}/interim_outputs"
}
if ( params.interim_dir = false ) {
    // TODO: Don't publish/create the copies of intermediate results
}


////////////////////////////////////////////////////
/* --         STAGE REPORTING CONFIG           -- */
////////////////////////////////////////////////////
// Stage config files for documentation/reporting
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)


////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
// TODO nf-core: Report custom parameters here
summary['Input']            = params.input
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
////////////////////////////////////////////////////
/* --         COLLECT SUMMARY FOR LOG          -- */
////////////////////////////////////////////////////
Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'autometa-nextflow-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'Autometa Workflow Summary'
    section_href: 'https://github.com/KwanLab/Autometa'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

////////////////////////////////////////////////////
/* --         COLLECT SOFTWARE VERSIONS        -- */
////////////////////////////////////////////////////
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml', emit: ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    python --version > v_python.txt
    docker --version > v_docker.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}
/*
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: "${params.publish_dir_mode}"

    input:
    path output_docs 

    output:
    path 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
 */

//log.info summary

/*
 * Retrieve the main 'AUTOMETA' worflow
 */
include { AUTOMETA } from './modules/core/autometa.nf' addParams(single_db_dir: params.single_db_dir)
// include { DOWNLOAD_DATABASES } from './modules/core/download-databases.nf' addParams(single_db_dir: params.single_db_dir)

workflow {
  Channel
    .fromPath(params.input, checkIfExists: true, type: 'file')
    .set{unfiltered_metagenome_ch}
  //DOWNLOAD_DATABASES()
  AUTOMETA(unfiltered_metagenome_ch)
}


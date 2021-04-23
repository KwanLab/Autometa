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

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-example-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/example Workflow Summary'
    section_href: 'https://github.com/nf-core/example'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

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


log.info """

 Autometa - Automated Extraction of Genomes from Shotgun Metagenomes
 =====================================================
 projectDir                         : ${workflow.projectDir}
 -----------------------------------------------------
 Data
 -----------------------------------------------------
 metagenome                         : ${params.input}
 interim                            : ${params.interim_dir}
 processed                          : ${params.outdir}
 -----------------------------------------------------
 Parameters
 -----------------------------------------------------
 cpus                               : ${params.cpus}
 length_cutoff                      : ${params.length_cutoff}
 kmer_size                          : ${params.kmer_size}
 norm_method                        : ${params.norm_method}
 pca_dimensions                     : ${params.pca_dimensions}
 embedding_method                   : ${params.embedding_method}
 embedding_dimensions               : ${params.embedding_dimensions}
 clustering_method                  : ${params.clustering_method}
 classification_kmer_pca_dimensions : ${params.classification_kmer_pca_dimensions}
 classification_method              : ${params.classification_method}
 completeness                       : ${params.completeness}
 purity                             : ${params.purity}
 gc_stddev_limit                    : ${params.gc_stddev_limit}
 cov_stddev_limit                   : ${params.cov_stddev_limit}
 kingdom                            : ${params.kingdom}
 -----------------------------------------------------
 Databases
 -----------------------------------------------------
 ncbi_database                      : ${params.ncbi_database}
 -----------------------------------------------------
"""

/*
 * Retrieve the main 'AUTOMETA' worflow
 */
include { AUTOMETA } from './modules/core/autometa.nf'


workflow {
  Channel
    .fromPath(params.input, checkIfExists: true, type: 'file')
    .set{unfiltered_metagenome_ch}

  AUTOMETA(unfiltered_metagenome_ch)
}

/*
 * completion handler
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[autometa] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[autometa] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[autometa] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[autometa] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[autometa]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[autometa]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

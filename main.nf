#!/usr/bin/env nextflow


/*
========================================================================================
    Autometa
========================================================================================
    Autometa's Nextflow Analysis Pipeline
    Github         : https://github.com/KwanLab/Autometa
    Documentation  : https://autometa.readthedocs.io/en/latest/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)


////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////


if (params.use_run_name){
    params.interim_dir_internal = "${params.interim_dir}/autometa_interim_dir/${workflow.runName}/${workflow.sessionId}" // Intermediate results directory
    params.outdir_internal = "${params.outdir}/autometa_outdir/${workflow.runName}/${workflow.sessionId}"           // Final results directory
} else {
    params.interim_dir_internal = "${params.interim_dir}/autometa_interim_dir/${workflow.sessionId}" // Intermediate results directory
    params.outdir_internal = "${params.outdir}/autometa_outdir/${workflow.sessionId}"           // Final results directory
}
println "\n"
println "Intermediate results directory: $params.interim_dir_internal"
println "Binning results directory: $params.outdir_internal"
println "\n"


/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { AUTOMETA } from './workflows/autometa.nf' addParams(single_db_dir: params.single_db_dir)

/*
========================================================================================
    Run Autometa
========================================================================================
*/

workflow {
    AUTOMETA()
}

/*
========================================================================================
    THE END
========================================================================================
*/

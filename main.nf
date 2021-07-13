#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/autometa
========================================================================================
    Github : https://github.com/nf-core/autometa
    Website: https://nf-co.re/autometa
    Slack  : https://nfcore.slack.com/channels/autometa
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { AUTOMETA } from './workflows/autometa'

//
// WORKFLOW: Run main nf-core/autometa analysis pipeline
//
workflow NFCORE_AUTOMETA {
    AUTOMETA ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_AUTOMETA ()
}

/*
========================================================================================
    THE END
========================================================================================
*/

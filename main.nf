#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
****************************************************************************************************************************
                         RNAinn  is a comprehensive  RNA-seq analysis pipeline Nextflow implemented
                         for the analysis of RNA-seq data. It is designed to be highly flexible and
                         can be run on a wide range of computing environments,from a single laptop,
                         to a computing cluster or cloud computing environments. Based on <nf-core>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  Github : https://github.com/jblancoheredia/rnainn
                                                Author : blancoj@mskcc.org
****************************************************************************************************************************
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                               IMPORT WORKFLOW & SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNAINN                                                } from './workflows/rnainn'
include { getGenomeAttribute                                    } from './subworkflows/local/utils_nfcore_rnainn_pipeline'
include { PIPELINE_COMPLETION                                   } from './subworkflows/local/utils_nfcore_rnainn_pipeline'
include { PIPELINE_INITIALISATION                               } from './subworkflows/local/utils_nfcore_rnainn_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                          WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//params.refs_dir             = getGenomeAttribute('refs_dir')

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CTI_RNAINN {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    RNAINN (
        samplesheet
    )

    emit:
    multiqc_report = RNAINN.out.multiqc_report // channel: /path/to/multiqc_report.html

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                      RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.help,
        params.input,
        params.outdir,
        params.version,
        params.monochrome_logs,
        params.validate_params,
        args
    )

    //
    // WORKFLOW: Run main workflow
    //
    CTI_RNAINN (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    channel_multiqc_report = CTI_RNAINN.out.multiqc_report

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.outdir,
        params.hook_url,
        params.email_on_fail,
        channel_multiqc_report,
        params.monochrome_logs,
        params.plaintext_email
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                           THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

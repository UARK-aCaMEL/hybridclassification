#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    aCaMEL/hybridclassification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/aCaMEL/hybridclassification
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HYBRIDCLASSIFICATION  } from './workflows/hybridclassification'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hybridclassification_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hybridclassification_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow ACAMEL_HYBRIDCLASSIFICATION {

    take:
    vcf
    tbi
    popmap
    speciesmap
    site_coords
    geo_data
    geo_data_dir
    combinations

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HYBRIDCLASSIFICATION (
        vcf,
        tbi,
        popmap,
        speciesmap,
        site_coords,
        geo_data,
        geo_data_dir,
        combinations
    )

    emit:
    multiqc_report = HYBRIDCLASSIFICATION.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.popmap,
        params.speciesmap,
        params.site_coords,
        params.geo_data_config,
        params.geo_data_dir,
        params.combinations
    )

    //
    // WORKFLOW: Run main workflow
    //
    ACAMEL_HYBRIDCLASSIFICATION (
        PIPELINE_INITIALISATION.out.vcf,
        PIPELINE_INITIALISATION.out.tbi,
        PIPELINE_INITIALISATION.out.popmap,
        PIPELINE_INITIALISATION.out.speciesmap,
        PIPELINE_INITIALISATION.out.site_coords,
        PIPELINE_INITIALISATION.out.geo_data,
        PIPELINE_INITIALISATION.out.geo_data_dir,
        PIPELINE_INITIALISATION.out.combinations
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        ACAMEL_HYBRIDCLASSIFICATION.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

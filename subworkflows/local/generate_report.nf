include { MULTIQC                } from '../../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_hybridclassification_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

include { PLOT_ADMIXTURE } from '../../modules/local/report/plot_admixture.nf'
include { NH_PLOT_TRACE } from '../../modules/local/report/newhybrids_trace.nf'

workflow GENERATE_REPORT {
    take:
    inds
    pops
    k2_clumpp
    bestK_clumpp
    nh_results
    nh_trace
    nh_map
    popmap
    ch_versions

    main:

    ch_multiqc_files = Channel.empty()


    //Admixture barplots
    admix_inputs = k2_clumpp
                    .join(inds)
                    .join(pops)
    PLOT_ADMIXTURE(
        admix_inputs.map { m, k, i, p -> tuple(m, k) },
        admix_inputs.map { m, k, i, p -> tuple(m, i) },
        admix_inputs.map { m, k, i, p -> tuple(m, p) }
    )
    ch_multiqc_files = PLOT_ADMIXTURE.out.admixture_html
    ch_versions = ch_versions.mix( PLOT_ADMIXTURE.out.versions )

    //Admixture barplots
    NH_PLOT_TRACE(
        nh_trace
    )
    ch_multiqc_files = ch_multiqc_files.join(NH_PLOT_TRACE.out.plot_html)
    ch_versions = ch_versions.mix( NH_PLOT_TRACE.out.versions )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // Prepare MultiQC files
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    //
    // Merge MultiQC files with report plots
    //
    ch_multiqc_files = ch_multiqc_files
        .combine(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .combine(ch_collated_versions)
        .combine(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
    )

    ch_multiqc_files.view()
    MULTIQC (
        ch_multiqc_files,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report  = MULTIQC.out.report.toList()
    versions        = ch_collated_versions
}

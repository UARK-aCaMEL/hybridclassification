/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_hybridclassification_pipeline'

include { ADMIXPIPE              } from '../subworkflows/local/admixpipe.nf'
include { NEWHYBRIDS             } from '../subworkflows/local/newhybrids.nf'
include { SNPIO_FILTER           } from '../modules/local/snpio/pre_filter.nf'
include { SNPIO_SELECT           } from '../modules/local/snpio/select_pops.nf'
include { FIND_CANDIDATES        } from '../modules/local/find_candidates.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HYBRIDCLASSIFICATION {

    take:
    ch_vcf     // [meta, vcf]
    ch_tbi     // [meta, tbi]
    ch_popmap  // [meta, popmap]
    ch_speciesmap // [meta, speciesmap]
    ch_site_coords // [meta, site_coords]
    ch_species_meta // [meta, species_meta]
    ch_combinations // [meta] where meta = [id: "${pop1}_${pop2}", pop1: pop1, pop2: pop2]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Create multi-item channels with combination meta for all inputs
    //
    ch_combo_popmap = ch_combinations
        .combine(ch_popmap)
        .map { combo_meta, orig_meta, popmap ->
            [combo_meta, popmap]
        }

    ch_combo_speciesmap = ch_combinations
        .combine(ch_speciesmap)
        .map { combo_meta, orig_meta, speciesmap ->
            [combo_meta, speciesmap]
        }

    ch_combo_site_coords = ch_combinations
        .combine(ch_site_coords)
        .map { combo_meta, orig_meta, site_coords ->
            [combo_meta, site_coords]
        }

    ch_combo_species_meta = ch_combinations
        .combine(ch_species_meta)
        .map { combo_meta, orig_meta, species_meta ->
            [combo_meta, species_meta]
        }

    //
    // Generate subset VCFs for each test combination (keep original working approach)
    //
    ch_combinations
        .combine(ch_vcf)
        .combine(ch_tbi)
        .combine(ch_speciesmap)
        .map { pops, meta_vcf, vcf, meta_tbi, tbi, meta_speciesmap, speciesmap ->
            [pops, [meta_vcf, vcf], [meta_tbi, tbi], [meta_speciesmap, speciesmap]]
        }
        .set { ch_snpio_input }

    SNPIO_SELECT(
        ch_snpio_input.map { it[1] }, // vcf with meta
        ch_snpio_input.map { it[2] }, // tbi with meta
        ch_snpio_input.map { it[3] }, // speciesmap with meta
        ch_snpio_input.map { it[0] }  // pops combination
    )
    ch_versions = ch_versions.mix(SNPIO_SELECT.out.versions)

    //
    // VCF pre-processing (use combo channels for natural matching)
    //
    SNPIO_FILTER(
        SNPIO_SELECT.out.filtered_vcf
            .join(ch_combo_popmap, by: 0)
            .map { meta, vcf, popmap ->
                [meta, vcf, popmap]
            }
    )
    ch_versions = ch_versions.mix(SNPIO_FILTER.out.versions)

    //
    // Run admixture pipeline - natural meta matching
    //
    ADMIXPIPE(
        SNPIO_FILTER.out.filtered_vcf,
        ch_combo_popmap
    )
    ch_versions = ch_versions.mix(ADMIXPIPE.out.versions)

    //
    // Generate inputs for newhybrids - natural meta matching
    //
    FIND_CANDIDATES(
        ADMIXPIPE.out.k2_clumpp,
        ADMIXPIPE.out.inds,
        ch_combo_popmap,
        ch_combo_speciesmap
    )
    ch_versions = ch_versions.mix( FIND_CANDIDATES.out.versions )

    //
    // NewHybrids subworkflow - natural meta matching
    //
    NEWHYBRIDS(
        SNPIO_FILTER.out.filtered_vcf,
        SNPIO_FILTER.out.filtered_tbi,
        FIND_CANDIDATES.out.popmap
    )
    ch_versions = ch_versions.mix( NEWHYBRIDS.out.versions )

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
    // MODULE: MultiQC
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

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/